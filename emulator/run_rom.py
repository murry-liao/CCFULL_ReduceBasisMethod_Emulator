#!/usr/bin/env python3
"""Online CCFULL-ROM predictions, controlled by emulator_rom.in."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import dataclass
from datetime import datetime
import json
import math
import os
from pathlib import Path
import time

# The online ROM is many dense solves of a few milliseconds each.  Threaded
# BLAS is slower than one thread for them, and the worker processes would
# compete for its threads.  BLAS reads these variables when NumPy is imported,
# so they are set first; a value chosen by the user is kept.
for _variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                  "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_variable, "1")

import numpy as np

from emulator.rom import (PARAMETER_NAMES, PROJECT_ROOT, ROM, latin_hypercube, project_path,
                          project_relative, read_manifest)
from physical_model.inputs import data_lines, load_dataset, numbers, parse_input, reaction_label


DEFAULT_INPUT = PROJECT_ROOT / "emulator" / "emulator_rom.in"


@dataclass(frozen=True)
class RomInput:
    """emulator_rom.in, completed from the manifest of its model; its lines are explained at the end."""

    source: Path
    model_id: str                        # line 1
    manifest: Path                       # line 2: training_root/operator_eim/<model_id>/manifest.json
    workers: int                         # line 3
    mode: str                            # line 4: "single", "latin_hypercube", or "table"
    points: np.ndarray                   # lines 5-7: (point, parameter), PARAMETER_NAMES order
    box: str                             # line 8: "training" or "custom"
    bounds: np.ndarray                   # the Latin-hypercube box; lines 9-10 if custom
    energy_source: str                   # line 11: "experimental", "grid", or "selected"
    energies: np.ndarray                 # lines 12-13, MeV
    include_experiment: bool             # line 14
    exact_numerov: bool
    report_speedup: bool
    output_dir: Path                     # line 15
    save_cross_sections_csv: bool        # line 16
    save_summary_json: bool
    save_point_timings_csv: bool
    reaction: str                        # the rest is recorded in the manifest of the model
    rank: int
    experimental_file: Path
    trained_bounds: np.ndarray
    energy_range: tuple[float, float]


def read_input(path: str | Path) -> RomInput:
    """Read emulator_rom.in line by line; the trained box, energies, and data come from the model."""
    source = Path(path).resolve()
    (model_id, root, workers, mode, single, lhs, table, box, deformations, potential,
     energy_source, grid, selected, switches, output_dir, saves) = data_lines(source)[:16]
    manifest = project_path(root) / "operator_eim" / model_id / "manifest.json"
    trained = read_manifest(manifest)
    settings = trained["settings"]
    trained_bounds = np.asarray(settings["parameter_bounds"], dtype=float)
    experimental_file = project_path(settings["experimental_file"])
    bounds = (trained_bounds if box == "training"
              else np.reshape(numbers(deformations) + numbers(potential), (6, 2)))
    if mode == "single":
        points = np.asarray([numbers(single)])
    elif mode == "latin_hypercube":
        count, seed = map(int, numbers(lhs))
        points = latin_hypercube(bounds, count, seed)
    else:
        points = np.loadtxt(project_path(table), delimiter=",", skiprows=1, ndmin=2)
    if energy_source == "experimental":
        energies = [row["energy_cm_mev"] for row in load_dataset(experimental_file)]
    elif energy_source == "grid":
        first, last, step = numbers(grid)
        count = int(math.floor((last - first) / step + 1.0e-9)) + 1
        energies = [round(first + index * step, 10) for index in range(count)]
    else:
        energies = numbers(selected)
    include_experiment, exact_numerov, report_speedup = (value != 0 for value in numbers(switches))
    save_csv, save_json, save_timings = (value != 0 for value in numbers(saves))
    return RomInput(
        source=source, model_id=model_id, manifest=manifest, workers=int(numbers(workers)[0]),
        mode=mode, points=points, box=box, bounds=bounds,
        energy_source=energy_source, energies=np.asarray(energies, dtype=float),
        include_experiment=include_experiment, exact_numerov=exact_numerov,
        report_speedup=report_speedup and exact_numerov,      # the speedup needs the exact check
        output_dir=project_path(output_dir), save_cross_sections_csv=save_csv,
        save_summary_json=save_json, save_point_timings_csv=save_timings,
        reaction=reaction_label(parse_input(project_path(settings["input_file"]))),
        rank=int(trained["rank"]), experimental_file=experimental_file,
        trained_bounds=trained_bounds, energy_range=tuple(settings["energy_range_mev"]))


_MODEL: ROM | None = None      # the model of this process, with its own CCFULL library


def _load_model(manifest: str) -> None:
    global _MODEL
    _MODEL = ROM.load(manifest)


def _rom_point(task) -> tuple[np.ndarray, float]:
    point, energies = task
    started = time.perf_counter()
    values = _MODEL.cross_sections(point, energies)
    return values, time.perf_counter() - started


def _exact_point(task) -> tuple[np.ndarray, float]:
    point, energies = task
    started = time.perf_counter()
    values = _MODEL.exact_cross_sections(point, energies)
    return values, time.perf_counter() - started


def _evaluate(function, tasks: list, pool, stage: str) -> tuple[np.ndarray, np.ndarray, float]:
    """Run ``function`` for every point, in the worker pool or in this process."""
    started = time.perf_counter()
    results = map(function, tasks) if pool is None else pool.map(function, tasks)
    values, seconds = [], []
    step = max(1, len(tasks) // 20)
    for done, (value, elapsed) in enumerate(results, start=1):
        values.append(value)
        seconds.append(elapsed)
        if done % step == 0 or done == len(tasks):
            print(f"[{stage}] {done}/{len(tasks)} points ({round(100.0 * done / len(tasks)):3d}%)", flush=True)
    return np.asarray(values), np.asarray(seconds), time.perf_counter() - started


def _domain(bounds: np.ndarray) -> dict:
    return {name: [float(lower), float(upper)] for name, (lower, upper) in zip(PARAMETER_NAMES, bounds)}


def _outside(point: np.ndarray, bounds: np.ndarray) -> dict:
    """The parameters of ``point`` outside the trained box, with their trained range."""
    outside = {}
    for name, value, (lower, upper) in zip(PARAMETER_NAMES, map(float, point), bounds):
        if value < lower or value > upper:
            outside[name] = {"value": value, "trained_range": [float(lower), float(upper)]}
    return outside


def _accuracy(rom: np.ndarray, exact: np.ndarray, energies: np.ndarray) -> dict:
    """|sigma_ROM - sigma_exact| / |sigma_exact| over every point and energy."""
    relative = np.abs(rom - exact) / np.maximum(np.abs(exact), 1e-14)
    worst = np.unravel_index(int(np.argmax(relative)), relative.shape)
    return {
        "values": int(relative.size),
        "median": float(np.median(relative)),
        "p95": float(np.percentile(relative, 95)),
        "p99": float(np.percentile(relative, 99)),
        "max": float(relative.max()),
        "max_at": {"point_index": int(worst[0]) + 1, "energy_cm_mev": float(energies[worst[1]])},
        "fraction_above_1_percent": float(np.mean(relative > 0.01)),
    }


def _rows(energies: np.ndarray, rom: np.ndarray, exact: np.ndarray | None,
          experiment: list[dict]) -> list[dict]:
    """The cross_sections.csv rows of one point, one per energy."""
    measured = {round(float(row["energy_cm_mev"]), 8): row for row in experiment}
    rows = []
    for index, energy in enumerate(energies):
        row = {"energy_cm_mev": float(energy), "emulator_cross_section_mb": float(rom[index])}
        if exact is not None:
            row["exact_numerov_cross_section_mb"] = float(exact[index])
            row["rom_exact_relative_error"] = float(
                abs(rom[index] - exact[index]) / max(abs(exact[index]), 1e-14))
        data = measured.get(round(float(energy), 8))
        if data:
            row.update({
                "experimental_cross_section_mb": float(data["fusion_cross_section_mb"]),
                "experimental_uncertainty_lower_mb": float(data["uncertainty_lower_mb"]),
                "experimental_uncertainty_upper_mb": float(data["uncertainty_upper_mb"]),
            })
        rows.append(row)
    return rows


def _write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _point_timings(points: np.ndarray, rom_seconds: np.ndarray,
                   exact_seconds: np.ndarray | None) -> list[dict]:
    """ROM and exact time of every point, each measured in its worker process over all energies.

    The first point of every worker also builds the outer Coulomb maps of the
    energies, which its later points reuse.
    """
    rows = []
    for index, point in enumerate(points):
        row = {"point_index": index + 1, **dict(zip(PARAMETER_NAMES, map(float, point))),
               "rom_seconds": float(rom_seconds[index])}
        if exact_seconds is not None:
            row["exact_seconds"] = float(exact_seconds[index])
            row["speedup"] = row["exact_seconds"] / max(row["rom_seconds"], 1.0e-300)
        rows.append(row)
    return rows


def _save(config: RomInput, rom: np.ndarray, exact: np.ndarray | None, experiment: list[dict],
          accuracy: dict | None, timings: dict, point_timings: list[dict]) -> Path:
    """Write cross_sections.csv, point_timings.csv, and run_summary.json into a new run directory."""
    created = datetime.now().astimezone()
    destination = config.output_dir / f"{created.strftime('%Y%m%d_%H%M%S_%f')}_{config.model_id}"
    destination.mkdir(parents=True, exist_ok=False)
    rows, points = [], []
    for index, point in enumerate(config.points):
        parameters = dict(zip(PARAMETER_NAMES, map(float, point)))
        rows.extend({"point_index": index + 1, **parameters, **row} for row in _rows(
            config.energies, rom[index], None if exact is None else exact[index], experiment))
        points.append({"point_index": index + 1, "parameters": parameters,
                       "outside_training_domain": _outside(point, config.trained_bounds)})
    if config.save_cross_sections_csv:
        _write_csv(destination / "cross_sections.csv", rows)
    if config.save_point_timings_csv:
        _write_csv(destination / "point_timings.csv", point_timings)
    speedup = None
    if config.report_speedup:
        per_point = np.asarray([row["speedup"] for row in point_timings])
        speedup = {
            "exact_over_rom_wall_time": timings["exact_numerov"] / max(timings["rom_end_to_end"], 1.0e-12),
            "per_point": {"mean": float(per_point.mean()), "median": float(np.median(per_point)),
                          "min": float(per_point.min()), "max": float(per_point.max())},
        }
    if config.save_summary_json:
        summary = {
            "format": "ccfull_rom_run_v2",
            "created_at": created.isoformat(),
            "configuration": str(config.source),
            "model_id": config.model_id,
            "model_manifest": project_relative(config.manifest),
            "reaction_system": config.reaction,
            "basis_rank": config.rank,
            "prediction_mode": config.mode,
            "parameter_domain": _domain(config.bounds),
            "training_parameter_domain": _domain(config.trained_bounds),
            "training_energy_range_mev": list(config.energy_range),
            "energies": {"source": config.energy_source,
                         "values_mev": [float(value) for value in config.energies]},
            "workers": config.workers,
            "points_outside_training_domain": sum(bool(item["outside_training_domain"]) for item in points),
            "points": points,
            "accuracy": accuracy,
            "timings_seconds": timings,
            "speedup": speedup,
        }
        (destination / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return destination


def _report(config: RomInput) -> None:
    """What each line of the input was read as."""
    names = ("beta2", "beta4", "beta6", "V0", "r0", "a0")
    yes = {True: "yes", False: "no"}
    points, energies = config.points, config.energies
    single = (": " + ", ".join(f"{name}={value:g}" for name, value in zip(names, points[0]))
              if config.mode == "single" else "")
    print(f"ROM input {config.source}")
    print(f"  lines 1-2    model {config.model_id} ({config.reaction}), basis rank {config.rank}: "
          f"{project_relative(config.manifest)}")
    print(f"  line 3       {config.workers} worker process(es)")
    print(f"  lines 4-7    {config.mode}, {len(points)} point(s){single}")
    print(f"  lines 8-10   {config.box} box: " + ", ".join(
        f"{name} [{lower:g}, {upper:g}]" for name, (lower, upper) in zip(names, config.bounds)))
    print(f"  lines 11-13  {config.energy_source} energies: {energies.size} between "
          f"{energies.min():g} and {energies.max():g} MeV")
    print(f"  line 14      experiment {yes[config.include_experiment]}, exact Numerov check "
          f"{yes[config.exact_numerov]}, speedup {yes[config.report_speedup]}")
    print(f"  line 15      output directory {project_relative(config.output_dir)}")
    print(f"  line 16      write cross_sections.csv {yes[config.save_cross_sections_csv]}, "
          f"run_summary.json {yes[config.save_summary_json]}, "
          f"point_timings.csv {yes[config.save_point_timings_csv]}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--check", action="store_true",
                        help="print what each line of the input was read as, without running")
    args = parser.parse_args()
    config = read_input(args.config)
    if args.check:
        _report(config)
        return
    points, energies = config.points, config.energies
    print(f"[setup] Loading trained model: {config.manifest}", flush=True)
    started = time.perf_counter()
    if config.workers == 1:
        _load_model(str(config.manifest))
    loading = time.perf_counter() - started
    outside = sum(bool(_outside(point, config.trained_bounds)) for point in points)
    if outside:
        print(f"Extrapolation: {outside}/{len(points)} point(s) are outside the trained parameter domain.")
    tasks = [(point, energies) for point in points]
    pool = (None if config.workers == 1 else ProcessPoolExecutor(
        max_workers=config.workers, initializer=_load_model, initargs=(str(config.manifest),)))
    try:
        print(f"[ROM] Evaluating {len(points)} parameter point(s) at {energies.size} energies, "
              f"basis rank {config.rank}, {config.workers} worker process(es).", flush=True)
        rom, rom_seconds, rom_wall = _evaluate(_rom_point, tasks, pool, "ROM")
        exact = exact_seconds = exact_wall = None
        if config.exact_numerov:
            print(f"[exact Numerov] Checking {len(points)} parameter point(s).", flush=True)
            exact, exact_seconds, exact_wall = _evaluate(_exact_point, tasks, pool, "exact Numerov")
    finally:
        if pool is not None:
            pool.shutdown()
    accuracy = None if exact is None else _accuracy(rom, exact, energies)
    experiment = load_dataset(config.experimental_file) if config.include_experiment else []
    timings = {
        "model_loading": loading,
        "rom_end_to_end": rom_wall,
        "rom_per_point_mean": float(rom_seconds.mean()),
        "exact_numerov": exact_wall,
        "exact_per_point_mean": None if exact_seconds is None else float(exact_seconds.mean()),
    }
    print("[output] Writing the run directory.", flush=True)
    destination = _save(config, rom, exact, experiment, accuracy, timings,
                        _point_timings(points, rom_seconds, exact_seconds))
    print(f"[complete] Saved {len(points)} prediction(s) to {destination}", flush=True)
    print(f"ROM: {rom_wall:.1f} s wall, {rom_seconds.mean():.2f} s per point in one process")
    if exact is not None:
        print(f"Exact Numerov: {exact_wall:.1f} s wall, {exact_seconds.mean():.2f} s per point in one process")
        print(f"Accuracy, |sigma_ROM/sigma_exact - 1| over {accuracy['values']} values: "
              f"median {accuracy['median']:.1e}, 95th percentile {accuracy['p95']:.1e}, "
              f"99th percentile {accuracy['p99']:.1e}, max {accuracy['max']:.1e} "
              f"(point {accuracy['max_at']['point_index']}, "
              f"E={accuracy['max_at']['energy_cm_mev']:.5g} MeV)")
        if config.report_speedup:
            print(f"Exact Numerov / ROM wall-time ratio: {exact_wall / max(rom_wall, 1.0e-12):.3f}x")


if __name__ == "__main__":
    main()
