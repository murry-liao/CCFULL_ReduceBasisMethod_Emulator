"""Run a trained reduced model at parameter points, and compare it with exact CCFULL.

    python -m emulator.run_rom                            # reads emulator/emulator_rom.in
    python -m emulator.run_rom --config <file> --check    # only prints how each line was read

The lines of emulator_rom.in are explained at its end.  Every run writes a new directory in
the output folder with cross_sections.csv (one row per point and energy), point_timings.csv
(the ROM and the exact time of every point), and run_summary.json.
"""

import argparse
import csv
from datetime import datetime
import json
import math
import os
import time

# The ROM is many small dense solves: one BLAS thread is the fastest, and it keeps the
# timings comparable.  BLAS reads these variables when numpy is imported, so they come first.
for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(variable, "1")

import numpy as np

from emulator import rom
from emulator.ccfull_run import (PROJECT, data_lines, exact_cross_sections, numbers, project_path,
                                 read_experimental)


def write_csv(path, rows):
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(key for row in rows for key in row)))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Run a trained reduced model of CCFULL.")
    parser.add_argument("--config", default=PROJECT / "emulator" / "emulator_rom.in")
    parser.add_argument("--check", action="store_true", help="print how each line was read, and stop")
    args = parser.parse_args()

    # 1. Read emulator_rom.in (its lines are explained at its end).
    lines = data_lines(args.config)
    model_name, training_root = lines[0], lines[1]
    model = rom.load_model(project_path(training_root) / model_name, int(numbers(lines[2])[0]))
    mode = lines[3]
    box = model["box"] if lines[7] == "training" else np.reshape(numbers(lines[8]) + numbers(lines[9]), (6, 2))
    if mode == "single":
        points = np.array([numbers(lines[4])])
    elif mode == "latin_hypercube":
        count, seed = (int(value) for value in numbers(lines[5]))
        points = rom.latin_hypercube(box, count, seed)
    else:
        points = np.loadtxt(project_path(lines[6]), delimiter=",", skiprows=1, ndmin=2)
    measured_energy, measured_sigma, measured_lower, measured_upper = read_experimental(
        project_path(model["experimental_file"]))
    if lines[10] == "experimental":
        energies = measured_energy
    elif lines[10] == "grid":
        first, last, step = numbers(lines[11])
        energies = np.array([round(first + k * step, 10) for k in range(int(math.floor((last - first) / step + 1e-9)) + 1)])
    else:
        energies = np.array(numbers(lines[12]))
    show_data, check_exact, report_speedup = (value != 0 for value in numbers(lines[13]))
    output_folder = project_path(lines[14])
    save_cross_sections, save_summary, save_timings = (value != 0 for value in numbers(lines[15]))
    system = model["system"]
    reaction = f"A={system['apro']}, Z={system['zpro']} + A={system['atar']}, Z={system['ztar']}"

    # 2. Say how each line was read.
    print(f"ROM input {args.config}")
    print(f"  lines 1-2    {model['type']} model {model_name} ({reaction})")
    print(f"  line 3       basis rank {model['rank']} of the trained {model['trained_rank']}"
          + (" per (E, J)" if model["type"] == "local" else ""))
    print(f"  lines 4-7    {mode}, {len(points)} point(s)")
    print(f"  lines 8-10   {lines[7]} box: " + ", ".join(
        f"{name} [{low:g}, {high:g}]" for name, (low, high) in zip(rom.PARAMETER_NAMES, box)))
    print(f"  lines 11-13  {lines[10]} energies: {len(energies)} from {min(energies):g} to {max(energies):g} MeV")
    print(f"  line 14      data {show_data}, exact CCFULL check {check_exact}, speedup {report_speedup}")
    print(f"  lines 15-16  output {output_folder}; files: cross_sections.csv {save_cross_sections}, "
          f"run_summary.json {save_summary}, point_timings.csv {save_timings}")
    if args.check:
        return

    # 3. Prepare the model at these energies: the local blocks, or the outer boundaries of
    #    the global basis.  The timed points do not include it.
    start = time.perf_counter()
    if model["type"] == "local":
        rom.load_blocks(model, energies)
    else:
        rom.prepare_global(model, energies)
    preparation = time.perf_counter() - start
    print(f"[setup] {model_name} at rank {model['rank']}, prepared in {preparation:.1f} s", flush=True)

    # 4. The reduced model at every point.
    print(f"[ROM] {len(points)} point(s) at {len(energies)} energies", flush=True)
    sigma_rom, rom_seconds = [], []
    start_all = time.perf_counter()
    for k, theta in enumerate(points, start=1):
        start = time.perf_counter()
        sigma_rom.append(rom.cross_sections(model, theta, energies))
        rom_seconds.append(time.perf_counter() - start)
        if k % max(1, len(points) // 20) == 0 or k == len(points):
            print(f"[ROM] {k}/{len(points)} points ({round(100 * k / len(points)):3d}%)", flush=True)
    rom_wall = time.perf_counter() - start_all

    # 5. Exact CCFULL (the clean program) at the same points, with the r_min of the model.
    sigma_exact, exact_seconds, exact_wall = None, None, None
    if check_exact:
        print(f"[exact CCFULL] {len(points)} point(s)", flush=True)
        sigma_exact, exact_seconds = [], []
        start_all = time.perf_counter()
        for k, theta in enumerate(points, start=1):
            start = time.perf_counter()
            sigma_exact.append(exact_cross_sections(model["records"], theta, energies, model["radius"][0]))
            exact_seconds.append(time.perf_counter() - start)
            if k % max(1, len(points) // 20) == 0 or k == len(points):
                print(f"[exact CCFULL] {k}/{len(points)} points ({round(100 * k / len(points)):3d}%)", flush=True)
        exact_wall = time.perf_counter() - start_all

    # 6. Accuracy: |sigma_ROM / sigma_exact - 1| over every point and energy.
    accuracy = None
    if check_exact:
        relative = np.abs(np.array(sigma_rom) - np.array(sigma_exact)) / np.maximum(np.abs(sigma_exact), 1e-14)
        worst = np.unravel_index(int(np.argmax(relative)), relative.shape)
        accuracy = {"values": int(relative.size), "median": float(np.median(relative)),
                    "p95": float(np.percentile(relative, 95)), "p99": float(np.percentile(relative, 99)),
                    "max": float(relative.max()),
                    "max_at": {"point_index": int(worst[0]) + 1, "energy_cm_mev": float(energies[worst[1]])},
                    "fraction_above_1_percent": float(np.mean(relative > 0.01))}

    # 7. The output directory.
    created = datetime.now().astimezone()
    suffix = "" if model["rank"] == model["trained_rank"] else f"_r{model['rank']}"
    folder = output_folder / f"{created.strftime('%Y%m%d_%H%M%S_%f')}_{model_name}{suffix}"
    folder.mkdir(parents=True)
    measured = {round(float(e), 8): (s, lo, up) for e, s, lo, up in
                zip(measured_energy, measured_sigma, measured_lower, measured_upper)}

    rows = []
    for k, theta in enumerate(points):
        for i, energy in enumerate(energies):
            row = {"point_index": k + 1, **dict(zip(rom.PARAMETER_NAMES, map(float, theta))),
                   "energy_cm_mev": float(energy), "emulator_cross_section_mb": float(sigma_rom[k][i])}
            if check_exact:
                row["exact_numerov_cross_section_mb"] = float(sigma_exact[k][i])
                row["rom_exact_relative_error"] = float(relative[k, i])
            if show_data and round(float(energy), 8) in measured:
                s, lo, up = measured[round(float(energy), 8)]
                row.update({"experimental_cross_section_mb": float(s), "experimental_uncertainty_lower_mb": float(lo),
                            "experimental_uncertainty_upper_mb": float(up)})
            rows.append(row)
    if save_cross_sections:
        write_csv(folder / "cross_sections.csv", rows)

    timing_rows = []
    for k, theta in enumerate(points):
        row = {"point_index": k + 1, **dict(zip(rom.PARAMETER_NAMES, map(float, theta))), "rom_seconds": rom_seconds[k]}
        if check_exact:
            row.update({"exact_seconds": exact_seconds[k], "speedup": exact_seconds[k] / rom_seconds[k]})
        timing_rows.append(row)
    if save_timings:
        write_csv(folder / "point_timings.csv", timing_rows)

    if save_summary:
        outside = lambda theta: {name: {"value": float(value), "trained_range": [float(low), float(high)]}
                                 for name, value, (low, high) in zip(rom.PARAMETER_NAMES, theta, model["box"])
                                 if value < low or value > high}
        speedups = np.array([row["speedup"] for row in timing_rows]) if check_exact else None
        summary = {
            "format": "ccfull_rom_run_v3", "created_at": created.isoformat(), "configuration": str(args.config),
            "model_id": model_name, "model_folder": str(model["folder"]), "reaction_system": reaction,
            "model_type": model["type"], "basis_rank": model["rank"], "trained_rank": model["trained_rank"],
            "prediction_mode": mode,
            "parameter_domain": {name: [float(low), float(high)] for name, (low, high) in zip(rom.PARAMETER_NAMES, box)},
            "training_parameter_domain": {name: [float(low), float(high)]
                                          for name, (low, high) in zip(rom.PARAMETER_NAMES, model["box"])},
            "training_energy_range_mev": list(model["energy_range"]),
            "energies": {"source": lines[10], "values_mev": [float(e) for e in energies]},
            "exact_model": "physical_model/ccfull_clean, one run per energy, r_min fixed at the model grid",
            "points_outside_training_domain": sum(bool(outside(theta)) for theta in points),
            "points": [{"point_index": k + 1, "parameters": dict(zip(rom.PARAMETER_NAMES, map(float, theta))),
                        "outside_training_domain": outside(theta)} for k, theta in enumerate(points)],
            "accuracy": accuracy,
            "timings_seconds": {"preparation": preparation, "rom_end_to_end": rom_wall,
                                "rom_per_point_mean": float(np.mean(rom_seconds)), "exact_numerov": exact_wall,
                                "exact_per_point_mean": float(np.mean(exact_seconds)) if check_exact else None},
            "speedup": None if not (check_exact and report_speedup) else {
                "exact_over_rom_wall_time": exact_wall / rom_wall,
                "per_point": {"mean": float(speedups.mean()), "median": float(np.median(speedups)),
                              "min": float(speedups.min()), "max": float(speedups.max())}},
        }
        (folder / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # 8. Report.
    print(f"[complete] {len(points)} point(s) saved to {folder}")
    print(f"ROM: {rom_wall:.1f} s wall, {np.mean(rom_seconds):.2f} s per point")
    if check_exact:
        print(f"Exact CCFULL: {exact_wall:.1f} s wall, {np.mean(exact_seconds):.2f} s per point")
        print(f"Accuracy, |sigma_ROM/sigma_exact - 1| over {accuracy['values']} values: median {accuracy['median']:.1e}, "
              f"95th percentile {accuracy['p95']:.1e}, 99th percentile {accuracy['p99']:.1e}, max {accuracy['max']:.1e} "
              f"(point {accuracy['max_at']['point_index']}, E={accuracy['max_at']['energy_cm_mev']:.5g} MeV)")
        if report_speedup:
            print(f"Exact CCFULL / ROM wall-time ratio: {exact_wall / rom_wall:.3f}x")


if __name__ == "__main__":
    main()
