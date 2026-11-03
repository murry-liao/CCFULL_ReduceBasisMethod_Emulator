"""The input files of CCFULL-ROM and the manifests of trained models.

``emulator_train.in`` is the CCFULL input of the reaction (nine records)
followed by the TOML training controls.  ``emulator_rom.in`` is TOML; it names
a trained model, whose manifest gives the parameter box, the energy range, and
the files the model was trained with.  Relative paths are read from the
project root.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import re
import tomllib

import numpy as np
from scipy.stats import qmc

from physical_model import ccfull
from physical_model.inputs import data_nuclei, energy_values, load_dataset, parse_input, records_sha256


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TRAINING_INPUT = PROJECT_ROOT / "emulator" / "configs" / "emulator_train.in"
DEFAULT_ROM_INPUT = PROJECT_ROOT / "emulator" / "configs" / "emulator_rom.in"
PARAMETER_NAMES = ("beta2", "beta4", "beta6", "v0", "r0", "a0")
MODEL_FORMAT = "global_projection_cc_eim_v1"
ENERGY_TOLERANCE = 1.0e-9    # MeV


def project_relative(path: str | Path) -> str:
    """``path`` as a manifest records it: relative to the project root, or absolute outside it."""
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def project_path(recorded: str | Path) -> Path:
    """A path recorded by ``project_relative``."""
    path = Path(recorded)
    return path if path.is_absolute() else PROJECT_ROOT / path


def latin_hypercube(bounds, count: int, seed: int) -> np.ndarray:
    """``count`` Latin-hypercube points inside ``bounds`` (n, 2); equal bounds fix a coordinate."""
    limits = np.asarray(bounds, dtype=float)
    points = np.repeat(limits[:, 0][None, :], count, axis=0)
    active = limits[:, 1] > limits[:, 0]
    if np.any(active):
        unit = qmc.LatinHypercube(int(np.sum(active)), seed=seed).random(count)
        points[:, active] = qmc.scale(unit, limits[active, 0], limits[active, 1])
    return points


def read_manifest(path: str | Path) -> dict:
    """A trained model's manifest; refused if ccfull.f90 or the nine CCFULL records changed since."""
    path = Path(path)
    manifest = json.loads(path.read_text())
    if manifest.get("format") != MODEL_FORMAT:
        raise ValueError(f"unsupported model format {manifest.get('format')!r} in {path}")
    settings = manifest["settings"]
    if settings.get("solver_sha256") != ccfull.source_sha256():
        raise ValueError(f"{path} was trained with another ccfull.f90. Move its model directory "
                         "out of training_data and train again.")
    input_file = project_path(settings["input_file"])
    if not input_file.is_file():
        raise FileNotFoundError(f"{path} needs its training input {input_file}")
    if settings.get("physical_input_sha256") != records_sha256(input_file):
        raise ValueError(f"{path} was trained with other CCFULL records than {input_file} now holds. "
                         "Restore the records, or move the model directory out of training_data "
                         "and train again.")
    return manifest


# Readers of single values; ``section`` names the TOML table in error messages.

def _existing(path: str | Path) -> Path:
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"input file does not exist: {source}")
    return source


def _table(payload: dict, name: str) -> dict:
    value = payload.get(name)
    if not isinstance(value, dict):
        raise ValueError(f"missing [{name}] table")
    return value


def _path(table: dict, key: str, section: str) -> Path:
    value = table.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{section}.{key} must be a path")
    path = Path(value).expanduser()
    return (path if path.is_absolute() else PROJECT_ROOT / path).resolve()


def _integer(table: dict, key: str, section: str, minimum: int = 1, default: int | None = None) -> int:
    value = table.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{section}.{key} must be an integer >= {minimum}")
    return value


def _positive(table: dict, key: str, section: str) -> float:
    value = table.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
        raise ValueError(f"{section}.{key} must be a positive number")
    return float(value)


def _boolean(table: dict, key: str, section: str, default: bool | None = None) -> bool:
    value = table.get(key, default)
    if not isinstance(value, bool):
        raise ValueError(f"{section}.{key} must be true or false")
    return value


def _model_id(value: object) -> str:
    model_id = str(value or "").strip()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", model_id):
        raise ValueError("model.model_id may contain only letters, digits, periods, underscores, "
                         "and hyphens")
    return model_id


def _bounds(table: dict, section: str) -> dict[str, tuple[float, float]]:
    """[lower, upper] of every parameter; equal endpoints fix a parameter."""
    bounds = {}
    for name in PARAMETER_NAMES:
        value = table.get(name)
        if not isinstance(value, list) or len(value) != 2:
            raise ValueError(f"{section}.{name} must be [lower, upper]")
        lower, upper = map(float, value)
        if not (math.isfinite(lower) and math.isfinite(upper) and lower <= upper):
            raise ValueError(f"{section}.{name} must be finite with lower <= upper")
        bounds[name] = (lower, upper)
    return bounds


def _reaction_label(physical: dict) -> str:
    return (f"A={physical['apro']}, Z={physical['zpro']} + "
            f"A={physical['atar']}, Z={physical['ztar']}")


@dataclass(frozen=True)
class TrainConfig:
    """emulator_train.in."""

    source: Path                                 # this file, also the CCFULL input of the model
    model_id: str
    reaction: str                                # for example "A=16, Z=8 + A=154, Z=62"
    experimental_file: Path
    training_root: Path
    bounds: dict[str, tuple[float, float]]       # in PARAMETER_NAMES order
    energy_source: str                           # "experimental", "ccfull", or "custom"
    energy_range: tuple[float, float]            # MeV
    training_points: int                         # snapshots, one exact CCFULL solve each
    rank: int                                    # POD modes kept for the basis
    seed: int
    workers: int                                 # CCFULL processes at a time
    eim_points: int                              # parameter points of the CC-EIM training
    eim_seed: int
    eim_tolerance: float
    eim_max_rank: int                            # per multipole

    @property
    def wave_pod_dir(self) -> Path:
        return self.training_root / "wave_pod" / self.model_id

    @property
    def cc_eim_dir(self) -> Path:
        return self.training_root / "operator_eim" / self.model_id


def _energy_range(value: object, physical: dict, experimental_file: Path) -> tuple[str, tuple[float, float]]:
    """The trained energy interval from energy.range: "experimental", "ccfull", or [lower, upper]."""
    if value == "experimental":
        # The measured energies, rounded outward to whole MeV.
        energies = [row["energy_cm_mev"] for row in load_dataset(experimental_file)]
        source, lower, upper = value, math.floor(min(energies)), math.ceil(max(energies))
    elif value == "ccfull":
        energies = energy_values(physical)
        source, lower, upper = value, min(energies), max(energies)
    elif isinstance(value, list) and len(value) == 2:
        source, (lower, upper) = "custom", value
    else:
        raise ValueError('energy.range must be "experimental", "ccfull", or [lower, upper] in MeV')
    lower, upper = float(lower), float(upper)
    if not (math.isfinite(lower) and math.isfinite(upper) and 0.0 < lower < upper):
        raise ValueError(f"the training energy range [{lower}, {upper}] MeV needs 0 < lower < upper")
    return source, (lower, upper)


def load_train_config(path: str | Path = DEFAULT_TRAINING_INPUT) -> TrainConfig:
    """Read and check emulator_train.in."""
    source = _existing(path)
    lines = source.read_text().splitlines()
    start = next((index for index, line in enumerate(lines) if line.lstrip().startswith("[")), None)
    if start is None:
        raise ValueError(f"{source} has no training tables after the nine CCFULL records")
    payload = tomllib.loads("\n".join(lines[start:]))
    physical = parse_input(source)
    if int(physical["target_mode"]) != 1:
        raise ValueError("CCFULL-ROM needs a rotational target (CCFULL target mode 1)")

    model, data = _table(payload, "model"), _table(payload, "data")
    requested = str(model.get("model_id", "auto")).strip()
    model_id = _model_id(
        f"a{physical['apro']}_z{physical['zpro']}__a{physical['atar']}_z{physical['ztar']}"
        if requested.lower() == "auto" else requested)
    experimental_file = _path(data, "experimental_file", "data")
    if not experimental_file.is_file():
        raise FileNotFoundError(f"data.experimental_file does not exist: {experimental_file}")
    reaction = ((physical["zpro"], physical["apro"]), (physical["ztar"], physical["atar"]))
    if data_nuclei(experimental_file) not in (reaction, reaction[::-1]):
        raise ValueError(f"{experimental_file.name} is not the reaction of the CCFULL records")

    bounds = _bounds(_table(payload, "parameter_bounds"), "parameter_bounds")
    target = physical["target1"]
    nominal = (target[1], target[2], target[3], physical["v0"], physical["r0"], physical["a0"])
    for name, value in zip(PARAMETER_NAMES, map(float, nominal)):
        lower, upper = bounds[name]
        if not lower <= value <= upper:
            raise ValueError(f"the CCFULL value {name} = {value} is outside its training range "
                             f"[{lower}, {upper}]")
    energy_source, energy_range = _energy_range(
        _table(payload, "energy").get("range", "experimental"), physical, experimental_file)

    wave, eim = _table(payload, "wave_pod"), _table(payload, "cc_eim")
    training_points = _integer(wave, "training_points", "wave_pod", 2)
    rank = _integer(wave, "rank", "wave_pod")
    if rank > training_points:
        raise ValueError("wave_pod.rank (basis vectors) cannot exceed "
                         "wave_pod.training_points (snapshots)")
    return TrainConfig(
        source=source, model_id=model_id, reaction=_reaction_label(physical),
        experimental_file=experimental_file, training_root=_path(data, "training_root", "data"),
        bounds=bounds, energy_source=energy_source, energy_range=energy_range,
        training_points=training_points, rank=rank,
        seed=_integer(wave, "seed", "wave_pod", 0),
        workers=_integer(wave, "workers", "wave_pod", default=1),
        eim_points=_integer(eim, "interaction_points", "cc_eim", 2),
        eim_seed=_integer(eim, "seed", "cc_eim", 0),
        eim_tolerance=_positive(eim, "tolerance", "cc_eim"),
        eim_max_rank=_integer(eim, "max_rank_per_multipole", "cc_eim"))


@dataclass(frozen=True)
class RomConfig:
    """emulator_rom.in, completed from the manifest of its model."""

    source: Path
    model_id: str
    manifest: Path                               # training_root/operator_eim/<model_id>/manifest.json
    workers: int                                 # worker processes
    reaction: str
    rank: int                                    # basis rank of the model
    experimental_file: Path                      # the model's
    trained_bounds: dict[str, tuple[float, float]]
    energy_range: tuple[float, float]            # trained, MeV
    range_source: str                            # Latin-hypercube bounds: "training" or "custom"
    bounds: dict[str, tuple[float, float]]
    mode: str                                    # "single", "latin_hypercube", or "table"
    points: np.ndarray                           # (point, parameter), in PARAMETER_NAMES order
    energy_source: str                           # "experimental", "grid", or "selected"
    energies: np.ndarray                         # MeV
    include_experiment: bool
    exact_numerov: bool
    report_speedup: bool
    output_dir: Path
    save_cross_sections_csv: bool
    save_summary_json: bool
    save_point_timings_csv: bool


def _table_points(path: Path) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(f"parameter table does not exist: {path}")
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        missing = [name for name in PARAMETER_NAMES if name not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"parameter table {path} has no column {', '.join(missing)}")
        points = [[float(row[name]) for name in PARAMETER_NAMES] for row in reader]
    if not points:
        raise ValueError(f"parameter table {path} has no points")
    return np.asarray(points, dtype=float)


def _points(prediction: dict, bounds: dict) -> tuple[str, np.ndarray]:
    mode = prediction.get("mode", "single")
    if mode == "single":
        single = _table(prediction, "single")
        missing = [name for name in PARAMETER_NAMES if name not in single]
        if missing:
            raise ValueError(f"prediction.single has no {', '.join(missing)}")
        return mode, np.asarray([[single[name] for name in PARAMETER_NAMES]], dtype=float)
    if mode == "latin_hypercube":
        lhs = _table(prediction, "latin_hypercube")
        section = "prediction.latin_hypercube"
        return mode, latin_hypercube([bounds[name] for name in PARAMETER_NAMES],
                                     _integer(lhs, "points", section), _integer(lhs, "seed", section, 0))
    if mode == "table":
        return mode, _table_points(_path(_table(prediction, "table"), "parameter_file", "prediction.table"))
    raise ValueError("prediction.mode must be single, latin_hypercube, or table")


def _energies(table: dict, experimental_file: Path, energy_range: tuple[float, float]) -> tuple[str, np.ndarray]:
    source = table.get("source", "experimental")
    if source == "experimental":
        values = [row["energy_cm_mev"] for row in load_dataset(experimental_file)]
    elif source == "grid":
        grid = table.get("grid")
        if not isinstance(grid, list) or len(grid) != 3:
            raise ValueError("energies.grid must be [first, last, step]")
        first, last, step = map(float, grid)
        if step <= 0.0 or last < first:
            raise ValueError("energies.grid needs step > 0 and last >= first")
        count = int(math.floor((last - first) / step + 1.0e-9)) + 1
        values = [round(first + index * step, 10) for index in range(count)]
    elif source == "selected":
        values = [float(value) for value in table.get("selected", [])]
        if not values:
            raise ValueError("energies.source = selected needs energies.selected")
    else:
        raise ValueError("energies.source must be experimental, grid, or selected")
    if len({round(value, 10) for value in values}) != len(values):
        raise ValueError("the output energies must be unique")
    lower, upper = energy_range
    outside = [value for value in values
               if value < lower - ENERGY_TOLERANCE or value > upper + ENERGY_TOLERANCE]
    if outside:
        raise ValueError(
            f"energies {', '.join(f'{value:g}' for value in outside[:5])} MeV lie outside the "
            f"trained range {lower:g}-{upper:g} MeV; widen energy.range of the training input "
            "and train again")
    return source, np.asarray(values, dtype=float)


def load_rom_config(path: str | Path = DEFAULT_ROM_INPUT) -> RomConfig:
    """Read and check emulator_rom.in; the trained model must exist."""
    source = _existing(path)
    with source.open("rb") as stream:
        payload = tomllib.load(stream)
    model = _table(payload, "model")
    model_id = _model_id(model.get("model_id"))
    manifest = _path(model, "training_root", "model") / "operator_eim" / model_id / "manifest.json"
    if not manifest.is_file():
        raise FileNotFoundError(f"model {model_id} is not trained: {manifest} does not exist. "
                                "Run python -m emulator.train first.")
    trained = read_manifest(manifest)
    settings = trained["settings"]
    if list(settings["parameter_names"]) != list(PARAMETER_NAMES):
        raise ValueError(f"{manifest} has the parameters {settings['parameter_names']}, "
                         f"not {list(PARAMETER_NAMES)}")
    trained_bounds = {name: (float(lower), float(upper))
                      for name, (lower, upper) in zip(PARAMETER_NAMES, settings["parameter_bounds"])}
    energy_range = (float(settings["energy_range_mev"][0]), float(settings["energy_range_mev"][1]))
    experimental_file = project_path(settings["experimental_file"])

    domain = _table(payload, "parameter_domain")
    range_source = domain.get("range_source", "training")
    if range_source == "training":
        bounds = trained_bounds
    elif range_source == "custom":
        bounds = _bounds(_table(domain, "custom"), "parameter_domain.custom")
    else:
        raise ValueError("parameter_domain.range_source must be training or custom")
    mode, points = _points(_table(payload, "prediction"), bounds)
    energy_source, energies = _energies(_table(payload, "energies"), experimental_file, energy_range)

    comparison, output = _table(payload, "comparison"), _table(payload, "output")
    exact_numerov = _boolean(comparison, "exact_numerov", "comparison")
    report_speedup = _boolean(comparison, "report_speedup", "comparison", default=False)
    if report_speedup and not exact_numerov:
        raise ValueError("comparison.report_speedup needs comparison.exact_numerov = true")
    return RomConfig(
        source=source, model_id=model_id, manifest=manifest,
        workers=_integer(model, "workers", "model", default=1),
        reaction=_reaction_label(parse_input(project_path(settings["input_file"]))),
        rank=int(trained["rank"]), experimental_file=experimental_file,
        trained_bounds=trained_bounds, energy_range=energy_range,
        range_source=range_source, bounds=bounds, mode=mode, points=points,
        energy_source=energy_source, energies=energies,
        include_experiment=_boolean(comparison, "include_experiment", "comparison"),
        exact_numerov=exact_numerov, report_speedup=report_speedup,
        output_dir=_path(output, "directory", "output"),
        save_cross_sections_csv=_boolean(output, "save_cross_sections_csv", "output"),
        save_summary_json=_boolean(output, "save_summary_json", "output"),
        save_point_timings_csv=_boolean(output, "save_point_timings_csv", "output", default=True))
