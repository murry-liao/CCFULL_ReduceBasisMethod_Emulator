#!/usr/bin/env python3
"""Offline training of CCFULL-ROM, controlled by emulator_train.in.

1. ``wave_pod``: exact CCFULL snapshots on a Latin hypercube in
   (beta2, beta4, beta6, V0, r0, a0, E, J) and one POD basis for all of them.
2. ``cc_eim``: the multipole CC-EIM of the coupling field and the Galerkin
   projection of every operator onto the basis.

Each stage writes a manifest with its settings.  A finished stage is skipped
while the input still matches; a directory trained with other settings is
refused instead of being overwritten.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from dataclasses import dataclass
import hashlib
from itertools import repeat
import json
import math
from pathlib import Path
import time

import numpy as np

from emulator.eim import (MultipoleEIM, RadialEIM, field_radii, rotor_geometry, rotor_static_field,
                          save_eim)
from emulator.rom import (MODEL_FORMAT, PARAMETER_NAMES, PROJECT_ROOT, ROM, Projector, Reaction,
                          latin_hypercube, project_path, project_relative)
from physical_model import ccfull
from physical_model.inputs import (data_lines, energy_values, load_dataset, numbers, parse_input,
                                   reaction_label, records_sha256)


DEFAULT_INPUT = PROJECT_ROOT / "emulator" / "emulator_train.in"
WAVE_POD_FORMAT = "global_wave_pod_v1"


@dataclass(frozen=True)
class TrainingInput:
    """emulator_train.in; its lines are explained at the end of the file."""

    source: Path                         # this file; lines 1-9 are the CCFULL input of the model
    reaction: str                        # for example "A=16, Z=8 + A=154, Z=62"
    model_id: str                        # line 10
    experimental_file: Path              # line 11
    training_root: Path                  # line 12
    bounds: np.ndarray                   # lines 13-14: [lower, upper] per parameter, PARAMETER_NAMES order
    energy_source: str                   # line 15: "experimental", "ccfull", or "custom"
    energy_range: tuple[float, float]    # MeV
    training_points: int                 # line 16: snapshots, one exact CCFULL solve each
    rank: int                            # POD modes kept for the basis
    seed: int
    workers: int                         # CCFULL processes at a time
    eim_points: int                      # line 17: parameter points of the CC-EIM training
    eim_seed: int
    eim_tolerance: float
    eim_max_rank: int                    # per multipole

    @property
    def wave_pod_dir(self) -> Path:
        return self.training_root / "wave_pod" / self.model_id

    @property
    def cc_eim_dir(self) -> Path:
        return self.training_root / "operator_eim" / self.model_id


def read_input(path: str | Path) -> TrainingInput:
    """Read emulator_train.in line by line."""
    source = Path(path).resolve()
    physical = parse_input(source)
    model_id, data_file, root, deformations, potential, energy, wave_pod, cc_eim = data_lines(source)[9:17]
    if model_id == "auto":
        model_id = f"a{physical['apro']}_z{physical['zpro']}__a{physical['atar']}_z{physical['ztar']}"
    experimental_file = project_path(data_file)
    if energy == "experimental":
        # The measured energies, rounded outward to whole MeV.
        measured = [row["energy_cm_mev"] for row in load_dataset(experimental_file)]
        energy_range = (float(math.floor(min(measured))), float(math.ceil(max(measured))))
    elif energy == "ccfull":
        grid = energy_values(physical)
        energy_range = (float(min(grid)), float(max(grid)))
    else:
        energy, energy_range = "custom", tuple(numbers(energy))
    training_points, rank, seed, workers = map(int, numbers(wave_pod))
    eim_points, eim_seed, eim_tolerance, eim_max_rank = numbers(cc_eim)
    return TrainingInput(
        source=source, reaction=reaction_label(physical), model_id=model_id,
        experimental_file=experimental_file, training_root=project_path(root),
        bounds=np.reshape(numbers(deformations) + numbers(potential), (6, 2)),
        energy_source=energy, energy_range=energy_range,
        training_points=training_points, rank=rank, seed=seed, workers=workers,
        eim_points=int(eim_points), eim_seed=int(eim_seed), eim_tolerance=eim_tolerance,
        eim_max_rank=int(eim_max_rank))


def _write_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    temporary.replace(path)


def _reaction(config: TrainingInput) -> Reaction:
    return Reaction(config.source, config.bounds, config.energy_range)


def _identity(config: TrainingInput) -> dict:
    """The settings that identify the physical model, the training box, and the energy range."""
    return {
        "model_id": config.model_id,
        "input_file": project_relative(config.source),
        "physical_input_sha256": records_sha256(config.source),
        "solver_sha256": ccfull.source_sha256(),
        "experimental_file": project_relative(config.experimental_file),
        "parameter_names": list(PARAMETER_NAMES),
        "parameter_bounds": config.bounds.tolist(),
        "energy_range_mev": list(config.energy_range),
    }


def _trained(manifest: Path, settings: dict) -> bool:
    """True if ``manifest`` was trained with ``settings``; an error if with other settings.

    Only the current settings are compared, so a key that an older manifest
    records but that is no longer a setting does not invalidate it.
    """
    if not manifest.is_file():
        return False
    recorded = json.loads(manifest.read_text()).get("settings", {})
    changed = sorted(key for key in settings if recorded.get(key) != settings[key])
    if not changed:
        return True
    raise ValueError(
        f"{manifest.parent} holds a model trained with other settings ({', '.join(changed)}). "
        "Move that directory out of training_data, or change model_id (line 10), and train again.")


def training_design(reaction: Reaction, count: int, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parameter points, energies, and partial waves of a Latin hypercube in (theta, E, J).

    The last coordinate is uniform on [0, 1) and picks J uniformly from the
    partial waves of the sample's energy.
    """
    bounds = np.vstack([reaction.bounds, [reaction.energy_range], [[0.0, 1.0]]])
    raw = latin_hypercube(bounds, int(count), int(seed))
    theta, energy = raw[:, :-2].copy(), raw[:, -2].copy()
    counts = np.array([reaction.partial_waves(value).size for value in energy])
    waves = np.minimum((raw[:, -1] * counts).astype(int), counts - 1)
    return theta, energy, waves


def _snapshot(config: dict, point: np.ndarray, energy: float, angular_momentum: int) -> tuple:
    """One exact CCFULL solve, normalized to unit norm and stored as complex64."""
    result = ccfull.snapshot(config, point, energy, angular_momentum)
    vector = np.asarray(result["wave"], dtype=np.complex128)
    norm = float(np.linalg.norm(vector))
    if not np.isfinite(norm) or norm == 0.0:
        raise RuntimeError(f"invalid snapshot at E = {energy} MeV, J = {angular_momentum}: norm {norm}")
    return (vector / norm).astype(np.complex64), result


def train_wave_pod(config: TrainingInput) -> Path:
    """Stage 1: exact snapshots on a Latin hypercube in (theta, E, J), and their POD basis."""
    directory = config.wave_pod_dir
    manifest = directory / "manifest.json"
    settings = {**_identity(config), "training_points": config.training_points,
                "rank": config.rank, "seed": config.seed}
    if _trained(manifest, settings):
        print(f"[wave_pod] already trained: {project_relative(manifest)}", flush=True)
        return manifest
    reaction = _reaction(config)
    count = config.training_points
    print(f"[wave_pod] {count} snapshots (exact CCFULL solves) in {config.workers} processes", flush=True)
    started = time.perf_counter()
    theta, energy, waves = training_design(reaction, count, config.seed)
    ccfull.snapshot_program()            # compile once, before the CCFULL processes start
    vectors, first = [], None
    step = max(1, count // 10)
    # Each thread waits on one CCFULL process.
    with ThreadPoolExecutor(max_workers=config.workers) as pool:
        solved = pool.map(_snapshot, repeat(reaction.config), theta, energy, waves)
        for done, (vector, result) in enumerate(solved, start=1):
            vectors.append(vector)
            first = first or result
            if done % step == 0 or done == count:
                print(f"[wave_pod] snapshots {done}/{count}", flush=True)
    radius, channels = first["radius"], int(vectors[0].shape[1])
    vectors = np.stack(vectors).reshape(count, -1)
    snapshot_seconds = time.perf_counter() - started
    if config.rank > min(vectors.shape):
        raise ValueError(f"the basis rank {config.rank} (line 16) exceeds the {min(vectors.shape)} "
                         "POD modes of these snapshots")
    started = time.perf_counter()
    left, values, _ = np.linalg.svd(vectors.T.astype(np.complex128), full_matrices=False)
    pod_seconds = time.perf_counter() - started
    tail = np.sqrt(np.maximum(1.0 - np.cumsum(values**2) / np.sum(values**2), 0.0))
    directory.mkdir(parents=True, exist_ok=True)
    np.savez(directory / "snapshots.npz", theta=theta, energy_mev=energy, j=waves, vectors=vectors)
    np.savez(directory / "basis.npz", basis=left[:, :config.rank], singular_values=values, radius=radius)
    with (directory / "pod_singular_values.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["mode", "singular_value", "relative_tail"])
        writer.writerows([mode, f"{value:.10e}", f"{rest:.10e}"]
                         for mode, (value, rest) in enumerate(zip(values, tail), start=1))
    _write_json(manifest, {
        "format": WAVE_POD_FORMAT,
        "settings": settings,
        "snapshots": "snapshots.npz",
        "basis": "basis.npz",
        "singular_values": "pod_singular_values.csv",
        "radial_points": int(radius.size),
        "channels": channels,
        "reduced_mass": first["reduced_mass"],
        "zz": first["zz"],
        "relative_pod_tail_at_rank": float(tail[config.rank - 1]),
        "seconds": {"snapshots": snapshot_seconds, "pod": pod_seconds},
    })
    print(f"[wave_pod] basis: {config.rank} of {values.size} POD modes, relative tail "
          f"{tail[config.rank - 1]:.1e}; snapshots {snapshot_seconds:.0f} s, POD {pod_seconds:.0f} s "
          f"-> {project_relative(manifest)}", flush=True)
    return manifest


def train_eim(reaction: Reaction, radius: np.ndarray, channels: int, config: TrainingInput) -> tuple:
    """The multipole CC-EIM from potential evaluations at parameter points; no wave solves."""
    points = latin_hypercube(reaction.bounds, config.eim_points, config.eim_seed)
    multipoles, geometries = rotor_geometry(reaction.config, channels)
    radii = field_radii(radius)
    static = rotor_static_field(reaction.config, radii.size, channels)
    inverse = np.linalg.pinv(geometries.reshape(len(multipoles), -1).T)
    factors = np.empty((len(points), radii.size, len(multipoles)), dtype=np.float64)
    # The complete field of the first points checks the multipole decomposition.
    errors = np.empty(min(16, len(points)))
    started = time.perf_counter()
    for index, point in enumerate(points):
        factors[index] = reaction.library.factors(point, radii, inverse)
        if index < errors.size:
            residual = reaction.library.fields(point, radii) - static
            reconstructed = np.einsum("rm,mij->rij", factors[index], geometries, optimize=True)
            errors[index] = (np.linalg.norm(reconstructed - residual)
                             / max(float(np.linalg.norm(residual)), 1.0e-300))
    radial = []
    for component in range(len(multipoles)):
        samples = factors[:, :, component]
        singular = np.linalg.svd(samples.T, compute_uv=False)
        kept = np.count_nonzero(singular > config.eim_tolerance * max(float(singular[0]), 1.0e-300))
        radial.append(RadialEIM.fit(samples, max(1, min(config.eim_max_rank, int(kept)))))
    eim = MultipoleEIM(multipoles, geometries, tuple(radial), static, inverse,
                       float(np.median(errors)), float(np.max(errors)))
    diagnostics = {
        "interaction_training_points": config.eim_points,
        "interaction_seed": config.eim_seed,
        "tolerance": config.eim_tolerance,
        "max_rank_per_multipole": config.eim_max_rank,
        "component_ranks": {str(label): rank for label, rank in zip(multipoles, eim.ranks)},
        "total_eim_rank": eim.rank,
        "selected_radial_nodes_total": int(eim.node_rows.size),
        "decomposition_median_relative_error": eim.decomposition_median_relative_error,
        "decomposition_max_relative_error": eim.decomposition_max_relative_error,
        "seconds": time.perf_counter() - started,
    }
    return eim, diagnostics


def train_cc_eim(config: TrainingInput) -> Path:
    """Stage 2: CC-EIM of the coupling field and the Galerkin projection of every operator."""
    parent = config.wave_pod_dir / "manifest.json"
    directory = config.cc_eim_dir
    manifest = directory / "manifest.json"
    settings = {**_identity(config),
                "parent_manifest": project_relative(parent),
                "parent_manifest_sha256": hashlib.sha256(parent.read_bytes()).hexdigest(),
                "interaction_points": config.eim_points,
                "interaction_seed": config.eim_seed,
                "tolerance": config.eim_tolerance,
                "max_rank_per_multipole": config.eim_max_rank}
    if _trained(manifest, settings):
        print(f"[cc_eim] already trained: {project_relative(manifest)}", flush=True)
        return manifest
    wave_pod = json.loads(parent.read_text())
    with np.load(parent.parent / wave_pod["basis"]) as data:
        basis, radius = data["basis"], data["radius"]
    reaction = _reaction(config)
    print(f"[cc_eim] coupling field at {config.eim_points} parameter points", flush=True)
    eim, diagnostics = train_eim(reaction, radius, wave_pod["channels"], config)
    directory.mkdir(parents=True, exist_ok=True)
    save_eim(eim, directory / "shared_interaction")
    (directory / "shared_interaction" / "diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n")
    print(f"[cc_eim] {eim.rank} modes over {len(eim.multipoles)} multipoles, "
          f"{diagnostics['selected_radial_nodes_total']} radial nodes", flush=True)
    started = time.perf_counter()
    projector = Projector(basis, radius, wave_pod["reduced_mass"], wave_pod["zz"], wave_pod["channels"])
    rom = ROM(reaction, eim, radius, wave_pod["reduced_mass"], wave_pod["zz"], projector.operators(eim))
    projection_seconds = time.perf_counter() - started
    # The projected CC-EIM modes must reproduce the projection of the reconstructed CC-EIM field.
    via_modes = rom.interaction(reaction.nominal)[0]
    via_field = projector.field_operator(eim.reconstruct(rom.coefficients(reaction.nominal)))
    check = float(np.linalg.norm(via_modes - via_field) / np.linalg.norm(via_field))
    rom.save(directory / "reduced_operators.npz")
    _write_json(manifest, {
        "format": MODEL_FORMAT,
        "settings": settings,
        "shared_eim": "shared_interaction",
        "operators": "reduced_operators.npz",
        "rank": rom.rank,
        "interaction_training": diagnostics,
        "projection_check_relative_difference": check,
        "seconds": {"cc_eim": diagnostics["seconds"], "projection": projection_seconds},
    })
    print(f"[cc_eim] operators of {eim.rank} CC-EIM modes projected onto the rank-{rom.rank} basis "
          f"in {projection_seconds:.0f} s; projection check {check:.1e} -> {project_relative(manifest)}",
          flush=True)
    return manifest


def _report(config: TrainingInput) -> None:
    """What each line of the input was read as."""
    names = ("beta2", "beta4", "beta6", "V0", "r0", "a0")
    bounds = [f"{name} [{lower:g}, {upper:g}]" for name, (lower, upper) in zip(names, config.bounds)]
    lower, upper = config.energy_range
    print(f"Training input {config.source}")
    print(f"  lines 1-9  CCFULL input of {config.reaction}")
    print(f"  line 10    model_id {config.model_id}")
    print(f"  line 11    experimental data {project_relative(config.experimental_file)}")
    print(f"  line 12    training_root {project_relative(config.training_root)}")
    print(f"  line 13    {', '.join(bounds[:3])}")
    print(f"  line 14    {', '.join(bounds[3:])}")
    print(f"  line 15    energy range {lower:g}-{upper:g} MeV ({config.energy_source})")
    print(f"  line 16    {config.training_points} snapshots, basis rank {config.rank}, "
          f"seed {config.seed}, {config.workers} CCFULL process(es)")
    print(f"  line 17    CC-EIM: {config.eim_points} parameter points, seed {config.eim_seed}, "
          f"tolerance {config.eim_tolerance:g}, at most {config.eim_max_rank} modes per multipole")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--stage", choices=("all", "wave_pod", "cc_eim"), default="all")
    parser.add_argument("--check", action="store_true",
                        help="print what each line of the input was read as, without training")
    args = parser.parse_args()
    config = read_input(args.config)
    if args.check:
        _report(config)
        return
    if args.stage in ("all", "wave_pod"):
        train_wave_pod(config)
    if args.stage in ("all", "cc_eim"):
        train_cc_eim(config)


if __name__ == "__main__":
    main()
