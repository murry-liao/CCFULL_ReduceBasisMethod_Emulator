"""Train a local reduced model: one small basis for every (E, J) at chosen energies.

    python -m emulator.train_local                            # reads emulator/emulator_train.in
    python -m emulator.train_local --config <file> --check    # only prints how each line was read

Steps (the lines of emulator_train.in are explained at its end):
  1. Read the input and name the model.
  2. The radial grid: the interface program at the nominal point puts r_min at the bottom of
     its pocket; every later run keeps this r_min.
  3. The CC-EIM of the coupling field (line 17).
  4. Every energy of line 18: the interface program at the Latin-hypercube points of line 19,
     every point solving the same partial waves; the corners of the box set how many are kept.
     For every J, the POD of its snapshots and the Numerov equations projected onto the
     leading vectors; the outer boundary comes from a run at the nominal point.  One file per
     energy.
  5. The manifest, and copies of the CCFULL records and of the interface program.
"""

import argparse
import itertools
import json
import math
import os
import shutil
import time

for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(variable, "1")

import numpy as np

from emulator import rom
from emulator.ccfull_run import PROJECT, data_lines, numbers, project_path, read_experimental, wave_functions
from emulator.coupling import partial_wave_count, reaction


def main():
    parser = argparse.ArgumentParser(description="Train a local reduced model of CCFULL.")
    parser.add_argument("--config", default=PROJECT / "emulator" / "emulator_train.in")
    parser.add_argument("--check", action="store_true", help="print how each line was read, and stop")
    args = parser.parse_args()

    # 1. Read emulator_train.in.
    lines = data_lines(args.config)
    records = lines[:9]
    system = reaction(records)
    measured = read_experimental(project_path(lines[10]))[0]
    training_root = project_path(lines[11])
    box = np.reshape(numbers(lines[12]) + numbers(lines[13]), (6, 2))
    if lines[14] == "experimental":
        energy_range = [float(math.floor(measured.min())), float(math.ceil(measured.max()))]
    elif lines[14] == "ccfull":
        energy_range = numbers(records[7])[:2]
    else:
        energy_range = numbers(lines[14])
    eim_points, eim_seed, eim_modes = (int(value) for value in numbers(lines[16]))
    if lines[17] == "none":
        raise SystemExit("line 18 of the input is none: no local model to train")
    if lines[17] == "experimental":
        energies = measured
    else:
        first, last, step = numbers(lines[17])
        energies = np.array([round(first + k * step, 10) for k in range(int(math.floor((last - first) / step + 1e-9)) + 1)])
    rank, count, seed = (int(value) for value in numbers(lines[18]))
    name = f"{lines[9]}_Local_Channel{system['channels']}_Snap{count}_POD{rank}_EIM{eim_modes}"
    folder = training_root / name
    corners = np.unique(np.array(list(itertools.product(*box))), axis=0)

    print(f"Training input {args.config}")
    print(f"  lines 1-9    CCFULL records of A={system['apro']}, Z={system['zpro']} + A={system['atar']}, "
          f"Z={system['ztar']}, {system['channels']} channels")
    print(f"  lines 10-12  model {name} in {training_root}")
    print("  lines 13-14  box: " + ", ".join(f"{n} [{low:g}, {high:g}]" for n, (low, high) in zip(rom.PARAMETER_NAMES, box)))
    print(f"  line 17      CC-EIM: {eim_points} points, seed {eim_seed}, {eim_modes} modes")
    print(f"  line 18      {len(energies)} energies from {min(energies):g} to {max(energies):g} MeV")
    print(f"  line 19      rank {rank} per (E, J); {count} Latin-hypercube points (seed {seed}) "
          f"per energy; the {len(corners)} corners of the box set the partial waves kept")
    if args.check:
        return
    if folder.exists():
        raise SystemExit(f"{folder} exists: move it away, or change the system name (line 10)")
    started = time.perf_counter()

    # 2. The radial grid, from the interface program at the nominal point.
    target, potential = numbers(records[2]), numbers(records[6])
    nominal = np.array([target[1], target[2], target[3], potential[0], potential[1], potential[2]])
    first_run = wave_functions(records, nominal, energies[0], None)
    radius = first_run["radius"]
    r_min = radius[0]
    print(f"[grid] r_min {r_min:.6f} fm (the pocket at the nominal point), {radius.size} points, "
          f"dr {radius[1] - radius[0]:g} fm", flush=True)

    # 3. The CC-EIM, and the fields of its modes on the grid.
    eim, diagnostics = rom.fit_eim(system, radius, box, eim_points, eim_seed, eim_modes)
    fields = rom.mode_fields(eim)[:, :radius.size]
    grid = {"radius": radius, "reduced_mass": first_run["reduced_mass"], "zz": system["zpro"] * system["ztar"],
            "channels": system["channels"], "system": system, "static_field": eim["static_field"]}
    print(f"[CC-EIM] {eim['basis'].shape[2]} modes, {eim['node_rows'].size} node rows; error at 16 new points: "
          f"median {diagnostics['interpolation_relative_error']['median']:.1e}, "
          f"max {diagnostics['interpolation_relative_error']['max']:.1e}", flush=True)
    folder.mkdir(parents=True)
    rom.save_eim(eim, diagnostics, folder / "shared_interaction")
    np.savez(folder / "grid.npz", radius=radius, reduced_mass=grid["reduced_mass"], zz=grid["zz"])

    # 4. Every energy.  The partial waves are J = 0 .. n-1, n the largest count of CCFULL's J loop
    #    over the corners of the (V0, r0, a0) box, and every point solves all of them (line 11 of
    #    ccfull.inp), so that every block has a snapshot at every training point.  A block is kept
    #    for each J that the J sum reaches at a training point or at a corner of the box, and 5 more.
    points = rom.latin_hypercube(box, count, seed)
    potential_corners = list(itertools.product(*box[3:]))
    channels = system["channels"]
    files, tails, ccfull_seconds = [], [], 0.0
    for index, energy in enumerate(energies):
        solving = time.perf_counter()
        e = rom.ccfull_energy(energy)
        n = max(partial_wave_count(e, v0, r0, a0, system) for v0, r0, a0 in potential_corners)
        reached = 0
        for theta in corners:
            run = wave_functions(records, theta, energy, r_min, n - 1)
            reached = max(reached, rom.partial_waves_in_sum(run["P"], run["J"]))
        snapshots = np.empty((n, count, radius.size * channels), dtype=np.complex64)    # (J, point, r x channel)
        for k, theta in enumerate(points):
            run = wave_functions(records, theta, energy, r_min, n - 1)
            reached = max(reached, rom.partial_waves_in_sum(run["P"], run["J"]))
            vectors = run["waves"][:n].reshape(n, -1)
            snapshots[:, k] = vectors / np.linalg.norm(vectors, axis=1)[:, None]
        nominal_run = wave_functions(records, nominal, energy, r_min)
        ccfull_seconds += time.perf_counter() - solving
        kept = min(n, reached + 5)
        ratio, source = rom.outer_boundary(grid, nominal_run["matching"], kept)
        blocks = [rom.local_block(snapshots[J], J, e, rank, grid, fields, ratio[J], source[J]) for J in range(kept)]
        file = f"energy_{index:03d}.npz"
        np.savez(folder / file, energy=float(energy), waves=np.arange(kept),
                 a0=np.array([b["a0"] for b in blocks]),
                 a_modes=np.ascontiguousarray(np.array([b["a_modes"] for b in blocks]).transpose(1, 0, 2, 3),
                                              dtype=np.complex64),
                 x0=np.array([b["x0"] for b in blocks]), x1=np.array([b["x1"] for b in blocks]),
                 rhs=np.array([b["rhs"] for b in blocks]), tail=np.array([b["tail"] for b in blocks]))
        files.append({"energy_mev": float(energy), "file": file, "partial_waves": kept, "partial_waves_solved": n})
        tails.extend(b["tail"] for b in blocks)
        print(f"[local] {index + 1}/{len(energies)}  E = {energy:.4f} MeV: {kept} of {n} partial waves, largest "
              f"POD tail {max(b['tail'] for b in blocks):.1e}; {time.perf_counter() - started:.0f} s", flush=True)

    # 5. The manifest, and copies of the CCFULL records and of the interface program.
    (folder / "ccfull.inp").write_text("\n".join(records) + f"\n{r_min!r}\n")
    shutil.copy(PROJECT / "physical_model" / "ccfull_interface" / "ccfull.f90", folder / "ccfull_interface.f90")
    seconds = time.perf_counter() - started
    manifest = {
        "format": "local_rom_v4",
        "settings": {"model_id": name, "input_file": os.path.relpath(os.path.abspath(args.config), PROJECT),
                     "experimental_file": lines[10],
                     "parameter_names": list(rom.PARAMETER_NAMES), "parameter_bounds": box.tolist(),
                     "energy_range_mev": energy_range, "local_energies_mev": [float(e) for e in energies],
                     "local_rank": rank, "local_snapshots": count, "local_seed": seed,
                     "local_partial_wave_margin": 5,
                     "interaction_points": eim_points, "interaction_seed": eim_seed, "eim_modes": eim_modes},
        "shared_eim": "shared_interaction", "grid": "grid.npz", "rank": rank,
        "interaction_training": diagnostics,
        "relative_pod_tail_at_rank": {"median": float(np.median(tails)), "max": float(np.max(tails))},
        "blocks": len(tails), "files": files,
        "seconds": {"total": seconds, "ccfull": ccfull_seconds},
    }
    (folder / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"[local] {len(tails)} blocks; POD tail at rank {rank}: median {np.median(tails):.1e}, max {np.max(tails):.1e}; "
          f"CCFULL {ccfull_seconds:.0f} s, total {seconds:.0f} s -> {folder}", flush=True)


if __name__ == "__main__":
    main()
