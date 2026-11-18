"""Train a global reduced model: one basis for every (E, J) of an energy range.

    python -m emulator.train_global                           # reads emulator/emulator_train.in
    python -m emulator.train_global --config <file> --check   # only prints how each line was read

Steps (the lines of emulator_train.in are explained at its end):
  1. Read the input and name the model.
  2. The radial grid: the interface program at the nominal point puts r_min at the bottom of
     its pocket; every later run keeps this r_min.
  3. Snapshots: a Latin hypercube in (beta2, beta4, beta6, V0, r0, a0, E, u) (line 16).  At each
     point the partial wave J = [u n], of the n partial waves that CCFULL's J loop reaches at
     energy E at a corner of the (V0, r0, a0) box, and one run of the interface program at
     (theta, E) that solves it (line 11 of ccfull.inp).
  4. The POD basis of the snapshots (line 16).
  5. The CC-EIM of the coupling field (line 17).
  6. The Numerov equations projected onto the basis; the outer boundary is set online, for
     the energies of each run (rom.prepare_global).
  7. The files and the manifest, and copies of the CCFULL records and of the interface program.
"""

import argparse
import csv
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
    parser = argparse.ArgumentParser(description="Train a global reduced model of CCFULL.")
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
    count, rank, seed = (int(value) for value in numbers(lines[15]))
    eim_points, eim_seed, eim_modes = (int(value) for value in numbers(lines[16]))
    name = f"{lines[9]}_Global_Channel{system['channels']}_Snap{count}_POD{rank}_EIM{eim_modes}"
    folder = training_root / name

    print(f"Training input {args.config}")
    print(f"  lines 1-9    CCFULL records of A={system['apro']}, Z={system['zpro']} + A={system['atar']}, "
          f"Z={system['ztar']}, {system['channels']} channels")
    print(f"  lines 10-12  model {name} in {training_root}")
    print("  lines 13-14  box: " + ", ".join(f"{n} [{low:g}, {high:g}]" for n, (low, high) in zip(rom.PARAMETER_NAMES, box)))
    print(f"  line 15      energy range {energy_range[0]:g}-{energy_range[1]:g} MeV ({lines[14]})")
    print(f"  line 16      {count} snapshots, basis rank {rank}, seed {seed}")
    print(f"  line 17      CC-EIM: {eim_points} points, seed {eim_seed}, {eim_modes} modes")
    if args.check:
        return
    if folder.exists():
        raise SystemExit(f"{folder} exists: move it away, or change the system name (line 10)")
    started = time.perf_counter()

    # 2. The radial grid, from the interface program at the nominal point.
    target, potential = numbers(records[2]), numbers(records[6])
    nominal = np.array([target[1], target[2], target[3], potential[0], potential[1], potential[2]])
    first_run = wave_functions(records, nominal, energy_range[0], None)
    radius = first_run["radius"]
    r_min = radius[0]
    print(f"[grid] r_min {r_min:.6f} fm (the pocket at the nominal point), {radius.size} points, "
          f"dr {radius[1] - radius[0]:g} fm", flush=True)

    # 3. Snapshots: one run of the interface program per Latin-hypercube point.
    design = rom.latin_hypercube(np.vstack([box, [energy_range], [[0.0, 1.0]]]), count, seed)
    theta, energy, u = design[:, :6], design[:, 6], design[:, 7]
    potential_corners = list(itertools.product(*box[3:]))
    J = np.zeros(count, dtype=int)
    vectors = np.empty((count, radius.size * system["channels"]), dtype=np.complex64)
    for k in range(count):
        n = max(partial_wave_count(rom.ccfull_energy(energy[k]), v0, r0, a0, system) for v0, r0, a0 in potential_corners)
        J[k] = min(int(u[k] * n), n - 1)
        run = wave_functions(records, theta[k], energy[k], r_min, J[k])
        vector = run["waves"][J[k]].reshape(-1)
        vectors[k] = vector / np.linalg.norm(vector)
        if (k + 1) % max(1, count // 10) == 0:
            print(f"[snapshots] {k + 1}/{count}; {time.perf_counter() - started:.0f} s", flush=True)
    snapshot_seconds = time.perf_counter() - started

    # 4. The POD basis.
    left, singular, _ = np.linalg.svd(vectors.T.astype(np.complex128), full_matrices=False)
    tail = np.sqrt(np.maximum(1.0 - np.cumsum(singular ** 2) / np.sum(singular ** 2), 0.0))
    basis = left[:, :rank]
    print(f"[POD] rank {rank} of {singular.size}, relative tail {tail[rank - 1]:.1e}", flush=True)

    # 5. The CC-EIM.
    eim, diagnostics = rom.fit_eim(system, radius, box, eim_points, eim_seed, eim_modes)
    print(f"[CC-EIM] {eim['basis'].shape[2]} modes, {eim['node_rows'].size} node rows; error at 16 new points: "
          f"median {diagnostics['interpolation_relative_error']['median']:.1e}, "
          f"max {diagnostics['interpolation_relative_error']['max']:.1e}", flush=True)

    # 6. The projected Numerov equations.  Check: the projected modes summed at the nominal point
    #    against the projection of the CC-EIM field summed there.
    fields = rom.mode_fields(eim)[:, :radius.size]
    mu, zz = first_run["reduced_mass"], system["zpro"] * system["ztar"]
    operators = rom.project_numerov(basis, radius, mu, zz, system["channels"], eim["static_field"], fields)
    model = {"system": system, "channels": system["channels"], "eim": eim}
    c = rom.eim_coefficients(model, nominal)
    summed_field = eim["static_field"][:radius.size] + np.einsum("m,mrij->rij", c, fields)
    direct = rom.project_numerov(basis, radius, mu, zz, system["channels"], 0.0 * summed_field,
                                 [summed_field])["a_modes"][0]
    via_modes = operators["a_static"] + np.tensordot(c, operators["a_modes"], axes=1)
    check = float(np.linalg.norm(via_modes - direct) / np.linalg.norm(direct))
    print(f"[projection] {len(fields)} modes onto rank {rank}; check {check:.1e}", flush=True)

    # 7. Files and manifest.
    folder.mkdir(parents=True)
    rom.save_eim(eim, diagnostics, folder / "shared_interaction")
    np.savez(folder / "snapshots.npz", theta=theta, energy_mev=energy, j=J, vectors=vectors)
    np.savez(folder / "basis.npz", basis=basis, singular_values=singular, radius=radius)
    with open(folder / "pod_singular_values.csv", "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["mode", "singular_value", "relative_tail"])
        writer.writerows([k + 1, f"{s:.10e}", f"{t:.10e}"] for k, (s, t) in enumerate(zip(singular, tail)))
    np.savez(folder / "reduced_operators.npz", radius=radius, reduced_mass=mu, zz=zz, **operators)
    (folder / "ccfull.inp").write_text("\n".join(records) + f"\n{r_min!r}\n")
    shutil.copy(PROJECT / "physical_model" / "ccfull_interface" / "ccfull.f90", folder / "ccfull_interface.f90")
    seconds = time.perf_counter() - started
    manifest = {
        "format": "global_rom_v3",
        "settings": {"model_id": name, "input_file": os.path.relpath(os.path.abspath(args.config), PROJECT),
                     "experimental_file": lines[10], "parameter_names": list(rom.PARAMETER_NAMES),
                     "parameter_bounds": box.tolist(), "energy_range_mev": energy_range,
                     "snapshots": count, "rank": rank, "seed": seed,
                     "interaction_points": eim_points, "interaction_seed": eim_seed, "eim_modes": eim_modes},
        "shared_eim": "shared_interaction", "operators": "reduced_operators.npz", "rank": rank,
        "interaction_training": diagnostics, "relative_pod_tail_at_rank": float(tail[rank - 1]),
        "projection_check_relative_difference": check,
        "seconds": {"total": seconds, "snapshots": snapshot_seconds},
    }
    (folder / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"[global] snapshots {snapshot_seconds:.0f} s, total {seconds:.0f} s -> {folder}", flush=True)


if __name__ == "__main__":
    main()
