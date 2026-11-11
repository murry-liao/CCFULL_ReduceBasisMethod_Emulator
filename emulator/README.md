# Emulator

The offline and online stages are separate:

1. `training_data/wave_pod/` stores the exact CCFULL snapshots of a Latin
   hypercube in (β2, β4, β6, V0, r0, a0, E, J) and the single POD basis built
   from them.
2. `training_data/operator_eim/` stores the multipole CC-EIM of the coupling
   field and every operator projected onto that basis. This is the model used
   online.
3. `outputs/` stores timestamped cross sections for new parameter queries.
4. `parameter_analysis.ipynb` analyzes one run: the best parameter point by
   χ², and the excitation function and barrier distribution at that point.

Training snapshots are stored as `complex64`; the exact Numerov propagation
remains `complex128`.

## Input files

`emulator_train.in` is the single offline setup: lines 1–9 are the unchanged
CCFULL input of the reaction, and lines 10–17 hold the model name, the
experimental data, the directory of trained models, the training box, the
energy range, the snapshots and basis, and the CC-EIM. `emulator_rom.in` has
16 lines with the trained model, the workers, the parameter points, the output
energies, the accuracy check against exact CCFULL, and the output files; the
trained box, energy range, and experimental data come from the manifest of the
model.

Both files are read line by line, like a CCFULL input: each line holds one
item or a few comma-separated numbers, and blank lines and lines starting with
`#` are skipped. Every line is explained at the end of its file, below the
dashed line; nothing after the last data line is read. Relative paths are
resolved from the project root, `CCFULL_ROM-new`, independent of the launch
directory. The values are used as written, without a separate validation, so
print how each line was read before a long run (`--check` below).

Several models of one reaction, for example of different ranks, can be kept
side by side under distinct model_ids. Line 1 of `emulator_rom.in` chooses the
model, so switching models needs no change of `emulator_train.in`.

To change the reaction system, replace lines 1–9 of `emulator_train.in` with
the new CCFULL input, and set the experimental data (line 11), the training
box (lines 13–14), and the energy range (line 15). Keep `auto` on line 10, or
give a new model_id. After training, put that model_id on line 1 of
`emulator_rom.in`. No Python source edit is required.

## Code

- `train.py` reads `emulator_train.in` and runs the two offline stages:
  snapshots and POD basis, then the CC-EIM and the projected operators.
- `run_rom.py` reads `emulator_rom.in`, evaluates its parameter points, runs
  the accuracy check, and writes the outputs.
- `rom.py` is the reduced model: the Galerkin projection of the Numerov
  operators (offline) and the solve for σ_fus (online). It also holds what both
  programs share: the parameter names, project paths, Latin-hypercube sampling,
  and the reading of trained-model manifests.
- `eim.py` is the multipole CC-EIM of the coupling field.

CCFULL itself is called through `../physical_model/ccfull.py`.

## Configured-system workflow

Run these commands from the project root, `CCFULL_ROM-new`:

```bash
# 1. Print how each line of emulator_train.in is read, then train.
python -m emulator.train --check
python -m emulator.train                     # both offline stages
python -m emulator.train --stage wave_pod    # snapshots and basis only
python -m emulator.train --stage cc_eim      # CC-EIM and projection only

# 2. Print how each line of emulator_rom.in is read; run it after training exists.
python -m emulator.run_rom --check
python -m emulator.run_rom
```

Each run creates one timestamped directory in `outputs/`. All parameter points
are written to the same `cross_sections.csv`; `point_index` and the six
parameter columns identify each curve. With the exact Numerov check (line 14
of `emulator_rom.in`), the file also holds `exact_numerov_cross_section_mb`
and `rom_exact_relative_error`, and the command prints the median, 95th and
99th percentile, and maximum of |σ_ROM/σ_exact − 1|.

## Generated-directory contract

- `training_data/wave_pod/<model-id>/` contains `snapshots.npz` (the design
  and the unit-norm snapshots), `basis.npz` (the basis, all singular values,
  and the radial grid), `pod_singular_values.csv` (each singular value and the
  relative POD tail after it), and `manifest.json`.
- `training_data/operator_eim/<model-id>/` contains `shared_interaction/` (the
  multipole CC-EIM), `reduced_operators.npz` (the projected operators and the
  basis rows at the radial boundaries), and `manifest.json`.
- `outputs/<timestamp>_<model-id>/` contains one combined
  `cross_sections.csv`, `run_summary.json`, and `point_timings.csv`. Batch
  points are distinguished by `point_index` and their six parameter columns.

Trained manifests record their input and parent files relative to the project
root, so the project directory may be moved. They also record the SHA-256 of
the standard Numerov source and of the nine CCFULL records of
`emulator_train.in`, and a model is refused once either changes. Each training
stage also records its settings: a stage whose settings still match is
skipped, and a directory trained with other settings is refused instead of
being overwritten. Move it out of `training_data/`, or set another model_id on
line 10, and train again.

Generated training and result directories are not source inputs. Their
contents may be removed without deleting the physical code, input files,
experimental data, or input templates.
