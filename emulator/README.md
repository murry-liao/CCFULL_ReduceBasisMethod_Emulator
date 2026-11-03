# Emulator

The offline and online stages are separate:

1. `training_data/wave_pod/` stores the exact CCFULL snapshots of a Latin
   hypercube in (β2, β4, β6, V0, r0, a0, E, J) and the single POD basis built
   from them.
2. `training_data/operator_eim/` stores the multipole CC-EIM of the coupling
   field and every operator projected onto that basis. This is the model used
   online.
3. `outputs/` stores timestamped cross sections for new parameter queries.
4. `analysis/` stores the parameter-analysis notebook and system-specific
   analysis scripts.

`configs/emulator_train.in` is the single offline setup for the physical input,
training box, energy range, snapshots and basis, and CC-EIM stage.
`configs/emulator_rom.in` names the trained model and controls the online
workers, allowed query ranges, default point, output energies, accuracy check,
and output directory.

Training snapshots are stored as `complex64`; the exact Numerov propagation
remains `complex128`.

## Code

- `train.py` runs the two offline stages: snapshots and POD basis, then the
  CC-EIM and the projected operators.
- `run_rom.py` evaluates the parameter points of `emulator_rom.in`, runs the
  accuracy check, and writes the outputs.
- `config.py` reads the two input files and the manifests of trained models.
- `rom.py` is the reduced model: the Galerkin projection of the Numerov
  operators (offline) and the solve for σ_fus (online).
- `eim.py` is the multipole CC-EIM of the coupling field.

CCFULL itself is called through `../physical_model/ccfull.py`.

## Configured-system workflow

Run these commands from the project root, `CCFULL_ROM-new`:

```bash
# 1. Validate, then run the configured offline stages.
python -m emulator.train --check
python -m emulator.train

# 2. Validate the online query; run it after training exists.
python -m emulator.run_rom --check
python -m emulator.run_rom
```

System-specific fitting and comparison scripts belong under `analysis/`,
rather than in the reusable emulator package.

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
`configs/emulator_train.in`; a model is refused once either differs from its
training run.

Generated training and result directories are not source inputs. Their
contents may be removed without deleting the physical code, configuration
files, experimental data, or input templates.
