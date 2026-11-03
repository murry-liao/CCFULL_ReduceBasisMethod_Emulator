# Emulator input files

Relative paths are resolved from the project root, `CCFULL_ROM-new`,
independent of the launch directory.

- `emulator_train.in` starts with the nine unchanged CCFULL records and then
  appends TOML sections for the model identity, experimental data,
  six-parameter training box, trained energy range, snapshots and basis, and
  every multipole CC-EIM construction setting.
- `emulator_rom.in` names the trained model and selects the worker count,
  parameter domain, single or batch prediction mode, output energies, the
  accuracy check against exact CCFULL, speedup reporting, experimental overlay,
  and output files.

## Snapshots and basis

`[wave_pod]` holds the two numbers that set the size of the model.

| Key | Meaning |
|---|---|
| `training_points` | Number of snapshots, the training data. Each is a point of a Latin hypercube in (β2, β4, β6, V0, r0, a0, E, J) with one exact CCFULL solve. |
| `rank` | Number of basis vectors: the POD modes kept from the snapshots, at most `training_points`. Every online (E, J) solve is a `rank` × `rank` system, so the online time grows with the rank. |

Each snapshot draws J uniformly from the partial waves that CCFULL needs at
its energy: J = 0 … n − 1, where n is the largest CCFULL partial-wave count
over the corners of the (V0, r0, a0) box.

## Energy range

The energy is a model input, so one model covers a whole interval.
`[energy] range` selects it:

- `"experimental"`: the measured energies, rounded outward to whole MeV
  (52–100 MeV for the default data);
- `"ccfull"`: Emin and Emax of CCFULL record 8;
- an interval in MeV, for example `range = [52.0, 100.0]`.

Online, any energy inside the trained range can be requested. An energy
outside it is rejected by `--check`.

## Trained models

A trained model records the SHA-256 of the nine CCFULL records in
`emulator_train.in` and of the standard Numerov source, and is refused once
either changes. Each stage also records its settings. A stage whose settings
still match the configuration is skipped. A directory trained with other
settings is refused instead of being overwritten: move it out of
`training_data/`, or set a distinct `model_id`, and train again.

For any parameter bound, equal endpoints freeze that coordinate. For example,
`beta6 = [0.6, 0.6]` keeps `beta6` equal to `0.6` at every training or online
Latin-hypercube point. The matching value in the nine-record CCFULL input and
the online single-point block must also be `0.6`. Parameters with unequal
endpoints continue to vary normally.

## Online queries

`prediction.mode` accepts `single`, `latin_hypercube`, or `table`. The table
template is `data/input_templates/rom_parameter_points.csv`.

Online extrapolation in the parameters is allowed. A `prediction.single` point
or a row in the prediction table may be outside the training and configured
sampling intervals. For Latin-hypercube extrapolation, set
`parameter_domain.range_source = "custom"` and give the desired sampling
interval. The command line and the saved `run_summary.json` identify points
outside the trained parameter domain.

`[energies] source` chooses the output energies:

- `"experimental"`: the energies of the experimental data file;
- `"grid"`: `grid = [first, last, step]` in MeV;
- `"selected"`: the list in `selected`.

`[comparison] exact_numerov = true` is the accuracy check. Exact Numerov
CCFULL runs at the same points and energies. `cross_sections.csv` then holds
`exact_numerov_cross_section_mb` and `rom_exact_relative_error`. The command
prints the median, 95th and 99th percentile, and maximum of
|σ_ROM/σ_exact − 1|, and `run_summary.json` stores them under `accuracy`.
`report_speedup = true` requires the check and reports exact Numerov wall
time divided by ROM wall time.

`[output] save_point_timings_csv = true`, the default, also writes
`point_timings.csv` with one row per parameter point: the six parameters,
`rom_seconds`, and, with the accuracy check, `exact_seconds` and `speedup`.
Both times are measured in the worker process, over all energies of the point.

## Commands

Validate the files without starting an expensive calculation:

```bash
python -m emulator.train --check
python -m emulator.run_rom --check
```

Run both offline stages, or either stage separately:

```bash
python -m emulator.train
python -m emulator.train --stage wave_pod
python -m emulator.train --stage cc_eim
```

After training, run the query in `emulator_rom.in`:

```bash
python -m emulator.run_rom
```

Each run creates one timestamped directory in `emulator/outputs/`. All
requested parameter points are written to the same `cross_sections.csv`;
`point_index` and the six parameter columns identify each curve.

With `model_id = "auto"`, the first CCFULL input record generates a directory
name such as `a16_z8__a154_z62`. To keep several trained variants of one
reaction, for example different ranks, set a distinct `model_id` for each.

`[model]` of `emulator_rom.in` names the model to use: its `model_id` and the
`training_root` it was trained into. The trained parameter box, energy range,
and experimental file come from that model's manifest, so switching models
needs no edit of `emulator_train.in`, and `run_rom --check` needs the trained
model.

## Change the reaction system

No Python source edit is required. In `emulator_train.in`:

1. replace the first nine records with the new CCFULL input;
2. select its matching `data.experimental_file`;
3. set the six training bounds and the energy range;
4. leave `model_id = "auto"`, or assign a distinct variant name.

After training, set `model_id` in `emulator_rom.in` to the new model, and
update its custom ranges and prediction point.
Training artifacts and online output names follow the new model ID
automatically. Configuration loading rejects an experimental file whose
projectile and target do not match the CCFULL input.
