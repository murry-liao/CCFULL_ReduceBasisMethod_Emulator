# Experimental data and query templates

This directory contains inputs of the physical model and emulator. It contains
no executable code.

- `experimental/4-final-ngt4/`: curated fusion excitation-function data.
- `input_templates/`: the batch-query parameter table.

The active experimental library contains 410 excitation-function files for
369 reaction systems and 7,311 measured points. Each row stores
`E_cm, sigma_fusion, uncertainty_lower, uncertainty_upper, E_lab`. The default
emulator dataset is
`experimental/4-final-ngt4/fusion-exp-8.16.62.154-tb-1995-611.txt`.

`input_templates/rom_parameter_points.csv` is the six-column batch-query
template used by `emulator_rom.in`. The CCFULL input of a trained system is
the first nine records of `../emulator/emulator_train.in`.

Generated wave bases, emulator models, analyses, and online results belong
under `../emulator/`.
