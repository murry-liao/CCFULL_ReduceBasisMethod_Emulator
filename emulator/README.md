# Emulator

A reduced model of CCFULL, trained once and then evaluated at many parameter
points θ = (β2, β4, β6, V0, r0, a0).

## Programs

Each program reads one input file line by line, as CCFULL reads `ccfull.inp`.
Every line is explained at the end of its file, and `--check` prints how each
line was read and stops.

| Command | Input | What it does |
|---|---|---|
| `python -m emulator.train_local` | `emulator_train.in` | local model: one small basis per (E, J) of chosen energies |
| `python -m emulator.train_global` | `emulator_train.in` | global model: one basis for every (E, J) of an energy range |
| `python -m emulator.run_rom` | `emulator_rom.in` | the model at parameter points, checked against exact CCFULL |
| `python -m emulator.mcmc --config <file>` | `emulator_mcmc_*.in` | Bayesian calibration of the six parameters |

Each program is written from top to bottom in numbered steps that match its
docstring. They share four modules:

- `ccfull_run.py` talks to CCFULL and reads files. It writes `ccfull.inp`,
  runs `../physical_model/ccfull_clean` (exact cross sections) or
  `ccfull_interface` (wave functions), and reads their output. It also reads
  the emulator's input files and the experimental data.
- `coupling.py` holds CCFULL's coupling matrix in Python. The reduced model
  needs it at the CC-EIM nodes for every new θ. It covers a rotational target
  with an inert projectile, and it agrees with CCFULL's matrices to 5e-13 MeV.
  It also holds CCFULL's rule for the last partial wave (PotShape).
- `rom.py` is the reduced model: loading a trained model, the CC-EIM
  coefficients, the boundary conditions, the solves, and σ. Its last section
  holds what the training needs: the CC-EIM fit, the Galerkin projection, and
  the local blocks.
- The notebook `parameter_analysis.ipynb` reads the output files.

## The reduced model

CCFULL solves the coupled equations of every partial wave J at energy E on a
radial grid with the standard Numerov step. The reduced model writes the
discrete equations of all grid points as one linear system A u = b and
projects it onto a small basis W: (W^H A W) x = W^H b, with u = W x.

A depends on θ only through the coupling field. The CC-EIM interpolates all
independent elements of the coupling matrix together, which writes the field
as V(r; θ) = V_static(r) + Σ_m c_m(θ) B_m(r). Every piece of W^H A W is
therefore projected once in the training.

Online, a point costs:
- the coupling matrix at the CC-EIM node radii (`coupling.py`), which gives
  c_m(θ);
- for every (E, J), the inner boundary rows (CCFULL's Runge–Kutta start at
  r_min), the sum over the modes, and one small solve;
- then P_J and CCFULL's J sum.

There are two kinds of model:
- **local**: one basis per (E, J) of the training energies, rank 24–80. It is
  the fast one, for the measured energies of one data set.
- **global**: one basis for all (E, J), rank about 300. It works at any energy
  of its range. Its outer boundary is set for the energies of each run, from
  one run of `ccfull_interface` at the nominal point.

The model uses CCFULL's own conventions:
- **Energy.** E is held in single precision, as CCFULL declares it.
- **π.** CCFULL's digits of π.
- **Outer boundary.** The Coulomb functions of CCFULL's matching, which
  `ccfull_interface` writes for J = 0 and 1. The higher J follow by the
  standard recurrence, which agrees with CCFULL to 3e-12.

## Training

`emulator_train.in`:
- lines 1–9 are the CCFULL records of the reaction, and they define the
  nominal point;
- lines 10–19 set the name, the data, the box, and the sizes.

**The grid.** CCFULL puts r_min at the bottom of the pocket of the nominal
point. Every later run keeps this r_min (line 10 of `ccfull.inp`), so that all
snapshots share one radial grid.

**Partial waves.** At each energy they are J = 0 … n−1, where n is the
largest count of CCFULL's J loop over the corners of the (V0, r0, a0) box.
Every training point solves all of them (line 11 of `ccfull.inp`), even
beyond the end of its own J sum. Without this, the blocks of the highest J
would see only the points with a low barrier, and would give unusable P_J
elsewhere.

**Local model.**
- At each energy, `ccfull_interface` runs at the Latin-hypercube points and
  at the 64 corners of the box.
- A block is kept for every J that the J sum reaches at a point or a corner,
  plus 5 more.
- Each block is the POD of its 300 snapshots and the Numerov equations
  projected onto its leading vectors.
- One file is written per energy.

**Global model.** 4000 Latin-hypercube points in (θ, E, u). At each one,
`ccfull_interface` solves the partial wave J = ⌊u n⌋. Then follow the POD,
the CC-EIM, and the projection.

A trained model keeps copies of the CCFULL records (`ccfull.inp`) and of the
interface program (`ccfull_interface.f90`). `run_rom` and `mcmc` refuse it
once either has changed. A training stops if its model folder already exists.

## Files

`training_data/<name>/` holds one model; the name tells how it was trained,
for example `34S_168Er_Local_Channel6_Snap300_POD80_EIM40`:
- `manifest.json` — settings and diagnostics.
- `shared_interaction/` — the CC-EIM.
- Local model: `grid.npz` and one `energy_NNN.npz` per energy. These hold the
  blocks: `a0`, `a_modes` (single precision, nearly the whole model), `x0`,
  `x1`, `rhs`, and the POD tail.
- Global model: `reduced_operators.npz`, `basis.npz`, `snapshots.npz`, and
  `pod_singular_values.csv`.

`outputs/<time>_<name>[_r<rank>]/` holds one run of `run_rom`:
`cross_sections.csv`, `point_timings.csv`, and `run_summary.json`.

`outputs/<time>_<name>[_r<rank>]_mcmc/` holds one calibration:
- the copied input, `progress.log`, and `mcmc_summary.json`;
- `scan.csv`, `best_fit.csv`, and `best_fit.npz`;
- the chains `pilot_chain.npz` and `chain.npz`;
- `emulator_error.csv` and `emulator_error_points.csv`;
- `posterior_samples.csv`, `posterior_check.csv`, and
  `posterior_predictive.csv`.

## Timing

`run_rom` times the reduced model per point, without loading or preparing the
model. It times exact CCFULL as `ccfull_clean` running once per measured
energy, so this time includes 41 program starts.

The ROM keeps BLAS at one thread, and every program runs in a single process.
