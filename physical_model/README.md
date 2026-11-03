# Physical model

- `solvers/standard_numerov/ccfull.f90` is the high-fidelity model and
  wave-snapshot source. Trained models record its SHA-256, so editing it
  requires retraining.
- `solvers/standard_numerov/snapshot_driver.f90` and `library_interface.f90`
  replace the main program of `ccfull.f90` in its two builds: a program that
  solves one (θ, E, J) and writes the wave function, and a shared library that
  stays loaded in Python for coupling potentials, multipole form factors,
  Coulomb functions, partial-wave counts, and exact penetrabilities.
- `ccfull.py` compiles both with the gfortran that generates code for the
  running Python's architecture, reuses a build while its source is
  unchanged, and calls them. Numerov always propagates in `complex128`; only
  the exported training snapshots are `complex64`.
- `inputs.py` reads the nine CCFULL records and the experimental data files.
- The CCFULL records used for training are the first nine records of
  `../emulator/configs/emulator_train.in`.
- `../data/experimental/4-final-ngt4/` is the active experimental library: 410
  curated fusion datasets covering 369 reaction systems.
- `build/` is generated and may be recreated from source.
