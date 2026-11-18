# Physical model

The exact model of the emulator is the F90 CCFULL program of
`CCfull_Real_fortran90_numerov/ccfull.f90`, in two versions:

- **`ccfull_clean/`** — the exact model. It differs from the F90 program in
  three places, each explained in its `Readme.dat`:
  - the standard Numerov step instead of Hagino's modified one, so that the
    discrete equations are linear in the potential;
  - a stabilization every 0.5 fm besides the one at the barrier;
  - an optional line 10 of `ccfull.inp` that fixes r_min.

  The unused snapshot routines of the original file are removed. It reads
  `ccfull.inp` and writes `sigma.dat` and `Output.dat`, like the original.
- **`ccfull_interface/`** — the same program with an output interface. Its
  `sigma.dat` is byte-identical to that of `ccfull_clean`. Every added line is
  marked `! interface` (`diff ccfull_clean/ccfull.f90 ccfull_interface/ccfull.f90`).
  It also writes:
  - the wave function and P_J of every partial wave it solves
    (`waves_info.txt`, `waves.bin`);
  - the coupling matrices on the grid (`coupling.bin`);
  - the Coulomb functions of its matching (`matching.txt`).

  An optional line 11 of `ccfull.inp` makes it go on solving partial waves
  after the J sum has ended, up to that J, for their wave functions only.

The emulator (`../emulator`) uses these programs by writing `ccfull.inp`:

- lines 3 and 7: the six parameters;
- line 8: the energy;
- line 10: the r_min of its grid;
- line 11: the last J.

It then reads the output files.

Build each program with `make` in its folder. On Apple silicon, use a
gfortran that makes arm64 code, for example
`make FC=/opt/homebrew/bin/gfortran`; check it with `gfortran -dumpmachine`.

The comparison of the standard and modified Numerov steps, and of one
stabilization against one every 0.5 fm, is in
`../../Plot/CalData/34S+168Er_numerov/README.md`.
