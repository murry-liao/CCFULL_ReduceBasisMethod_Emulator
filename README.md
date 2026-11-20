# CCFULL-ROM

**A reduced-order emulator for coupled-channels heavy-ion fusion calculations**

CCFULL-ROM replaces repeated runs of a CCFULL-type coupled-channels
calculation with a small projected model that is built once. For a rotational
target it returns the fusion cross section σ_fus as a function of:

- the deformations β₂, β₄, β₆;
- the Woods–Saxon depth, radius, and diffuseness V₀, r₀, a₀;
- the centre-of-mass energy E.

Two reduced models are built:

- **Global model.** One reduced basis serves every parameter point, energy,
  and partial wave, so σ_fus can be evaluated at any energy of the trained
  range
- **Local model.** Trained separately from its own CCFULL snapshots, every
  (E, J) of chosen energies (for example the measured ones) gets its own small
  basis. It is several times faster.

χ² scans, optimisation, and Bayesian inference then no longer need a
coupled-channels run for every parameter point.

The emulator reduces the discretised coupled equations themselves rather than
fitting the cross section:

- **Basis.** Wave functions of many parameter points, energies, and partial
  waves are compressed into one basis by proper orthogonal decomposition
  (POD).
- **Projection.** The Numerov equations are projected onto the basis with the
  Galerkin method.
- **Coupling potential.** The non-affine coupling potential is made affine in
  the parameters by an empirical interpolation method that treats all
  independent elements of the coupling matrix together.

The cross section follows from the reduced wave solution and the original
CCFULL definitions.

## How it works

The numbers below belong to the `34S + 168Er` example in this repository (seven
channels).

```text
 Offline (once per training box)                Online (per parameter point θ)
 ---------------------------------              ---------------------------------
 1. exact Numerov snapshots: 4000 points        1. potential at the EIM nodes only
    of a Latin hypercube in                        (42 radial points)
    (β2, β4, β6, V0, r0, a0, E, J)              2. EIM coefficients (48 numbers)
 2. one POD basis for all snapshots, rank 300   3. for every energy E and partial
 3. joint CC-EIM of the coupling field:            wave J: assemble the reduced
    28 matrix elements, 48 modes                   matrix (300x300 global, 80x80
 4. Galerkin projection of all operators           local), solve -> P_J
 5. local model (optional, separate): for       4. sigma_fus = (pi hbar^2 / 2 mu E)
    every (E, J) of the measured energies, a          * sum_J (2J+1) P_J
    rank-80 POD basis of 300 exact snapshots
```

The reduced operator of the block (E, J) is

```text
A(θ, E, J) = P_d2 − t [P_coul − E·P_one + (ħ²/2μ) J(J+1) P_cent]
             + A_static + Σ_m c_m(θ) A_m + boundary rows at r_min and r_max
```

- t = dr²·(2μ/ħ²)/12, and c_m(θ) are the CC-EIM coefficients.
- E and J(J+1) enter linearly, so every operator is projected offline.
- Online, only the two boundary rows are built for each (θ, E, J). The outer
  Coulomb map depends on (E, J) only and is cached per energy.
- The local model projects the same equations of one (E, J) onto the basis of
  that block, so E and J are fixed offline and only the CC-EIM sum and the
  inner boundary rows remain to be built online.


## Installation

**Requirements**

- Python ≥ 3.11, NumPy ≥ 2.0, SciPy ≥ 1.13
- `gfortran`, to build the two CCFULL programs: run `make` in
  `physical_model/ccfull_clean` and in `physical_model/ccfull_interface`. On
  Apple silicon use a gfortran that makes arm64 code, e.g.
  `make FC=/opt/homebrew/bin/gfortran`.
- emcee and pandas for the calibration (`emulator/mcmc.py`).

Developed and tested on macOS 13 (Apple silicon) with Python 3.12 and GNU
Fortran 15.1. Linux should work but has not been tested.

**Install**

```bash
conda create -n ccfull-rom python=3.12 -y     # or any virtual environment
conda activate ccfull-rom
python -m pip install -e .
```

## Quick start

Run all commands from the repository root. The repository ships without trained
models, so the first step is training.

```bash
make -C physical_model/ccfull_clean       # the exact model (add FC=... on Apple silicon)
make -C physical_model/ccfull_interface   # the same program, which also writes the wave functions

python -m emulator.train_local --check    # print how each line of emulator/emulator_train.in is read
python -u -m emulator.train_local         # the local model
python -u -m emulator.train_global        # the global model

python -m emulator.run_rom --check        # print how each line of emulator/emulator_rom.in is read
python -m emulator.run_rom                # emulate and check against exact CCFULL; results in emulator/outputs/

python -m emulator.mcmc --config emulator/emulator_mcmc_34S_168Er.in    # Bayesian calibration
```

## Repository layout

```text
.
├── emulator/
│   ├── train_local.py               reads emulator_train.in; the local model
│   ├── train_global.py              reads emulator_train.in; the global model
│   ├── run_rom.py                   reads emulator_rom.in; predictions and the accuracy check
│   ├── ccfull_run.py                writes ccfull.inp, runs CCFULL, reads its output and the data files
│   ├── coupling.py                  CCFULL's coupling matrix in Python (rotational target)
│   ├── rom.py                       the reduced model, and the CC-EIM fit and projection of the training
│   ├── emulator_*.in                inputs, each explained at its end
│   ├── parameter_analysis.ipynb     best parameters, excitation function, barrier distribution
│   ├── training_data/               generated models (not tracked)
│   └── outputs/                     generated query results (not tracked)
├── physical_model/
│   ├── ccfull_clean/                the exact model: F90 CCFULL, standard Numerov (see its Readme.dat)
│   └── ccfull_interface/            the same program with an output interface for the emulator
├── data/
│   ├── experimental/                measured fusion excitation functions
│   └── input_templates/             example batch-query table
└── pyproject.toml
```



## Related work

1. K. Hagino, N. Rowley and A. T. Kruppa, *A program for coupled-channel
   calculations with all order couplings for heavy-ion fusion reactions*,
   Comput. Phys. Commun. **123**, 143 (1999). The CCFULL code.
2. Z. Liao, K. Hagino, L. Zhu, S. Yoshida and K. Uzawa, *Determination of
   nuclear deformations with an emulator for sub-barrier fusion reactions*,
   Phys. Rev. C **113**, 044611 (2026). Eigenvector-continuation emulator for
   CCFULL fusion calculations.
3. M. Catacora-Rios, K. Beyer, J. Lei and F. M. Nunes, *Emulation of
   coupled-channels scattering for strongly deformed nuclei*,
   arXiv:2610.03910. Reduced-basis emulator with a multipole coupled-channels
   EIM.
4. M. Catacora-Rios, K. Beyer, P. Giuliani, K. Godbey, R. J. Furnstahl and
   F. M. Nunes, *Wavefunction-based emulation of coupled-channels scattering
   with non-affinely parametrized interactions*, arXiv:2512.08097.
5. D. Odell et al., *ROSE: A reduced-order scattering emulator for optical
   models*, Phys. Rev. C **109**, 044612 (2024).
6. J. Lei, *Reduced basis emulator for elastic scattering in
   continuum-discretized coupled-channels calculations*, arXiv:2512.17687.
7. J. Liu, J. Lei and Z. Ren, *A complex scaling method for efficient and
   accurate scattering emulation in nuclear reactions*, Phys. Lett. B **858**,
   139070 (2024).
8. J. A. Melendez et al., *Model reduction methods for nuclear emulators*,
   J. Phys. G **49**, 102001 (2022). Review of the reduced-basis framework.

Background on the β₆ deformation in fusion: T. Rumin, K. Hagino and
N. Takigawa, Phys. Rev. C **61**, 014605 (2000); C. R. Morton et al.,
Phys. Rev. C **64**, 034604 (2001).

If you use this code in published work, please cite CCFULL (1) and the emulator
papers (2) that are relevant to your use.
