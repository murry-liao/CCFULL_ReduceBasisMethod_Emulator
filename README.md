# CCFULL-ROM

**A reduced-order emulator for coupled-channels heavy-ion fusion calculations**

CCFULL-ROM replaces repeated runs of a CCFULL-type coupled-channels
calculation with a small projected model that is built once. For a rotational
target it returns the fusion cross section σ_fus as a function of:

- the deformations β₂, β₄, β₆;
- the Woods–Saxon depth, radius, and diffuseness V₀, r₀, a₀;
- the centre-of-mass energy E.

One reduced basis serves every parameter point, energy, and partial wave.
σ_fus can therefore be evaluated at any energy of the trained range: at the
measured energies for χ² scans, or on a fine grid for the barrier
distribution. χ² scans, optimisation, and Bayesian inference then no longer
need a coupled-channels run for every parameter point.

The emulator reduces the discretised coupled equations themselves rather than
fitting the cross section:

- **Basis.** Wave functions of many parameter points, energies, and partial
  waves are compressed into one basis by proper orthogonal decomposition
  (POD).
- **Projection.** The Numerov equations are projected onto the basis with the
  Galerkin method.
- **Coupling potential.** The non-affine coupling potential is made affine in
  the parameters by a multipole-wise empirical interpolation method (CC-EIM).

The cross section follows from the reduced wave solution and the original
CCFULL definitions.

## How it works

The numbers below belong to the `16O + 154Sm` example in this repository.

```text
 Offline (once per training box)                Online (per parameter point θ)
 ---------------------------------              ---------------------------------
 1. exact Numerov snapshots: 4000 points        1. potential at the EIM nodes only
    of a Latin hypercube in                        (72 radial points)
    (β2, β4, β6, V0, r0, a0, E, J)              2. EIM coefficients (120 numbers)
 2. one POD basis for all snapshots, rank 300   3. for every energy E and partial
 3. multipole CC-EIM of the coupling field         wave J: assemble the 300x300
    13 multipoles, 120 modes                       reduced matrix, solve -> P_J
 4. Galerkin projection of all operators        4. sigma_fus = (pi hbar^2 / 2 mu E)
                                                      * sum_J (2J+1) P_J
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

## Installation

**Requirements**

- Python ≥ 3.11, NumPy ≥ 2.0, SciPy ≥ 1.13
- `gfortran` and BLAS/LAPACK that link as `-lblas -llapack`. The Fortran solver
  is compiled automatically on first use into `physical_model/build/`.
- macOS or Linux (POSIX file locking is used by the build step).

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
python -m emulator.train --check       # print how each line of emulator/emulator_train.in is read
python -u -m emulator.train            # stage 1 (snapshots, basis), then stage 2 (CC-EIM, projection)

python -m emulator.run_rom --check     # print how each line of emulator/emulator_rom.in is read
python -m emulator.run_rom             # emulate and check against exact CCFULL; results in emulator/outputs/
```

Both input files are written like a CCFULL input: one item, or a few
comma-separated numbers, per line, with every line explained at the end of the
file. See [`emulator/README.md`](emulator/README.md).

**Analysis notebook.** `emulator/parameter_analysis.ipynb` reads a run
directory from `emulator/outputs/`. It finds the best parameter point by χ²
against experiment, and plots the excitation function and the barrier
distribution D(E) = d²(Eσ)/dE² at that point. Set `RUN` at the top of the
notebook to your run.

```bash
python -m pip install -e ".[notebook]"
python -m jupyter lab emulator/parameter_analysis.ipynb
```

## Repository layout

```text
.
├── emulator/
│   ├── train.py                     reads emulator_train.in; snapshots, POD, CC-EIM, projection
│   ├── run_rom.py                   reads emulator_rom.in; predictions and the accuracy check
│   ├── emulator_train.in            training input: CCFULL records + training lines, explained at the end
│   ├── emulator_rom.in              prediction input, explained at the end
│   ├── rom.py                       Galerkin projection, the online solve, trained-model manifests
│   ├── eim.py                       multipole CC-EIM of the coupling field
│   ├── parameter_analysis.ipynb     best parameters, excitation function, barrier distribution
│   ├── training_data/               generated models (not tracked)
│   └── outputs/                     generated query results (not tracked)
├── physical_model/
│   ├── ccfull.py                    compiles and calls CCFULL
│   ├── inputs.py                    CCFULL records, line-based inputs, experimental data files
│   ├── solvers/standard_numerov/    Fortran 90 CCFULL with standard Numerov, and its two drivers
│   └── build/                       generated binaries (not tracked)
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
