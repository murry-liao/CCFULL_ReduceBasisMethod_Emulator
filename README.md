# CCFULL-DBM-EC Emulator

Fortran research code for emulating coupled-channels calculations of sub-barrier heavy-ion fusion with eigenvector continuation (EC), a reduced-basis method (RBM). This CCFULL is solved by discrete basis method (DBM)

The implementation accompanies:

> Z. Liao, K. Hagino, L. Zhu, S. Yoshida, and K. Uzawa,  
> "Determination of nuclear deformations with an emulator for sub-barrier fusion reactions,"  
> *Physical Review C* **113**, 044611 (2026).  
> DOI: [10.1103/8ycb-lrqq](https://doi.org/10.1103/8ycb-lrqq)

The paper is also included in this repository as [`PhysRevC.113. 044611 (2026).pdf`](PhysRevC.113.%20044611%20%282026%29.pdf).

## What this code does

The program combines a CCFULL-like coupled-channels model with a reduced-basis emulator. It calculates fusion penetrabilities and cross sections while varying nuclear deformation parameters.

The current source snapshot is configured for the two-dimensional rotational-deformation study of the `16O + 186W` system in Sec. III D of the paper:

- emulator parameters: target quadrupole deformation `beta2T` and hexadecapole deformation `beta4T`;
- training grid: `beta2T = 0.25, 0.30, 0.35` and `beta4T = -0.025, 0.000, 0.025`;
- evaluation grid: `beta2T = 0.25 ... 0.35` and `beta4T = -0.05 ... 0.05`, both in steps of `0.01`;
- nine high-fidelity snapshots and 121 evaluation points;
- partial waves from `J = 0` to `200`;
- modified Numerov propagation and incoming-wave boundary conditions;
- direct and emulated results produced together for comparison.

## Calculation flow

For every center-of-mass energy and partial wave, the program follows this sequence:

1. Read the reaction, excitation, potential, energy-grid, and radial-grid settings from `ccfull.inp`.
2. Locate the Coulomb-barrier top and the pocket minimum. The pocket minimum defines the inner radial boundary.
3. Build the intrinsic-state and coupling matrices for vibrational or rotational excitations.
4. At each of the nine training points, construct the modified-Numerov discrete Hamiltonian and solve the full coupled-channels problem. The resulting scattering wave functions are stored as high-fidelity snapshots.
5. At each of the 121 evaluation points:
   - solve the full problem with `Numerov_Ec` to obtain the reference penetrability;
   - project the parameter-dependent problem onto the snapshot space with `Ec_emulator`;
   - solve the reduced system and obtain the emulated penetrability.
6. Sum `(2J + 1) P_J` over partial waves to obtain the fusion cross section.
7. Write the direct result, emulator result, relative error, and timing information.

The reduced system has dimension `NumEmulator + Nlevel`, rather than the full radial-grid dimension. Training is performed independently for each energy and partial wave, consistent with the procedure described in the paper.

## Requirements

- A Fortran compiler supporting the source's legacy extensions, such as GNU Fortran (`gfortran`)
- BLAS
- LAPACK

The code calls `ZGEMM`, `ZGEMV`, `ZGESV`, `ZGETRF`, and `ZGETRI`.

## Build

From the repository root:

```bash
gfortran -O3 -o ccfull_emulator ccfull_emulator_real.f90 -llapack -lblas
```

This command has been checked with GNU Fortran 15.1.0. The source compiles successfully, although strict warning flags report several implicit precision conversions and unused variables typical of this research-code snapshot.

For a diagnostic build:

```bash
gfortran -O0 -g -Wall -Wextra -fcheck=all \
  -o ccfull_emulator_debug ccfull_emulator_real.f90 -llapack -lblas
```

## Input

The executable always reads a file named `ccfull.inp` from the current working directory. The parser consumes one case only, using the following record order:

| Record | Values | Meaning |
| --- | --- | --- |
| 1 | `Apro, Zpro, Atar, Ztar` | Projectile and target mass and charge numbers. |
| 2 | `R0P, IVIBROTP, R0T, IVIBROTT` | Projectile/target radius parameters and excitation modes. Mode `-1` is inert, `0` is vibrational, and `1` is rotational. |
| 3-4 | Target excitation data | Two vibrational-mode records, or one rotational record followed by one placeholder record. |
| 5 | Projectile excitation data | Vibrational or rotational data, or a placeholder record for an inert projectile. |
| 6 | `Qtrans, Ftr` | Transfer-channel parameters used by the current parser. Any additional value on this record is ignored. |
| 7 | `V0, R0, A0` | Woods-Saxon depth, radius, and diffuseness parameters in MeV and fm. |
| 8 | `Emin, Emax, dE` | Center-of-mass energy grid in MeV. |
| 9 | `R_max, dr` | Outer radial boundary and mesh spacing in fm. |

For a rotational target, record 3 is:

```text
E2T, Beta2T, Beta4T, NrotT
```

For a vibrational target, records 3 and 4 are:

```text
OmegaT,  BetaT,  LambdaT,  NphononT
OmegaT2, BetaT2, LambdaT2, NphononT2
```

The active block supplied in `ccfull.inp` begins with:

```text
16,8,186,74
1.2,-1,1.06,1
0.,0.29,-0.03,6
1.66,0.11,2,0
6.13,0.733,3,1
0,0.,0.3
188,0.95,1.05
60.,80.,0.5
30,0.05
```

Only these first nine records are read. To run another case, move its complete block to the beginning of the file and ensure that the hard-coded emulator parameter mapping and grids in `Initialize_CCFull` and `Prameterize_H` match the intended physics case.

## Run

```bash
./ccfull_emulator
```

The default calculation is computationally demanding: it repeats training and direct validation for every energy, partial wave, and point in the 11 x 11 deformation grid.

## Output

For every deformation pair, the program appends data to a file named like:

```text
CCfull_sigma_16O+186W_Ec_Beta2T_0.290_Beta4T_-0.020.dat
```

Each row contains:

```text
E  P_exact(J=0)  P_emulator(J=0)  sigma_exact  sigma_emulator  relative_sigma_error
```

Cross sections are written in millibarns. Detailed per-partial-wave validation and timing data are written implicitly to `fort.999` with columns:

```text
E  J  parameter_index  P_exact  P_emulator  relative_P_error  exact_time  emulator_time
```

Output files are opened in append mode. Move or remove results from an earlier run before restarting if duplicate energy rows are not desired.

## Implementation notes and current limitations

- The deformation training and evaluation grids are hard-coded in `Initialize_CCFull`.
- `Prameterize_H` currently maps its two emulator parameters to `Beta2T` and `Beta4T`; changing to the vibrational `144Sm` case requires a source edit.
- The program uses fixed limits `Nlevelmax = 30` and `J <= 200`.
- Several module variables rely on compiler/static initialization in the current snapshot, including values associated with optional second-mode and transfer channels. Explicit initialization is recommended before extending those paths.
- The transfer-input record and the `Ntrans`/`iq` variables should be reviewed before enabling transfer coupling.
- The direct reference solve remains active at every evaluation point so that emulator errors can be measured. For production-only emulation, that validation path can be disabled after verification.
- The experimental-data comparison, chi-squared surface, posterior distribution, and parameter extraction shown in the paper are not implemented in this Fortran file; they must be performed as post-processing of the generated cross sections.
- The repository is a reproducibility-oriented research snapshot, not a general-purpose replacement for the original CCFULL distribution.

## Citation

If this code contributes to published work, please cite:

```bibtex
@article{Liao2026NuclearDeformationEmulator,
  author  = {Liao, Zehong and Hagino, K. and Zhu, Long and Yoshida, S. and Uzawa, K.},
  title   = {Determination of nuclear deformations with an emulator for sub-barrier fusion reactions},
  journal = {Physical Review C},
  volume  = {113},
  pages   = {044611},
  year    = {2026},
  doi     = {10.1103/8ycb-lrqq}
}
```

## Contact

- Zehong Liao: `liaozh26@mail2.sysu.edu.cn`
