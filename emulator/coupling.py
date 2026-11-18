"""CCFULL's potentials in Python, for the emulator's online stage.

The reduced model is rebuilt at every parameter point, so it needs CCFULL's coupling
matrix at a few radii (the CC-EIM nodes) for any (beta2, beta4, beta6, V0, r0, a0), and
CCFULL's rule for where the J loop ends.  This file writes out what CCFULL computes there,
for a rotational target and an inert projectile (line 2 of ccfull.inp: IVIBROTP = -1,
IVIBROTT = 1).  The function names follow CCFULL.  The interface program writes CCFULL's
own matrices to coupling.bin, against which this file is checked.
"""

import math

import numpy as np

from emulator.ccfull_run import numbers


HBAR = 197.329          # MeV fm, Hbar in CCFULL
PI = 3.141592653        # Pi in CCFULL, with its digits
AMU = 938.0             # MeV, Amu in CCFULL


def reaction(records):
    """The constants of the reaction, from the nine CCFULL records."""
    apro, zpro, atar, ztar = (int(value) for value in numbers(records[0]))
    r0p, ivibrotp, r0t, ivibrott = numbers(records[1])
    if ivibrotp != -1 or ivibrott != 1:
        raise ValueError("coupling.py writes out a rotational target with an inert projectile only")
    target = numbers(records[2])                      # E2T, beta2, beta4, beta6, NrotT
    return {"apro": apro, "zpro": zpro, "atar": atar, "ztar": ztar,
            "rtar": r0t * atar ** (1.0 / 3.0),       # Rtar
            "e2": target[0], "channels": int(target[-1]) + 1,
            "reduced_mass": apro * atar / (apro + atar) * AMU}


def clebsch_gordan(j1, m1, j2, m2, j3, m3):
    """<j1 m1 j2 m2 | j3 m3> for integer angular momenta (CG in CCFULL, Racah's formula)."""
    if m1 + m2 != m3 or j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0
    f = lambda n: float(math.factorial(n))
    total = 0.0
    for n in range(max(j1 + j2 - j3, j1 - m1, j2 + m2) + 1):
        k = (j1 + j2 - j3 - n, j3 - j2 + m1 + n, j3 - j1 - m2 + n, j1 - m1 - n, j2 + m2 - n)
        if min(k) < 0:
            continue
        total += (-1.0) ** n / (f(n) * f(k[0]) * f(k[1]) * f(k[2]) * f(k[3]) * f(k[4]))
    total *= math.sqrt(f(j1 + m1) * f(j1 - m1)) * math.sqrt(f(j2 + m2) * f(j2 - m2))
    total *= math.sqrt(f(j3 + m3) * f(j3 - m3))
    triangle = math.sqrt(f(j1 + j2 - j3) * f(j3 + j1 - j2) * f(j2 + j3 - j1) / f(j1 + j2 + j3 + 1))
    return total * math.sqrt(2.0 * j3 + 1.0) * triangle


def rotor_geometry(channels, multipole):
    """sqrt((2I+1)(2 lambda+1)(2I'+1)/4 pi) <I 0 lambda 0|I' 0>^2 / (2I'+1) between the rotor
    states I = 0, 2, 4, ... (coupled_matrix0 and coupled_matrix take I <= I' and copy it)."""
    geometry = np.zeros((channels, channels))
    for it in range(channels):
        for jt in range(it, channels):
            spin, spin2 = 2 * it, 2 * jt
            geometry[it, jt] = (math.sqrt((2 * spin + 1) * (2.0 * multipole + 1) * (2 * spin2 + 1) / (4.0 * PI))
                                * clebsch_gordan(spin, 0, multipole, 0, spin2, 0) ** 2 / (2.0 * spin2 + 1))
            geometry[jt, it] = geometry[it, jt]
    return geometry


def coupling_matrix(radii, theta, system):
    """Everything in CCFULL's channel potential that depends on the parameter point, at radii:
    (radius, channel, channel).  It equals CPOT + Vn(r) - the rotational energies.

    Nuclear part (coupled_matrix0, coupled_matrix): the deformation operator
    R_T (beta2 Y20 + beta4 Y40 + beta6 Y60) in the rotor states is diagonalized, and the
    Woods-Saxon VnCC is taken at its eigenvalues, to all orders.
    Coulomb part (Fct2, Fct4, Fct6): multipoles 2, 4, 6, with the beta2^2 terms of CCFULL.
    """
    beta2, beta4, beta6, v0, r0, a0 = (float(value) for value in theta)
    r = np.atleast_1d(np.asarray(radii, dtype=float))
    rt, channels = system["rtar"], system["channels"]
    g2, g4, g6 = (rotor_geometry(channels, multipole) for multipole in (2, 4, 6))

    deformation = rt * beta2 * g2 + rt * beta4 * g4 + rt * beta6 * g6
    eigenvalues, eigenvectors = np.linalg.eigh(deformation)
    r12 = r0 * (system["apro"] ** (1.0 / 3.0) + system["atar"] ** (1.0 / 3.0))
    vncc = -v0 / (1.0 + np.exp((r[:, None] - r12 - eigenvalues[None, :]) / a0))   # (radius, eigenchannel)
    nuclear = np.einsum("ik,rk,jk->rij", eigenvectors, vncc, eigenvectors)

    def multipole_form(lam):                      # Fct2, Fct4, Fct6 without their beta factors
        zz = system["zpro"] * system["ztar"]
        outside = 3.0 / (2.0 * lam + 1.0) * zz / 137.0 * HBAR / r * (rt / r) ** lam
        inside = 3.0 / (2.0 * lam + 1.0) * zz / 137.0 * HBAR / rt * (r / rt) ** lam
        return np.where(r > rt, outside, inside)[:, None, None]

    coulomb = (g2 * (beta2 + 2.0 * math.sqrt(5.0 / PI) * beta2 ** 2 / 7.0) * multipole_form(2)
               + g4 * (beta4 + 9.0 * beta2 ** 2 / 7.0 / math.sqrt(PI)) * multipole_form(4)
               + g6 * beta6 * multipole_form(6))
    return nuclear + coulomb


# ----------------------------------------------------------------- where CCFULL's J loop ends

def bare_potential(r, angular_momentum, v0, r0, a0, system):
    """V(r, L) of CCFULL: Woods-Saxon + point Coulomb + centrifugal."""
    r12 = r0 * (system["apro"] ** (1.0 / 3.0) + system["atar"] ** (1.0 / 3.0))
    h2m = HBAR ** 2 / (2 * system["reduced_mass"])
    return (-v0 / (1.0 + math.exp((r - r12) / a0)) + system["zpro"] * system["ztar"] * HBAR / 137.0 / r
            + h2m * angular_momentum * (angular_momentum + 1.0) / r ** 2)


def bare_slope(r, angular_momentum, v0, r0, a0, system):
    """dV/dr of CCFULL (dV_dr)."""
    r12 = r0 * (system["apro"] ** (1.0 / 3.0) + system["atar"] ** (1.0 / 3.0))
    growth = math.exp((r - r12) / a0)
    return (v0 / a0 * growth / (1.0 + growth) ** 2 - system["zpro"] * system["ztar"] * HBAR / 137.0 / r ** 2
            - 2.0 * angular_momentum * (angular_momentum + 1.0) * HBAR ** 2 / (2.0 * system["reduced_mass"] * r ** 3))


def pocket_bottom(angular_momentum, v0, r0, a0, system):
    """The bottom of the pocket of V(r, L), as PotShape finds it; None without a barrier."""
    slope = lambda r: bare_slope(r, angular_momentum, v0, r0, a0, system)
    # the barrier: step inwards from 50.5 fm until dV/dr changes sign, then bisect
    r = 50.5
    u0 = slope(r)
    while True:
        r -= 1.0
        if r < 0.0:
            return None
        u1 = slope(r)
        if u0 * u1 <= 0.0:
            break
        u0 = u1
    inner, outer = r, r + 1.0
    for _ in range(int(math.log10(abs(inner - outer) / 1.0e-6) / math.log10(2.0) + 0.5)):
        r = (outer + inner) / 2.0
        if u0 * slope(r) < 0.0:
            inner = r
        else:
            outer = r
    barrier = r
    # the pocket: bisect between 0.5 fm and half a fermi inside the barrier
    outer, inner = barrier - 0.5, 0.5
    u0 = slope(outer)
    for _ in range(int(math.log10(abs(barrier - outer) / 1.0e-6) / math.log10(2.0) + 0.5)):
        r = (outer + inner) / 2.0
        if u0 * slope(r) < 0.0:
            inner = r
        else:
            outer = r
    return r


def partial_wave_count(energy, v0, r0, a0, system, lmax=200):
    """How many partial waves J = 0, 1, ... have a pocket below E: CCFULL's J loop stops at
    the first one without (it also stops earlier, once the cross section converges)."""
    count = 0
    for angular_momentum in range(lmax + 1):
        bottom = pocket_bottom(angular_momentum, v0, r0, a0, system)
        if bottom is None or bottom < 0.0 or bare_potential(bottom, angular_momentum, v0, r0, a0, system) > energy:
            break
        count += 1
    return count
