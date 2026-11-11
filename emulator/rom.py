"""The global-basis ROM: one reduced basis for every parameter point, energy, and partial wave.

The Numerov equations of the full radial grid are projected onto one basis Q
of rank R.  The reduced operator of the block (E, J) is

    A = P_d2 - t [P_coul - E P_one + (hbar^2 / 2 mu) J(J+1) P_cent]
        + A_static + sum_m c_m(theta) A_m + boundary rows at r_min and r_max,

with t = dr^2 (2 mu / hbar^2) / 12 and c_m(theta) the CC-EIM coefficients.
E and J(J+1) enter linearly, so every operator is projected offline.  Online,
the CC-EIM coefficients are computed once per parameter point, and each
(E, J) needs its two boundary rows and one R x R solve.  The outer Coulomb
map depends on (E, J) only and is cached per energy.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from emulator.eim import MultipoleEIM, field_radii, load_eim
from physical_model import ccfull
from physical_model.ccfull import HBARC
from physical_model.inputs import parse_input, records_sha256


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PARAMETER_NAMES = ("beta2", "beta4", "beta6", "v0", "r0", "a0")
MODEL_FORMAT = "global_projection_cc_eim_v1"
# The arrays of reduced_operators.npz besides radius, reduced_mass, and zz.
ARRAYS = ("x0", "x1", "xn", "xnm2", "p_d2", "p_one", "p_coul", "p_cent", "a_static", "a_modes")


def project_path(path: str | Path) -> Path:
    """A path of an input file or a manifest: relative to the project root, or absolute."""
    return (PROJECT_ROOT / Path(path).expanduser()).resolve()


def project_relative(path: str | Path) -> str:
    """``path`` as a manifest records it: relative to the project root, or absolute outside it."""
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def latin_hypercube(bounds, count: int, seed: int) -> np.ndarray:
    """``count`` Latin-hypercube points inside ``bounds`` (n, 2); equal bounds fix a coordinate."""
    limits = np.asarray(bounds, dtype=float)
    points = np.repeat(limits[:, 0][None, :], count, axis=0)
    active = limits[:, 1] > limits[:, 0]
    if np.any(active):
        unit = qmc.LatinHypercube(int(np.sum(active)), seed=seed).random(count)
        points[:, active] = qmc.scale(unit, limits[active, 0], limits[active, 1])
    return points


def read_manifest(path: str | Path) -> dict:
    """A trained model's manifest; refused if ccfull.f90 or the nine CCFULL records changed since."""
    manifest = json.loads(Path(path).read_text())
    settings = manifest["settings"]
    if settings["solver_sha256"] != ccfull.source_sha256():
        raise ValueError(f"{path} was trained with another ccfull.f90; train the model again")
    if settings["physical_input_sha256"] != records_sha256(project_path(settings["input_file"])):
        raise ValueError(f"{path} was trained with other CCFULL records than {settings['input_file']} "
                         "now holds; restore them or train the model again")
    return manifest


def fusion_cross_section(penetrabilities, waves, energy: float, reduced_mass: float) -> float:
    """Fusion cross section in mb.

    The J sum stops once an increment is below 1e-4 of the running sum.  The
    ROM and the exact reference both call this function, so they share the rule.
    """
    prefactor = np.pi * HBARC**2 / (2.0 * reduced_mass * energy) * 10.0
    total = 0.0
    for angular_momentum, value in zip(waves, penetrabilities):
        previous = total
        total += prefactor * (2 * int(angular_momentum) + 1) * float(value)
        if previous > 0.0 and total - previous < previous * 1.0e-4:
            break
    return total


class Reaction:
    """The CCFULL input of a model, its parameter box and energy range, and its CCFULL library."""

    def __init__(self, input_file: str | Path, bounds, energy_range):
        self.config = parse_input(input_file)
        self.bounds = np.asarray(bounds, dtype=float)       # (6, 2), in PARAMETER_NAMES order
        self.energy_range = tuple(map(float, energy_range))
        self._waves: dict[float, np.ndarray] = {}

    @classmethod
    def from_settings(cls, settings: dict) -> "Reaction":
        """The reaction recorded in a trained manifest."""
        bounds = dict(zip(settings["parameter_names"], settings["parameter_bounds"]))
        return cls(project_path(settings["input_file"]), [bounds[name] for name in PARAMETER_NAMES],
                   settings["energy_range_mev"])

    @property
    def library(self) -> ccfull.Library:
        return ccfull.library(self.config)

    @property
    def nominal(self) -> np.ndarray:
        """The six parameters written in the CCFULL input."""
        target = self.config["target1"]
        return np.array([target[1], target[2], target[3],
                         self.config["v0"], self.config["r0"], self.config["a0"]], dtype=float)

    def partial_waves(self, energy: float) -> np.ndarray:
        """J = 0 .. n-1, n the largest CCFULL partial-wave count over the (V0, r0, a0) box corners."""
        key = round(float(energy), 10)
        if key not in self._waves:
            corners = itertools.product(*self.bounds[3:])
            self._waves[key] = np.arange(max(self.library.partial_wave_count(energy, v0, r0, a0)
                                             for v0, r0, a0 in corners))
        return self._waves[key]

    def exact_cross_sections(self, point, energies, reduced_mass: float) -> np.ndarray:
        """Exact Numerov CCFULL cross sections with the partial waves and J sum of the ROM."""
        values = []
        for energy in energies:
            waves = self.partial_waves(energy)
            penetrabilities = self.library.penetrabilities(point, float(energy), waves)
            values.append(fusion_cross_section(penetrabilities, waves, float(energy), reduced_mass))
        return np.array(values)


class Projector:
    """Offline Galerkin projection of the Numerov operator pieces onto a basis."""

    def __init__(self, basis: np.ndarray, radius: np.ndarray, reduced_mass: float, zz: int,
                 channels: int):
        self.radius = np.asarray(radius, dtype=float)
        self.x = np.asarray(basis, dtype=np.complex128).reshape(self.radius.size, int(channels), -1)
        dr = float(self.radius[1] - self.radius[0])
        self.t = dr**2 * (2.0 * float(reduced_mass) / HBARC**2) / 12.0
        self.zz = int(zz)

    def project(self, residual: np.ndarray) -> np.ndarray:
        """Q^H residual."""
        return np.tensordot(self.x.conj(), residual, axes=([0, 1], [0, 1]))

    def weighted(self, f: np.ndarray, y: np.ndarray | None = None) -> np.ndarray:
        """Q^H of the Numerov-weighted product f[i+1] y[i+1] + 10 f[i] y[i] + f[i-1] y[i-1].

        Interior rows only.  ``f`` is a scalar or a channel matrix per radius;
        ``y`` is the basis unless given.
        """
        y = self.x if y is None else y
        residual = np.zeros_like(y)
        if np.ndim(f) == 1:
            residual[1:-1] = (f[2:, None, None] * y[2:] + 10.0 * f[1:-1, None, None] * y[1:-1]
                              + f[:-2, None, None] * y[:-2])
        else:
            residual[1:-1] = f[2:] @ y[2:] + 10.0 * (f[1:-1] @ y[1:-1]) + f[:-2] @ y[:-2]
        return self.project(residual)

    def operators(self, eim: MultipoleEIM) -> dict:
        """Every reduced operator, and the basis rows of the boundary conditions."""
        x, n = self.x, self.radius.size
        d2 = np.zeros_like(x)
        d2[1:-1] = x[2:] - 2.0 * x[1:-1] + x[:-2]
        operators = {
            # Basis rows at r_min, r_min + dr, r_max, and r_max - 2 dr.
            "x0": x[0], "x1": x[1], "xn": x[-1], "xnm2": x[-3],
            "p_d2": self.project(d2),
            "p_one": self.weighted(np.ones(n)),
            # Coulomb Zp Zt e^2 / r with e^2 = hbar c / 137, as in CCFULL.
            "p_coul": self.weighted(self.zz * HBARC / 137.0 / self.radius),
            "p_cent": self.weighted(1.0 / self.radius**2),
            "a_static": -self.t * self.weighted(np.asarray(eim.static_field[:n], dtype=float)),
        }
        # One operator per CC-EIM mode: its radial function times the channel matrix of its multipole.
        modes = []
        for geometry, item in zip(eim.geometries, eim.radial):
            gx = np.einsum("cd,rdk->rck", geometry, x)
            for k in range(item.basis.shape[1]):
                modes.append(-self.t * self.weighted(item.basis[:n, k], gx))
        operators["a_modes"] = np.stack(modes)
        return operators

    def field_operator(self, field: np.ndarray) -> np.ndarray:
        """The reduced coupling operator of a complete coupling field, for checks."""
        return -self.t * self.weighted(np.asarray(field, dtype=float)[:self.radius.size])


class ROM:
    """Fusion cross sections at any (theta, E) inside the trained box and energy range."""

    def __init__(self, reaction: Reaction, eim: MultipoleEIM, radius: np.ndarray,
                 reduced_mass: float, zz: int, arrays: dict):
        self.reaction = reaction
        self.eim = eim
        self.radius = np.asarray(radius, dtype=float)
        self.h = float(self.radius[1] - self.radius[0])
        self.reduced_mass = float(reduced_mass)
        self.kfac = 2.0 * self.reduced_mass / HBARC**2
        self.h2m = HBARC**2 / (2.0 * self.reduced_mass)
        self.t = self.h**2 * self.kfac / 12.0
        self.zz = int(zz)
        self.channels = int(eim.static_field.shape[1])
        # Basis rows at r_min, r_min + dr, r_max, and r_max - 2 dr, for the boundary conditions.
        self.x0, self.x1, self.xn, self.xnm2 = (
            np.asarray(arrays[name], dtype=np.complex128) for name in ("x0", "x1", "xn", "xnm2"))
        # Reduced operators, R x R; a_modes holds one per CC-EIM mode.
        self.p_d2, self.p_one, self.p_coul, self.p_cent, self.a_static, self.a_modes = (
            np.asarray(arrays[name], dtype=np.complex128)
            for name in ("p_d2", "p_one", "p_coul", "p_cent", "a_static", "a_modes"))
        # The CC-EIM nodes of each multipole, as positions among the nodes of all multipoles.
        node_rows = eim.node_rows
        self._node_radii = field_radii(self.radius)[node_rows]
        self._node_positions = [np.searchsorted(node_rows, item.nodes) for item in eim.radial]
        # The outer Coulomb matching uses the nominal coupling near r_max, so it depends on (E, J) only.
        self._outer_radii = self.radius[[-3, -1]]
        self._outer_potentials = np.array(
            [np.diag(field) for field in reaction.library.fields(reaction.nominal, self._outer_radii)])
        self._outer: dict[float, tuple] = {}

    @property
    def rank(self) -> int:
        return int(self.p_d2.shape[0])

    def coefficients(self, point) -> np.ndarray:
        """CC-EIM coefficients from the multipole form factors at the interpolation nodes."""
        factors = self.reaction.library.factors(point, self._node_radii, self.eim.decomposition_inverse)
        return self.eim.coefficients([factors[positions, component]
                                      for component, positions in enumerate(self._node_positions)])

    def interaction(self, point) -> tuple:
        """Reduced coupling operator, and the coupling matrices at r_min and r_min + dr/2."""
        c = self.coefficients(point)
        operator = self.a_static + np.tensordot(c, self.a_modes, axes=1)
        p0, ph = self.eim.reconstruct(c, rows=[0, self.radius.size])
        return operator, p0, ph

    def outer_maps(self, energy: float) -> tuple:
        """Partial waves of ``energy`` and their outer Coulomb maps."""
        key = round(float(energy), 10)
        if key not in self._outer:
            waves = self.reaction.partial_waves(energy)
            maps = [self.reaction.library.outer_map(float(energy), int(j), self._outer_potentials,
                                                    self._outer_radii)
                    for j in waves]
            self._outer[key] = (waves, np.array([m[0] for m in maps]).astype(np.complex128),
                                np.array([m[1] for m in maps]).astype(np.complex128))
        return self._outer[key]

    def _inner_maps(self, p0, ph, energy: float, waves: np.ndarray):
        """RK4 start of the incoming-wave boundary condition at r_min for every J, as CCFULL."""
        identity = np.eye(self.channels)
        radii = np.array([self.radius[0], self.radius[0] + self.h / 2.0])
        lam = waves * (waves + 1.0)
        scalar = self.zz * HBARC / 137.0 / radii[None, :] + self.h2m * lam[:, None] / radii[None, :]**2
        q0 = self.kfac * (p0[None] + (scalar[:, 0] - energy)[:, None, None] * identity)
        qh = self.kfac * (ph[None] + (scalar[:, 1] - energy)[:, None, None] * identity)
        channel_energy = energy - scalar[:, :1] - np.diag(p0)[None, :]
        root = np.sqrt(self.kfac * np.abs(channel_energy))
        derivative = np.where(channel_energy > 0.0, -1j * root, root)
        d0 = derivative[:, :, None] * identity
        h = self.h
        k1y, k1v = h * d0, h * q0
        k2y = h * (d0 + 0.5 * k1v)
        k2v = h * qh @ (identity + 0.5 * k1y)
        k3y = h * (d0 + 0.5 * k2v)
        k3v = h * qh @ (identity + 0.5 * k2y)
        k4y = h * (d0 + k3v)
        return identity + (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0, channel_energy

    def penetrabilities(self, interaction: tuple, energy: float) -> tuple[np.ndarray, np.ndarray]:
        """P_J of every partial wave of ``energy``: one R x R solve per J."""
        operator, p0, ph = interaction
        waves, ratio, source = self.outer_maps(energy)
        inner, channel_energy = self._inner_maps(p0, ph, energy, waves)
        lam = waves * (waves + 1.0)
        base = self.p_d2 - self.t * (self.p_coul - energy * self.p_one) + operator
        matrices = base[None] - (self.t * self.h2m) * lam[:, None, None] * self.p_cent[None]
        x0, x1, xn, xnm2 = self.x0, self.x1, self.xn, self.xnm2
        matrices += x0.conj().T @ (x1[None] - inner @ x0)
        matrices += xn.conj().T @ (xn[None] - ratio @ xnm2)
        rhs = source @ xn.conj()
        solution = np.linalg.solve(matrices, rhs[..., None])[..., 0]
        amplitudes = solution @ x0.T
        weight = np.sqrt(np.maximum(channel_energy, 0.0) / energy) * (channel_energy > 0.0)
        return waves, np.sum(np.abs(amplitudes)**2 * weight, axis=1)

    def cross_sections(self, point, energies) -> np.ndarray:
        """sigma_fus(E) in mb at one parameter point, for energies inside the trained range."""
        energies = np.asarray(energies, dtype=float).reshape(-1)
        lower, upper = self.reaction.energy_range
        if np.any(energies < lower - 1.0e-9) or np.any(energies > upper + 1.0e-9):
            raise ValueError(f"energies must lie in the trained range {lower:g}-{upper:g} MeV")
        interaction = self.interaction(point)
        values = []
        for energy in energies:
            waves, penetrabilities = self.penetrabilities(interaction, float(energy))
            values.append(fusion_cross_section(penetrabilities, waves, float(energy), self.reduced_mass))
        return np.array(values)

    def exact_cross_sections(self, point, energies) -> np.ndarray:
        """Exact Numerov CCFULL at the same point and energies, for accuracy checks."""
        return self.reaction.exact_cross_sections(point, energies, self.reduced_mass)

    def save(self, path: Path) -> None:
        np.savez(path, radius=self.radius, reduced_mass=self.reduced_mass, zz=self.zz,
                 **{name: getattr(self, name) for name in ARRAYS})

    @classmethod
    def load(cls, manifest: str | Path) -> "ROM":
        path = Path(manifest)
        payload = read_manifest(path)
        eim = load_eim(path.parent / payload["shared_eim"])
        with np.load(path.parent / payload["operators"]) as data:
            arrays = {name: data[name] for name in ARRAYS}
            radius, reduced_mass, zz = data["radius"], float(data["reduced_mass"]), int(data["zz"])
        return cls(Reaction.from_settings(payload["settings"]), eim, radius, reduced_mass, zz, arrays)
