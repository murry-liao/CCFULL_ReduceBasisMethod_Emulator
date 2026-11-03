"""Multipole CC-EIM of the CCFULL coupling field.

The coupling matrix of the axial rotor target is written as

    V(r; theta) = V_static(r) + sum_lambda G_lambda f_lambda(r; theta),

with a fixed channel matrix G_lambda for every even multipole lambda and a
scalar radial form factor f_lambda (Catacora-Rios et al., Eqs. (27)-(30) and
(B15)-(B28)).  Each f_lambda has its own empirical interpolation: POD modes in
r and MaxVol nodes, so a parameter point needs f_lambda at the nodes only.

The field has one row per radial grid point and a last row at r_min + dr/2,
where the RK4 start of CCFULL needs it.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import qr


FORMAT = "ccfull_multipole_eim_v1"


def field_radii(radius: np.ndarray) -> np.ndarray:
    """Radii of the coupling-field rows: the radial grid, then r_min + dr/2."""
    return np.append(radius, radius[0] + 0.5 * (radius[1] - radius[0]))


def _maxvol(basis: np.ndarray, max_iterations: int = 100, tolerance: float = 1.0e-10) -> np.ndarray:
    """Interpolation rows of ``basis``: pivoted QR, improved by MaxVol row swaps."""
    rank = basis.shape[1]
    _, _, pivots = qr(basis.T, pivoting=True, mode="economic")
    rows = np.asarray(pivots[:rank], dtype=int).copy()
    for _ in range(max_iterations):
        coefficients = basis @ np.linalg.inv(basis[rows, :])
        row, column = np.unravel_index(np.argmax(np.abs(coefficients)), coefficients.shape)
        if abs(coefficients[row, column]) <= 1.0 + tolerance:
            break
        rows[column] = row
    return rows


@dataclass
class RadialEIM:
    """Empirical interpolation of one radial form factor."""

    basis: np.ndarray              # POD modes, (radius, rank)
    singular_values: np.ndarray
    nodes: np.ndarray              # interpolation rows of basis
    inverse: np.ndarray            # inverse of basis[nodes]

    @classmethod
    def fit(cls, samples: np.ndarray, rank: int) -> "RadialEIM":
        """``rank`` POD modes of ``samples`` (sample, radius) and their MaxVol nodes."""
        left, singular_values, _ = np.linalg.svd(samples.T, full_matrices=False)
        basis = left[:, :rank]
        nodes = _maxvol(basis)
        return cls(basis, singular_values, nodes, np.linalg.inv(basis[nodes, :]))


@dataclass
class MultipoleEIM:
    """CC-EIM of the coupling field: fixed channel matrices and a radial EIM per multipole."""

    multipoles: tuple[int, ...]
    geometries: np.ndarray             # G_lambda, (multipole, channel, channel)
    radial: tuple[RadialEIM, ...]
    static_field: np.ndarray           # V_static on every field row, (radius + 1, channel, channel)
    decomposition_inverse: np.ndarray  # coupling matrix, flattened -> f_lambda: (multipole, channel**2)
    decomposition_median_relative_error: float = 0.0
    decomposition_max_relative_error: float = 0.0

    @property
    def ranks(self) -> list[int]:
        return [item.basis.shape[1] for item in self.radial]

    @property
    def rank(self) -> int:
        return sum(self.ranks)

    @property
    def node_rows(self) -> np.ndarray:
        """Field rows where at least one multipole has an interpolation node."""
        return np.unique(np.concatenate([item.nodes for item in self.radial])).astype(int)

    def coefficients(self, samples: list[np.ndarray]) -> np.ndarray:
        """EIM coefficients from each f_lambda at the nodes of its multipole, in multipole order."""
        return np.concatenate([item.inverse @ values for item, values in zip(self.radial, samples)])

    def reconstruct(self, coefficients: np.ndarray, rows=None) -> np.ndarray:
        """The interpolated coupling field on every field row, or on ``rows`` only."""
        static = self.static_field if rows is None else self.static_field[rows]
        field = static.astype(np.result_type(static.dtype, coefficients.dtype), copy=True)
        parts = np.split(coefficients, np.cumsum(self.ranks)[:-1])
        for geometry, item, part in zip(self.geometries, self.radial, parts):
            basis = item.basis if rows is None else item.basis[rows]
            field += (basis @ part)[:, None, None] * geometry[None, :, :]
        return field


def _wigner_3j_zero(j1: int, j2: int, j3: int) -> float:
    """The Wigner 3j symbol (j1 j2 j3; 0 0 0) of integer spins."""
    total = int(j1 + j2 + j3)
    if total % 2 or j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0
    half = total // 2
    if min(half - j1, half - j2, half - j3) < 0:
        return 0.0
    prefactor = ((-1.0) ** half * math.factorial(half)
                 / (math.factorial(half-j1) * math.factorial(half-j2)
                    * math.factorial(half-j3)))
    root = math.sqrt(
        math.factorial(total-2*j1) * math.factorial(total-2*j2)
        * math.factorial(total-2*j3) / math.factorial(total+1))
    return prefactor * root


def rotor_geometry(config: dict, channels: int) -> tuple[tuple[int, ...], np.ndarray]:
    """The multipoles 0, 2, ..., 4 Nrot and their channel matrices G_lambda.

    The CCFULL rotor band has I = 0, 2, ..., 2 Nrot.  Every even multipole up to
    2 Imax is kept: that is all the operator content of this channel space,
    including the higher multipoles of the all-order deformed Woods-Saxon potential.
    """
    if int(config["target_mode"]) != 1 or int(config["projectile_mode"]) != -1:
        raise ValueError("the multipole CC-EIM supports an axial rotor target and an inert projectile")
    nrot = int(config["target1"][4])
    if channels != nrot + 1:
        raise ValueError("the channel count does not match the target rotational band")
    multipoles = tuple(range(0, 4*nrot + 1, 2))
    geometry = np.zeros((len(multipoles), channels, channels), dtype=float)
    for component, multipole in enumerate(multipoles):
        for row in range(channels):
            initial_spin = 2 * row
            for column in range(channels):
                final_spin = 2 * column
                symbol = _wigner_3j_zero(initial_spin, multipole, final_spin)
                geometry[component, row, column] = math.sqrt(
                    (2*initial_spin+1) * (2*multipole+1)
                    * (2*final_spin+1) / (4.0*math.pi)) * symbol**2
    return multipoles, geometry


def rotor_static_field(config: dict, rows: int, channels: int) -> np.ndarray:
    """The rotational energies I(I+1) E2/6 on the diagonal, the same on every field row."""
    e2 = float(config["target1"][0])
    energies = np.asarray([(2*index) * (2*index+1) * e2 / 6.0 for index in range(channels)], dtype=float)
    return np.broadcast_to(np.diag(energies), (int(rows), channels, channels)).copy()


def save_eim(eim: MultipoleEIM, directory: Path) -> None:
    """Write geometry.npz, one lambda_XX.npz per multipole, and metadata.json."""
    directory.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(directory / "geometry.npz",
                        geometries=np.asarray(eim.geometries, dtype=np.float64),
                        static_field=np.asarray(eim.static_field, dtype=np.float64),
                        decomposition_inverse=np.asarray(eim.decomposition_inverse, dtype=np.float64))
    components = []
    for label, item in zip(eim.multipoles, eim.radial):
        name = f"lambda_{label:02d}.npz"
        np.savez_compressed(directory / name,
                            basis=np.asarray(item.basis, dtype=np.float64),
                            singular_values=np.asarray(item.singular_values, dtype=np.float64),
                            interpolation_indices=np.asarray(item.nodes, dtype=np.int32),
                            interpolation_inverse=np.asarray(item.inverse, dtype=np.float64))
        components.append({"multipole": int(label), "artifact": name, "rank": item.basis.shape[1]})
    metadata = {
        "format": FORMAT,
        "multipoles": list(eim.multipoles),
        "total_rank": eim.rank,
        "components": components,
        "decomposition_median_relative_error": eim.decomposition_median_relative_error,
        "decomposition_max_relative_error": eim.decomposition_max_relative_error,
    }
    (directory / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def load_eim(directory: Path) -> MultipoleEIM:
    metadata = json.loads((directory / "metadata.json").read_text())
    if metadata.get("format") != FORMAT:
        raise ValueError(f"unsupported CC-EIM format in {directory}")
    with np.load(directory / "geometry.npz") as data:
        geometries, static_field = data["geometries"], data["static_field"]
        decomposition_inverse = data["decomposition_inverse"]
    radial = []
    for component in metadata["components"]:
        with np.load(directory / component["artifact"]) as data:
            radial.append(RadialEIM(data["basis"], data["singular_values"],
                                    data["interpolation_indices"], data["interpolation_inverse"]))
    return MultipoleEIM(
        tuple(int(value) for value in metadata["multipoles"]), geometries, tuple(radial),
        static_field, decomposition_inverse,
        float(metadata.get("decomposition_median_relative_error", 0.0)),
        float(metadata.get("decomposition_max_relative_error", 0.0)))
