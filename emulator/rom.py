"""The reduced model of CCFULL: what the online stage computes at a parameter point.

For every partial wave J at energy E, CCFULL solves the coupled equations on a radial grid
with the standard Numerov step.  The reduced model writes these discrete equations for all
grid points as one linear system A u = b and projects it onto a small basis W:

    (W^H A W) x = W^H b,   u = W x.

A depends on the parameter point theta only through the coupling field, which the CC-EIM
writes as
    V(r; theta) = V_static(r) + sum_m c_m(theta) B_m(r),
so W^H A W is a fixed matrix plus sum_m c_m(theta) times fixed matrices, all projected once
in the training.  Online, one parameter point costs
    1. the CC-EIM coefficients c_m: CCFULL's coupling matrix at the node radii (coupling.py);
    2. for every (E, J): the rows of the inner boundary (CCFULL's Runge-Kutta start at r_min),
       the sum over the modes, and one small linear solve;
    3. P_J from the solution at r_min, and sigma from CCFULL's J sum.

A local model has its own basis for every (E, J) of its training energies (rank r, 24-80).
A global model has one basis for all (E, J) (rank R, about 300) and works at any energy of
its range; it needs the outer boundary of each energy, from CCFULL, before the timed points.
"""

import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from emulator.ccfull_run import PROJECT, data_lines, project_path, wave_functions
from emulator.coupling import PI, coupling_matrix, partial_wave_count, reaction


HBARC = 197.329                               # MeV fm, as in CCFULL
PARAMETER_NAMES = ("beta2", "beta4", "beta6", "v0", "r0", "a0")


def latin_hypercube(box, count, seed):
    """count Latin-hypercube points in box (6 x 2); a parameter with equal bounds stays fixed."""
    box = np.asarray(box, dtype=float)
    points = np.repeat(box[:, 0][None, :], count, axis=0)
    free = box[:, 1] > box[:, 0]
    unit = qmc.LatinHypercube(int(free.sum()), seed=seed).random(count)
    points[:, free] = qmc.scale(unit, box[free, 0], box[free, 1])
    return points


def ccfull_energy(energy):
    """CCFULL keeps E in single precision (Real :: E), so the reduced model uses the same value."""
    return float(np.float32(energy))


# ----------------------------------------------------------------- reading a trained model

def load_model(folder, rank=0):
    """A trained model as a dict; rank > 0 keeps the leading rank vectors of each basis."""
    folder = Path(folder)
    manifest = json.loads((folder / "manifest.json").read_text())
    settings = manifest["settings"]
    model = {"name": folder.name, "folder": folder,
             "type": "local" if manifest["format"].startswith("local") else "global",
             "trained_rank": int(manifest["rank"]),
             "box": np.array(settings["parameter_bounds"], dtype=float),
             "energy_range": tuple(settings["energy_range_mev"]),
             "experimental_file": settings["experimental_file"]}
    model["rank"] = int(rank) or model["trained_rank"]
    if model["rank"] > model["trained_rank"]:
        raise ValueError(f"rank {model['rank']} is above the trained rank {model['trained_rank']} of {folder.name}")

    # The CCFULL records of the training.  A model trained by this code keeps a copy of them and
    # of the interface program; it is refused once either has changed.
    model["records"] = data_lines(project_path(settings["input_file"]))[:9]
    if (folder / "ccfull.inp").exists() and data_lines(folder / "ccfull.inp")[:9] != model["records"]:
        raise ValueError(f"the CCFULL records of {settings['input_file']} changed since {folder.name} "
                         "was trained; restore them or train the model again")
    source = PROJECT / "physical_model" / "ccfull_interface" / "ccfull.f90"
    if (folder / "ccfull_interface.f90").exists() and (folder / "ccfull_interface.f90").read_text() != source.read_text():
        raise ValueError(f"{source} changed since {folder.name} was trained; train the model again")
    model["system"] = reaction(model["records"])

    # The radial grid; its first point is the fixed r_min of every CCFULL run of the model.
    grid_file = folder / (manifest["grid"] if model["type"] == "local" else manifest["operators"])
    with np.load(grid_file) as data:
        model["radius"] = np.array(data["radius"])
        model["reduced_mass"] = float(data["reduced_mass"])
        model["zz"] = int(data["zz"])
    model["channels"] = model["system"]["channels"]

    # The CC-EIM: the coupling field on the grid and at r_min + dr/2 (the last field row).
    eim_folder = folder / manifest["shared_eim"]
    with np.load(eim_folder / "geometry.npz") as data:
        geometries, whitening, static_field = data["geometries"], data["whitening"], data["static_field"]
    with np.load(eim_folder / "joint.npz") as data:
        basis, nodes, inverse = data["basis"], data["nodes"].astype(int), data["inverse"]
    radius = model["radius"]
    field_radii = np.append(radius, radius[0] + 0.5 * (radius[1] - radius[0]))
    node_rows = np.unique(nodes[:, 0])
    model["eim"] = {"geometries": geometries, "whitening": whitening, "static_field": static_field,
                    "basis": basis, "nodes": nodes, "inverse": inverse, "node_rows": node_rows,
                    "node_radii": field_radii[node_rows]}

    r = model["rank"]
    if model["type"] == "local":
        # one file per energy, read before the timed points (load_blocks)
        model["files"] = {round(float(item["energy_mev"]), 6): folder / item["file"] for item in manifest["files"]}
        model["blocks"] = {}
    else:
        with np.load(grid_file) as data:
            model["operators"] = {name: np.ascontiguousarray(data[name][..., :r] if name in ("x0", "x1", "xn", "xnm2")
                                                             else data[name][..., :r, :r])
                                  for name in ("x0", "x1", "xn", "xnm2", "p_d2", "p_one", "p_coul",
                                               "p_cent", "a_static", "a_modes")}
        model["outer"] = {}
    return model


def load_blocks(model, energies):
    """Read the blocks of a local model at energies, keeping the leading rank vectors."""
    r = model["rank"]
    for energy in energies:
        key = round(float(energy), 6)
        if key not in model["files"]:
            raise ValueError(f"{energy} MeV is not an energy of the local model {model['name']}")
        with np.load(model["files"][key]) as data:
            model["blocks"][key] = {
                "J": np.array(data["waves"]),
                "a0": np.ascontiguousarray(data["a0"][:, :r, :r]),
                "a_modes": np.ascontiguousarray(data["a_modes"][:, :, :r, :r]),   # complex64
                "x0": np.ascontiguousarray(data["x0"][..., :r]),
                "x1": np.ascontiguousarray(data["x1"][..., :r]),
                "rhs": np.ascontiguousarray(data["rhs"][:, :r])}


def nominal_point(model):
    """The six parameters of the CCFULL records of the model (lines 3 and 7)."""
    target = [float(value) for value in model["records"][2].replace(",", " ").split()]
    potential = [float(value) for value in model["records"][6].replace(",", " ").split()]
    return np.array([target[1], target[2], target[3], *potential[:3]])


# ----------------------------------------------------------------- the outer boundary

def coulomb_all_waves(f0, g0, f1, g1, eta, rho, count):
    """F_L and G_L for L = 0 .. count-1 from L = 0 and 1, by the recurrence (Abramowitz and
    Stegun 14.2.3).  It is stable where r_max lies outside the turning point, as here."""
    f, g = [f0, f1], [g0, g1]
    for L in range(1, count - 1):
        a = L * math.sqrt((L + 1) ** 2 + eta ** 2)
        b = (2 * L + 1) * (eta + L * (L + 1) / rho)
        c = (L + 1) * math.sqrt(L ** 2 + eta ** 2)
        f.append((b * f[L] - c * f[L - 1]) / a)
        g.append((b * g[L] - c * g[L - 1]) / a)
    return np.array(f[:count]), np.array(g[:count])


def outer_boundary(model, matching, count):
    """CCFULL's matching to Coulomb waves at r1 = r_max - dr and r2 = r_max + dr, for
    J = 0 .. count-1: u_c(r2) - ratio_c u_c(r1) = source_c, with an incoming wave in the
    entrance channel only.  matching: the rows for J = 0 and 1 of the interface program at
    the nominal point (the channel energies there do not depend on J)."""
    radius, mu = model["radius"], model["reduced_mass"]
    r1, r2 = radius[-3], radius[-1]
    zz = model["system"]["zpro"] * model["system"]["ztar"]
    ratio = np.zeros((count, model["channels"]), dtype=complex)
    source = np.zeros((count, model["channels"]), dtype=complex)
    for channel in range(model["channels"]):
        row0 = matching[(matching[:, 1] == 0) & (matching[:, 2] == channel + 1)][0]
        row1 = matching[(matching[:, 1] == 1) & (matching[:, 2] == channel + 1)][0]
        ec1, ec2 = row0[3], row0[4]
        if ec1 < 0.0:                             # closed: the decaying exponential
            ratio[:, channel] = math.exp(-math.sqrt(2.0 * mu * abs(ec1)) / HBARC * (r2 - r1))
            continue
        waves = []
        for ec, r, k in ((ec1, r1, 5), (ec2, r2, 7)):     # columns of F and G at r1, then r2
            eta = zz / 137.0 * math.sqrt(mu / (2.0 * ec))
            rho = math.sqrt(2.0 * mu * ec) / HBARC * r
            waves.append(coulomb_all_waves(row0[k], row0[k + 1], row1[k], row1[k + 1], eta, rho, count))
        (f1, g1), (f2, g2) = waves
        ratio[:, channel] = (g2 + 1j * f2) / (g1 + 1j * f1)
        if channel == 0:
            source[:, 0] = (g2 - 1j * f2) - ratio[:, 0] * (g1 - 1j * f1)
    return ratio, source


def prepare_global(model, energies):
    """For a global model: the partial waves of every energy (J = 0 .. n-1, n the largest count
    of CCFULL's J loop over the corners of the (V0, r0, a0) box) and their outer boundary
    (from one run of the interface program at the nominal point)."""
    nominal = nominal_point(model)
    corners = list(itertools.product(*model["box"][3:]))
    for energy in energies:
        e = ccfull_energy(energy)
        count = max(partial_wave_count(e, v0, r0, a0, model["system"]) for v0, r0, a0 in corners)
        run = wave_functions(model["records"], nominal, energy, model["radius"][0])
        ratio, source = outer_boundary(model, run["matching"], count)
        model["outer"][round(float(energy), 6)] = (np.arange(count), ratio, source)


# ----------------------------------------------------------------- one parameter point

def eim_coefficients(model, theta):
    """c_m(theta): CCFULL's coupling matrix at the node radii, its elements i <= j, and the
    interpolation at the nodes."""
    eim = model["eim"]
    field = coupling_matrix(eim["node_radii"], theta, model["system"])
    upper = np.triu_indices(model["channels"])
    elements = field[:, upper[0], upper[1]] @ eim["whitening"]
    rows = np.searchsorted(eim["node_rows"], eim["nodes"][:, 0])
    return eim["inverse"] @ elements[rows, eim["nodes"][:, 1]]


def coupling_field(model, coefficients, rows):
    """The interpolated coupling field V_static + sum_m c_m B_m on the field rows rows."""
    eim = model["eim"]
    elements = (eim["basis"][rows] @ coefficients) @ np.linalg.inv(eim["whitening"])
    return eim["static_field"][rows] + np.einsum("rl,lij->rij", elements, eim["geometries"])


def inner_boundary(model, p0, ph, energy, J):
    """CCFULL's start at r_min (rkutta00): one Runge-Kutta step from r_min to r_min + dr of an
    incoming wave in every open channel, for every J; and sqrt(E_c/E) of the open channels,
    which turns the amplitudes at r_min into P_J.  p0, ph: the coupling field at r_min and
    r_min + dr/2."""
    channels, radius, mu = model["channels"], model["radius"], model["reduced_mass"]
    h = radius[1] - radius[0]
    kfac = 2.0 * mu / HBARC ** 2
    identity = np.eye(channels)
    radii = np.array([radius[0], radius[0] + h / 2.0])
    lam = J * (J + 1.0)
    scalar = model["zz"] * HBARC / 137.0 / radii[None, :] + HBARC ** 2 / (2.0 * mu) * lam[:, None] / radii[None, :] ** 2
    q0 = kfac * (p0[None] + (scalar[:, 0] - energy)[:, None, None] * identity)
    qh = kfac * (ph[None] + (scalar[:, 1] - energy)[:, None, None] * identity)
    channel_energy = energy - scalar[:, :1] - np.diag(p0)[None, :]
    root = np.sqrt(kfac * np.abs(channel_energy))
    d0 = np.where(channel_energy > 0.0, -1j * root, root)[:, :, None] * identity
    k1y, k1v = h * d0, h * q0
    k2y = h * (d0 + 0.5 * k1v)
    k2v = h * qh @ (identity + 0.5 * k1y)
    k3y = h * (d0 + 0.5 * k2v)
    k3v = h * qh @ (identity + 0.5 * k2y)
    k4y = h * (d0 + k3v)
    weight = np.sqrt(np.maximum(channel_energy, 0.0) / energy) * (channel_energy > 0.0)
    return identity + (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0, weight


def local_penetrabilities(model, coefficients, p0, ph, energy, block):
    """P_J of every partial wave of a local block: one rank x rank solve per J."""
    inner, weight = inner_boundary(model, p0, ph, energy, block["J"])
    modes = block["a_modes"]                  # summed in its own (single) precision
    matrices = block["a0"] + np.tensordot(coefficients.astype(modes.real.dtype), modes, axes=1)
    x0, x1 = block["x0"], block["x1"]
    matrices = matrices + x0.conj().transpose(0, 2, 1) @ (x1 - inner @ x0)
    solution = np.linalg.solve(matrices, block["rhs"][..., None])[..., 0]
    amplitudes = np.einsum("jr,jcr->jc", solution, x0)
    return np.sum(np.abs(amplitudes) ** 2 * weight, axis=1)


def global_penetrabilities(model, interaction, p0, ph, energy, outer):
    """P_J of every partial wave of energy with the global basis: one R x R solve per J."""
    J, ratio, source = outer
    o = model["operators"]
    radius, mu = model["radius"], model["reduced_mass"]
    t = (radius[1] - radius[0]) ** 2 * (2.0 * mu / HBARC ** 2) / 12.0
    inner, weight = inner_boundary(model, p0, ph, energy, J)
    base = o["p_d2"] - t * (o["p_coul"] - energy * o["p_one"]) + interaction
    matrices = base[None] - (t * HBARC ** 2 / (2.0 * mu)) * (J * (J + 1.0))[:, None, None] * o["p_cent"][None]
    matrices = matrices + o["x0"].conj().T @ (o["x1"][None] - inner @ o["x0"])
    matrices = matrices + o["xn"].conj().T @ (o["xn"][None] - ratio[:, :, None] * o["xnm2"][None])
    rhs = source @ o["xn"].conj()
    solution = np.linalg.solve(matrices, rhs[..., None])[..., 0]
    amplitudes = solution @ o["x0"].T
    return np.sum(np.abs(amplitudes) ** 2 * weight, axis=1)


def partial_waves_in_sum(P, J):
    """How many partial waves CCFULL's J sum takes: it stops at the first term below 1e-4 of
    the sum so far (that term included)."""
    total = 0.0
    for count, (angular_momentum, value) in enumerate(zip(J, P), start=1):
        term = (2.0 * angular_momentum + 1.0) * value
        if total > 0.0 and term < total * 1.0e-4:
            return count
        total += term
    return len(J)


def fusion_cross_section(P, J, energy, reduced_mass):
    """sigma (mb) = sum_J (2J+1) P_J pi hbar^2 / (2 mu E) x 10, with CCFULL's rule for the end
    of the sum: it stops once a term is below 1e-4 of the sum so far."""
    sigma = 0.0
    for angular_momentum, value in zip(J, P):
        term = (2.0 * angular_momentum + 1.0) * value * PI * HBARC ** 2 / (2.0 * reduced_mass * energy) * 10.0
        previous = sigma
        sigma += term
        if term < previous * 1.0e-4:
            break
    return sigma


def cross_sections(model, theta, energies):
    """sigma (mb) of the reduced model at one parameter point and the given energies."""
    coefficients = eim_coefficients(model, theta)
    p0, ph = coupling_field(model, coefficients, [0, model["radius"].size])
    if model["type"] == "global":
        o = model["operators"]
        interaction = o["a_static"] + np.tensordot(coefficients, o["a_modes"], axes=1)
    sigma = []
    for energy in energies:
        key, e = round(float(energy), 6), ccfull_energy(energy)
        if model["type"] == "local":
            block = model["blocks"][key]
            P, J = local_penetrabilities(model, coefficients, p0, ph, e, block), block["J"]
        else:
            P = global_penetrabilities(model, interaction, p0, ph, e, model["outer"][key])
            J = model["outer"][key][0]
        sigma.append(fusion_cross_section(P, J, e, model["reduced_mass"]))
    return np.array(sigma)


# ----------------------------------------------------------------- training (train_local.py, train_global.py)

def field_radii(radius):
    """The rows of the coupling field: the radial grid, then r_min + dr/2 for CCFULL's start."""
    return np.append(radius, radius[0] + 0.5 * (radius[1] - radius[0]))


def maxvol(basis, iterations=100, tolerance=1.0e-10):
    """The interpolation rows of basis: pivoted QR, improved by MaxVol row swaps."""
    from scipy.linalg import qr
    rank = basis.shape[1]
    rows = np.asarray(qr(basis.T, pivoting=True, mode="economic")[2][:rank], dtype=int).copy()
    for _ in range(iterations):
        coefficients = basis @ np.linalg.inv(basis[rows, :])
        row, column = np.unravel_index(np.argmax(np.abs(coefficients)), coefficients.shape)
        if abs(coefficients[row, column]) <= 1.0 + tolerance:
            break
        rows[column] = row
    return rows


def fit_eim(system, radius, box, count, seed, modes):
    """The joint CC-EIM of the coupling field, from coupling.py at count Latin-hypercube points.

    All independent matrix elements f_ij (i <= j) of the field on the field rows are interpolated
    together: one POD of the stacked elements, in the norm of the coupling matrix, and one set of
    MaxVol nodes, each a pair (field row, element).  Returns the CC-EIM (as load_model keeps it)
    and diagnostics: its error at 16 new points, and how many modes one interpolation per element
    would need for the same relative singular value.
    """
    channels = system["channels"]
    upper = np.triu_indices(channels)
    geometries = np.zeros((upper[0].size, channels, channels))
    for k, (i, j) in enumerate(zip(*upper)):
        geometries[k, i, j] = geometries[k, j, i] = 1.0
    whitening = np.linalg.cholesky(np.einsum("lij,kij->lk", geometries, geometries))
    radii = field_radii(radius)
    rotor = np.diag([(2 * k) * (2 * k + 1) * system["e2"] / 6.0 for k in range(channels)])
    static_field = np.broadcast_to(rotor, (radii.size, channels, channels)).copy()
    elements = lambda theta: coupling_matrix(radii, theta, system)[:, upper[0], upper[1]]

    samples = np.array([elements(theta) for theta in latin_hypercube(box, count, seed)])
    whitened = samples @ whitening
    left, singular, _ = np.linalg.svd(whitened.reshape(count, -1).T, full_matrices=False)
    basis = left[:, :min(modes, singular.size)]
    flat = maxvol(basis)
    nodes = np.column_stack(np.unravel_index(flat, whitened.shape[1:]))
    node_rows = np.unique(nodes[:, 0])
    eim = {"geometries": geometries, "whitening": whitening, "static_field": static_field,
           "basis": basis.reshape(*whitened.shape[1:], basis.shape[1]), "singular_values": singular,
           "nodes": nodes, "inverse": np.linalg.inv(basis[flat]), "node_rows": node_rows,
           "node_radii": radii[node_rows]}

    errors = []
    for theta in latin_hypercube(box, 16, seed + 1):
        exact = elements(theta) @ whitening
        rows = np.searchsorted(node_rows, nodes[:, 0])
        c = eim["inverse"] @ exact[node_rows][rows, nodes[:, 1]]
        errors.append(np.linalg.norm(eim["basis"] @ c - exact) / np.linalg.norm(exact))
    tolerance = singular[basis.shape[1] - 1] / singular[0]
    separate = {f"{i + 1}-{j + 1}": int(np.count_nonzero(np.linalg.svd(samples[:, :, k], compute_uv=False)
                                                         >= tolerance * np.linalg.svd(samples[:, :, k], compute_uv=False)[0]))
                for k, (i, j) in enumerate(zip(*upper))}
    diagnostics = {"interaction_training_points": count, "interaction_seed": seed, "components": "matrix elements",
                   "modes": int(basis.shape[1]), "relative_singular_value_of_last_mode": float(tolerance),
                   "node_rows": int(node_rows.size), "modes_one_interpolation_per_component": sum(separate.values()),
                   "modes_per_component_if_separate": separate,
                   "interpolation_relative_error": {"median": float(np.median(errors)), "max": float(np.max(errors))}}
    return eim, diagnostics


def save_eim(eim, diagnostics, folder):
    folder.mkdir(parents=True, exist_ok=True)
    np.savez(folder / "geometry.npz", geometries=eim["geometries"], whitening=eim["whitening"],
             static_field=eim["static_field"])
    np.savez(folder / "joint.npz", basis=eim["basis"], singular_values=eim["singular_values"],
             nodes=eim["nodes"].astype(np.int32), inverse=eim["inverse"])
    (folder / "metadata.json").write_text(json.dumps(diagnostics, indent=2) + "\n")


def mode_fields(eim):
    """The channel fields B_m of the CC-EIM modes, (mode, field row, channel, channel)."""
    elements = np.einsum("rkm,kl->mrl", eim["basis"], np.linalg.inv(eim["whitening"]))
    return np.einsum("mrl,lij->mrij", elements, eim["geometries"])


def project_numerov(basis, radius, reduced_mass, zz, channels, static_field, fields):
    """CCFULL's Numerov equations projected onto a basis W (columns over radius x channel).

    For u on the grid the standard Numerov equations are
        u(n+1) - 2u(n) + u(n-1) - [T(n+1)u(n+1) + 10 T(n)u(n) + T(n-1)u(n-1)] = 0,
    T = dr^2 (2 mu/hbar^2)(V - E)/12, for the interior points.  They are linear in V and E, so
    each piece of V (point Coulomb, centrifugal, the static field, each CC-EIM mode) and E gives
    one projected matrix W^H (...) W.  The rows of W at r_min, r_min + dr, r_max + dr, and
    r_max - dr carry the boundary conditions.
    """
    x = np.asarray(basis, dtype=np.complex128).reshape(radius.size, channels, -1)
    t = (radius[1] - radius[0]) ** 2 * (2.0 * reduced_mass / HBARC ** 2) / 12.0
    project = lambda residual: np.tensordot(x.conj(), residual, axes=([0, 1], [0, 1]))

    def weighted(f):                       # W^H of f(n+1)u(n+1) + 10 f(n)u(n) + f(n-1)u(n-1)
        residual = np.zeros_like(x)
        if np.ndim(f) == 1:
            residual[1:-1] = f[2:, None, None] * x[2:] + 10.0 * f[1:-1, None, None] * x[1:-1] + f[:-2, None, None] * x[:-2]
        else:
            residual[1:-1] = f[2:] @ x[2:] + 10.0 * (f[1:-1] @ x[1:-1]) + f[:-2] @ x[:-2]
        return project(residual)

    second = np.zeros_like(x)
    second[1:-1] = x[2:] - 2.0 * x[1:-1] + x[:-2]
    n = radius.size
    return {"x0": x[0], "x1": x[1], "xn": x[-1], "xnm2": x[-3],
            "p_d2": project(second), "p_one": weighted(np.ones(n)),
            "p_coul": weighted(zz * HBARC / 137.0 / radius), "p_cent": weighted(1.0 / radius ** 2),
            "a_static": -t * weighted(np.asarray(static_field[:n], dtype=float)),
            "a_modes": np.stack([-t * weighted(field[:n]) for field in fields])}


def local_block(snapshots, J, energy, rank, grid, fields, ratio, source):
    """The block of one (E, J) of a local model.

    snapshots: (count, radius * channel) unit-norm wave functions of this (E, J).  Their leading
    POD vectors are the basis W; with fewer snapshots than rank, the basis is smaller and the
    block is filled up with an identity that the solve leaves at zero.  Everything but the
    CC-EIM sum and the inner boundary rows is fixed at this (E, J) and kept in a0.
    """
    radius, mu, channels = grid["radius"], grid["reduced_mass"], grid["channels"]
    left, singular, _ = np.linalg.svd(np.asarray(snapshots, dtype=np.complex128).T, full_matrices=False)
    size = min(rank, left.shape[1])
    ops = project_numerov(left[:, :size], radius, mu, grid["zz"], channels, grid["static_field"], fields)
    t = (radius[1] - radius[0]) ** 2 * (2.0 * mu / HBARC ** 2) / 12.0
    xn, xnm2 = ops["xn"], ops["xnm2"]
    a0 = (ops["p_d2"] - t * (ops["p_coul"] - energy * ops["p_one"]) + ops["a_static"]
          - t * HBARC ** 2 / (2.0 * mu) * J * (J + 1.0) * ops["p_cent"]
          + xn.conj().T @ (xn - ratio[:, None] * xnm2))
    block = {"a0": np.eye(rank, dtype=complex), "a_modes": np.zeros((len(fields), rank, rank), dtype=complex),
             "x0": np.zeros((channels, rank), dtype=complex), "x1": np.zeros((channels, rank), dtype=complex),
             "rhs": np.zeros(rank, dtype=complex)}
    block["a0"][:size, :size] = a0
    block["a_modes"][:, :size, :size] = ops["a_modes"]
    block["x0"][:, :size], block["x1"][:, :size] = ops["x0"], ops["x1"]
    block["rhs"][:size] = source @ xn.conj()
    block["tail"] = float(np.sqrt(max(0.0, 1.0 - np.sum(singular[:size] ** 2) / np.sum(singular ** 2))))
    block["size"] = size
    return block
