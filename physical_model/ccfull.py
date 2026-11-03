"""The standard Numerov CCFULL, compiled for CCFULL-ROM.

``solvers/standard_numerov/ccfull.f90`` is used unchanged, without its main
program, in two builds:

* ``snapshot_driver.f90`` makes a program that solves one (theta, E, J) and
  writes the wave function.  Training runs it once per snapshot.
* ``library_interface.f90`` makes a shared library that stays loaded in the
  Python process: coupling potentials, multipole form factors, Coulomb
  functions, partial-wave counts, and exact penetrabilities.

A build is reused while the source it was compiled from is unchanged.
"""

from __future__ import annotations

import ctypes
import fcntl
from functools import cache
import hashlib
import os
from pathlib import Path
import platform
import subprocess
import tempfile

import numpy as np

from .inputs import input_text


HBARC = 197.329          # MeV fm, as in CCFULL
NUCLEON_MASS = 938.0     # MeV, as in CCFULL

SOLVER = Path(__file__).resolve().parent / "solvers" / "standard_numerov"
SOURCE = SOLVER / "ccfull.f90"
BUILD = Path(__file__).resolve().parent / "build" / "rom_backend"
_ARCHITECTURES = {"aarch64": "arm64", "arm64": "arm64", "amd64": "x86_64", "x86_64": "x86_64"}


def source_sha256() -> str:
    """SHA-256 of ccfull.f90, the solver behind every snapshot and exact reference."""
    return hashlib.sha256(SOURCE.read_bytes()).hexdigest()


def machine() -> str:
    """CPU architecture of this Python."""
    name = platform.machine().lower()
    return _ARCHITECTURES.get(name, name)


@cache
def gfortran() -> str:
    """The first gfortran on PATH that generates code for this Python's architecture.

    The library is loaded into this Python process, so it must match it.  On
    Apple silicon an Intel gfortran under /usr/local can come first on PATH;
    it is skipped.
    """
    found = []
    for folder in os.environ.get("PATH", "").split(os.pathsep):
        compiler = Path(folder or ".") / "gfortran"
        if not (compiler.is_file() and os.access(compiler, os.X_OK)):
            continue
        try:
            triple = subprocess.run([str(compiler), "-dumpmachine"], capture_output=True,
                                    text=True, timeout=30, check=True).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            triple = "unknown"
        target = triple.split("-", 1)[0].lower()
        if _ARCHITECTURES.get(target, target) == machine():
            return str(compiler)
        found.append(f"{compiler} -> {target}")
    raise RuntimeError(
        f"no gfortran on PATH generates {machine()} code for this Python "
        f"(found: {', '.join(found) or 'none'}); install a native gfortran or put it first on PATH")


def _ccfull_without_main() -> str:
    source = SOURCE.read_text()
    if "Program CCFull_Main" not in source:
        raise RuntimeError(f"no 'Program CCFull_Main' in {SOURCE}")
    return source.split("Program CCFull_Main", 1)[0]


def _compile(source: str, generated: Path, output: Path, flags: list[str]) -> Path:
    """Compile ``source`` into ``output``, unless ``output`` was built from the same source."""
    generated.parent.mkdir(parents=True, exist_ok=True)
    with (generated.parent / ".build.lock").open("w") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)      # processes that start together build once
        if output.exists() and generated.is_file() and generated.read_text() == source:
            return output
        output.unlink(missing_ok=True)
        generated.write_text(source)
        completed = subprocess.run(
            [gfortran(), *flags, "-o", str(output), str(generated), "-llapack", "-lblas"],
            cwd=generated.parent, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if completed.returncode:
            raise RuntimeError(f"failed to compile {generated}:\n{completed.stdout}")
    return output


@cache
def snapshot_program() -> Path:
    """The one-snapshot program, compiled on first use."""
    source = _ccfull_without_main() + (SOLVER / "snapshot_driver.f90").read_text()
    return _compile(source, BUILD / "ccfull_snapshot_driver.f90",
                    BUILD / f"ccfull_snapshot_driver_{machine()}", ["-O2", "-ffree-line-length-none"])


def snapshot(config: dict, point, energy: float, angular_momentum: int) -> dict:
    """One exact CCFULL solve at (theta, E, J), in its own process.

    Returns the wave function (radius, channel) as complex64, the radial grid,
    the reduced mass, and Zp * Zt.
    """
    run_config = {**config, "target1": list(config["target1"]),
                  "emin": float(energy), "emax": float(energy), "de": 1.0}
    run_config["target1"][1:4] = [float(value) for value in point[:3]]
    with tempfile.TemporaryDirectory(prefix="ccfull-snapshot-") as temporary:
        work = Path(temporary)
        (work / "ccfull.inp").write_text(input_text(run_config))
        (work / "snapshot.dat").write_text(
            " ".join(f"{value:.17g}" for value in [*point, energy]) + f" {int(angular_momentum)}\n")
        # CCFULL asks interactive questions while it starts; every answer is "n".
        completed = subprocess.run([str(snapshot_program())], cwd=work, input="n\n" * 128, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=1800)
        if completed.returncode:
            raise RuntimeError(f"CCFULL snapshot failed:\n{completed.stdout[-8000:]}")
        meta = (work / "snapshot_meta.txt").read_text().split()
        wave = np.fromfile(work / "snapshot_wave.bin", dtype=np.complex64)
    nlevel, nradial = int(meta[0]), int(meta[1])
    rmin, dr, reduced_mass = map(float, meta[2:5])
    if wave.size != nradial * nlevel:
        raise RuntimeError(f"CCFULL wrote {wave.size} wave-function values, expected {nradial * nlevel}")
    return {"wave": wave.reshape(nradial, nlevel),        # radius-major, channel-minor
            "radius": rmin + dr * np.arange(nradial),
            "reduced_mass": reduced_mass, "zz": int(meta[5]) * int(meta[6])}


_D, _I = ctypes.c_double, ctypes.c_int
_INT_OUT, _DOUBLE_OUT = ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double)
_F64 = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
_I32 = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS")
_THETA = [_D] * 6        # beta2, beta4, beta6, V0, r0, a0
_SIGNATURES = {
    "ccbrom_initialize": [_INT_OUT, _INT_OUT, _INT_OUT, _DOUBLE_OUT, _DOUBLE_OUT],
    "ccbrom_evaluate_fields": [*_THETA, _I, _F64, _F64, _INT_OUT],
    "ccbrom_evaluate_factors": [*_THETA, _I, _F64, _I, _F64, _I, _F64, _INT_OUT],
    "ccbrom_coulomb_values": [_I, _I, *[_F64] * 6, _INT_OUT],
    "ccbrom_exact_penetrabilities": [*_THETA, _D, _I, _I32, _F64, _INT_OUT],
    "ccbrom_active_partial_waves": [_D, _D, _D, _D, _I, _INT_OUT],
}


def _check(status: ctypes.c_int, what: str) -> None:
    if status.value:
        raise RuntimeError(f"CCFULL {what} failed with status {status.value}")


class Library:
    """CCFULL of one input, loaded into this process.

    CCFULL keeps its parameters in Fortran module variables, so every call sets
    all the parameters it uses, and a process loads each input once (``library``).
    """

    def __init__(self, config: dict):
        if int(config["target_mode"]) != 1 or int(config["projectile_mode"]) != -1:
            raise ValueError("the CCFULL library supports an axial rotor target and an inert projectile")
        self.config = config
        text = input_text(config)
        key = hashlib.sha256(text.encode()).hexdigest()[:16]
        directory = BUILD / f"linked-{key}"
        runtime = directory / "runtime"
        runtime.mkdir(parents=True, exist_ok=True)
        (runtime / "ccfull.inp").write_text(text)
        interactive = "    Call Output_information()"      # the interactive setup of the main program
        core = _ccfull_without_main()
        if interactive not in core:
            raise RuntimeError(f"no '{interactive.strip()}' in {SOURCE}")
        source = (core.replace(interactive, "    Call CCBROM_Noninteractive_Setup()", 1)
                  + (SOLVER / "library_interface.f90").read_text())
        path = _compile(source, directory / "ccfull_ccbrom_linked.f90",
                        directory / f"libccbrom_{key}_{machine()}.dylib",
                        ["-O3", "-fPIC", "-shared", "-ffree-line-length-none"])
        try:
            self._lib = ctypes.CDLL(str(path))
        except OSError as error:
            raise RuntimeError(f"cannot load {path} into this {machine()} Python") from error
        for name, argtypes in _SIGNATURES.items():
            function = getattr(self._lib, name)
            function.argtypes, function.restype = argtypes, None
        status, nlevel, nradial = ctypes.c_int(), ctypes.c_int(), ctypes.c_int()
        rmin, dr = ctypes.c_double(), ctypes.c_double()
        previous = Path.cwd()
        os.chdir(runtime)                                     # CCFULL reads ccfull.inp from here
        try:
            self._lib.ccbrom_initialize(ctypes.byref(status), ctypes.byref(nlevel), ctypes.byref(nradial),
                                        ctypes.byref(rmin), ctypes.byref(dr))
        finally:
            os.chdir(previous)
        _check(status, "initialization")
        self.nlevel = nlevel.value

    def fields(self, point, radii) -> np.ndarray:
        """Coupling matrices (radius, channel, channel), the bare nuclear potential included."""
        beta2, beta4, beta6, v0, r0, a0 = map(float, point)
        radii = np.ascontiguousarray(radii, dtype=np.float64)
        values = np.empty(radii.size * self.nlevel * self.nlevel, dtype=np.float64)
        status = ctypes.c_int()
        self._lib.ccbrom_evaluate_fields(beta2, beta4, beta6, v0, r0, a0, radii.size, radii,
                                         values, ctypes.byref(status))
        _check(status, "coupling field")
        fields = values.reshape(radii.size, self.nlevel, self.nlevel)
        size = float(self.config["apro"])**(1.0/3.0) + float(self.config["atar"])**(1.0/3.0)
        nuclear = -v0/(1.0 + np.exp((radii-r0*size)/a0))
        diagonal = np.arange(self.nlevel)
        fields[:, diagonal, diagonal] += nuclear[:, None]
        return fields

    def factors(self, point, radii, decomposition_inverse) -> np.ndarray:
        """Multipole form factors (radius, multipole) of the coupling field, nuclear part included.

        ``decomposition_inverse`` (multipole, channel**2) maps a row-major coupling
        matrix to the form factors; Fortran applies it, so no matrix is returned.
        """
        beta2, beta4, beta6, v0, r0, a0 = map(float, point)
        radii = np.ascontiguousarray(radii, dtype=np.float64)
        inverse = np.ascontiguousarray(decomposition_inverse, dtype=np.float64)
        if inverse.ndim != 2 or inverse.shape[1] != self.nlevel * self.nlevel:
            raise ValueError("the multipole decomposition does not match the CCFULL channels")
        values = np.empty((radii.size, inverse.shape[0]), dtype=np.float64)
        status = ctypes.c_int()
        self._lib.ccbrom_evaluate_factors(beta2, beta4, beta6, v0, r0, a0, radii.size, radii,
                                          inverse.shape[0], inverse.reshape(-1), 1,
                                          values.reshape(-1), ctypes.byref(status))
        _check(status, "multipole factors")
        return values

    def outer_map(self, energy: float, angular_momentum: int, potentials: np.ndarray,
                  radii: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """CCFULL's two-point Coulomb matching of one (E, J) at r1 = r_max - 2 dr and r2 = r_max.

        ``potentials`` are the channel potentials (2, channel) at ``radii`` = (r1, r2).
        Returns the diagonal ratio matrix and the source vector of the entrance channel.
        """
        r1, r2 = radii
        ec1 = float(energy) - potentials[0]
        ec2 = float(energy) - potentials[1]
        ratio = np.empty(ec1.size, dtype=np.complex128)
        source = np.zeros(ec1.size, dtype=np.complex128)
        open_mask = (ec1 > 0.0) & (ec2 > 0.0)
        closed = ~open_mask
        apro, atar = float(self.config["apro"]), float(self.config["atar"])
        reduced_mass = apro*atar/(apro+atar)*NUCLEON_MASS
        ratio[closed] = np.exp(-np.sqrt(2.0*reduced_mass/HBARC**2*np.abs(ec1[closed]))*(r2-r1))
        if np.any(open_mask):
            channels = np.flatnonzero(open_mask)
            energies = np.ascontiguousarray(
                np.column_stack((ec1[channels], ec2[channels])).reshape(-1), dtype=np.float64)
            points = np.ascontiguousarray(np.tile(np.array([r1, r2], dtype=np.float64), channels.size))
            up_real, up_imag, down_real, down_imag = (np.empty_like(energies) for _ in range(4))
            status = ctypes.c_int()
            self._lib.ccbrom_coulomb_values(int(angular_momentum), energies.size, energies, points,
                                            up_real, up_imag, down_real, down_imag, ctypes.byref(status))
            _check(status, "Coulomb functions")
            up = (up_real + 1j*up_imag).reshape(-1, 2)
            down = (down_real + 1j*down_imag).reshape(-1, 2)
            for local, channel in enumerate(channels):
                ratio[channel] = up[local, 1]/up[local, 0]
                if channel == 0:
                    source[channel] = down[local, 1] - ratio[channel]*down[local, 0]
        return np.diag(ratio), source

    def partial_wave_count(self, energy: float, v0: float, r0: float, a0: float, lmax: int = 200) -> int:
        """The number of partial waves J = 0, 1, ... that CCFULL main solves at ``energy``."""
        count = ctypes.c_int()
        self._lib.ccbrom_active_partial_waves(float(energy), float(v0), float(r0), float(a0),
                                              int(lmax), ctypes.byref(count))
        return count.value

    def penetrabilities(self, point, energy: float, waves) -> np.ndarray:
        """Exact Numerov P_J at one point and energy, with CCFULL's own boundary conditions."""
        beta2, beta4, beta6, v0, r0, a0 = map(float, point)
        waves = np.ascontiguousarray(waves, dtype=np.int32).reshape(-1)
        values = np.empty(waves.size, dtype=np.float64)
        status = ctypes.c_int()
        self._lib.ccbrom_exact_penetrabilities(beta2, beta4, beta6, v0, r0, a0, float(energy),
                                               waves.size, waves, values, ctypes.byref(status))
        _check(status, "exact Numerov")
        return values


_LIBRARIES: dict[str, Library] = {}


def library(config: dict) -> Library:
    """The CCFULL library of ``config``, built if needed and loaded once per process."""
    text = input_text(config)
    if text not in _LIBRARIES:
        _LIBRARIES[text] = Library(config)
    return _LIBRARIES[text]
