"""CCFULL input records and experimental fusion cross-section files."""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
import re


def _numbers(line: str) -> list[float]:
    return [float(item) for item in re.split(r"[,\s]+", line.strip()) if item]


def parse_input(path: str | Path) -> dict:
    """The nine CCFULL records of ``path``; comment lines and anything after them are ignored."""
    lines = [line.strip() for line in Path(path).read_text().splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    if len(lines) < 9:
        raise ValueError(f"{path}: a CCFULL input needs nine data records")
    records = [_numbers(line) for line in lines[:9]]
    if any(len(record) < size for record, size in zip(records, (4, 4, 4, 4, 4, 2, 3, 3, 2))):
        raise ValueError(f"{path}: a CCFULL record has too few values")
    reaction, geometry, target1, target2, projectile, transfer, potential, energy, radial = records
    target_mode = int(geometry[3])
    if target_mode == 1:
        # Rotor record E2T,beta2,beta4,beta6,NrotT; the original four-value record has no beta6.
        target1 = target1[:5] if len(target1) >= 5 else [*target1[:3], 0.0, target1[3]]
    else:
        target1 = target1[:4]
    return {
        "apro": int(reaction[0]), "zpro": int(reaction[1]),
        "atar": int(reaction[2]), "ztar": int(reaction[3]),
        "r0p": geometry[0], "projectile_mode": int(geometry[1]),
        "r0t": geometry[2], "target_mode": target_mode,
        "target1": target1, "target2": target2[:4], "projectile": projectile[:4],
        "qtrans": transfer[0], "ftr": transfer[1],
        "v0": potential[0], "r0": potential[1], "a0": potential[2],
        "emin": energy[0], "emax": energy[1], "de": energy[2],
        "rmax": radial[0], "dr": radial[1],
    }


def _validate(config: dict) -> None:
    if min(int(config[key]) for key in ("apro", "zpro", "atar", "ztar")) <= 0:
        raise ValueError("mass and charge numbers must be positive")
    if int(config["zpro"]) > int(config["apro"]) or int(config["ztar"]) > int(config["atar"]):
        raise ValueError("a charge number cannot exceed its mass number")
    if (int(config["projectile_mode"]) not in (-1, 0, 1)
            or int(config["target_mode"]) not in (-1, 0, 1)):
        raise ValueError("a coupling mode must be -1 (inert), 0 (vibration), or 1 (rotation)")
    for key in ("r0p", "r0t", "v0", "r0", "a0", "de", "rmax", "dr"):
        if float(config[key]) <= 0:
            raise ValueError(f"{key} must be positive")
    if float(config["emax"]) < float(config["emin"]):
        raise ValueError("emax must be at least emin")


def input_text(config: dict) -> str:
    """The nine records of ``config`` as CCFULL reads them from ccfull.inp."""
    _validate(config)
    rows = [
        [int(config["apro"]), int(config["zpro"]), int(config["atar"]), int(config["ztar"])],
        [config["r0p"], int(config["projectile_mode"]), config["r0t"], int(config["target_mode"])],
        config["target1"], config["target2"], config["projectile"],
        [config["qtrans"], config["ftr"]],
        [config["v0"], config["r0"], config["a0"]],
        [config["emin"], config["emax"], config["de"]],
        [config["rmax"], config["dr"]],
    ]
    return "\n".join(
        ",".join(str(value) if isinstance(value, int) else f"{float(value):.12g}" for value in row)
        for row in rows) + "\n"


def records_sha256(path: str | Path) -> str:
    """SHA-256 of the nine CCFULL records of ``path``; comments and the training section are excluded."""
    return hashlib.sha256(input_text(parse_input(path)).encode()).hexdigest()


def energy_values(config: dict) -> list[float]:
    """Emin, Emin + dE, ... up to Emax (record 8)."""
    first, last, step = map(float, (config["emin"], config["emax"], config["de"]))
    if step <= 0 or last < first:
        raise ValueError("the CCFULL energy grid needs dE > 0 and Emax >= Emin")
    count = int(math.floor((last - first) / step + 1e-9))
    return [first + index * step for index in range(count + 1)]


_DATA_NAME = re.compile(r"^fusion-exp-(\d+)\.(\d+)\.(\d+)\.(\d+)-[^-]+-\d+-\d+\.txt$")


def data_nuclei(path: str | Path) -> tuple[tuple[int, int], tuple[int, int]]:
    """(Z, A) of both nuclei of an Expdata file, named fusion-exp-Z1.A1.Z2.A2-kind-year-set.txt."""
    match = _DATA_NAME.match(Path(path).name)
    if not match:
        raise ValueError(f"unrecognized Expdata file name: {Path(path).name}")
    z1, a1, z2, a2 = map(int, match.groups())
    return (z1, a1), (z2, a2)


def load_dataset(path: str | Path) -> list[dict]:
    """The rows of an Expdata file: E_cm, sigma_fus, and its lower and upper uncertainty."""
    rows = []
    for number, line in enumerate(Path(path).read_text().splitlines(), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 5:
            raise ValueError(f"{path}:{number}: expected five columns")
        energy, cross_section, lower, upper, _ = map(float, fields)     # the fifth is E_lab
        rows.append({"energy_cm_mev": energy, "fusion_cross_section_mb": cross_section,
                     "uncertainty_lower_mb": lower, "uncertainty_upper_mb": upper})
    if not rows:
        raise ValueError(f"experimental dataset is empty: {path}")
    return rows
