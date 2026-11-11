"""CCFULL input records, line-based input files, and experimental fusion cross-section files."""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
import re


def numbers(line: str) -> list[float]:
    """The numbers of an input line, separated by commas or blanks."""
    return [float(item) for item in re.split(r"[,\s]+", line.strip()) if item]


def data_lines(path: str | Path) -> list[str]:
    """The lines of an input file that hold data: not blank and not starting with '#'."""
    lines = (line.strip() for line in Path(path).read_text().splitlines())
    return [line for line in lines if line and not line.startswith("#")]


def parse_input(path: str | Path) -> dict:
    """The nine CCFULL records at the top of ``path``; the lines after them are ignored."""
    (reaction, geometry, target1, target2, projectile, transfer, potential, energy,
     radial) = (numbers(line) for line in data_lines(path)[:9])
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


def input_text(config: dict) -> str:
    """The nine records of ``config`` as CCFULL reads them from ccfull.inp."""
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
    """SHA-256 of the nine CCFULL records of ``path``; comments and the lines after them are excluded."""
    return hashlib.sha256(input_text(parse_input(path)).encode()).hexdigest()


def reaction_label(config: dict) -> str:
    """For example "A=16, Z=8 + A=154, Z=62"."""
    return f"A={config['apro']}, Z={config['zpro']} + A={config['atar']}, Z={config['ztar']}"


def energy_values(config: dict) -> list[float]:
    """Emin, Emin + dE, ... up to Emax (record 8)."""
    first, last, step = map(float, (config["emin"], config["emax"], config["de"]))
    count = int(math.floor((last - first) / step + 1e-9))
    return [first + index * step for index in range(count + 1)]


def load_dataset(path: str | Path) -> list[dict]:
    """The rows of an Expdata file: E_cm, sigma_fus, its lower and upper uncertainty, and E_lab."""
    rows = []
    for line in data_lines(path):
        energy, cross_section, lower, upper, _ = numbers(line)
        rows.append({"energy_cm_mev": energy, "fusion_cross_section_mb": cross_section,
                     "uncertainty_lower_mb": lower, "uncertainty_upper_mb": upper})
    return rows
