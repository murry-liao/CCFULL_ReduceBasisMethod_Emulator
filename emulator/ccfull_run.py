"""Everything the emulator does with files: CCFULL's input and output, the emulator's own
input files, and the experimental data.

The two CCFULL programs are in physical_model/ (build each with make):
    ccfull_clean       the exact model; it writes sigma.dat
    ccfull_interface   the same program, which also writes the wave functions (the snapshots)
The emulator only changes ccfull.inp: the six parameters on lines 3 and 7, the energy on
line 8, and a fixed r_min on line 10, so that every parameter point has the same grid.
"""

from pathlib import Path
import subprocess
import tempfile

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
CLEAN = PROJECT / "physical_model" / "ccfull_clean" / "ccfull"
INTERFACE = PROJECT / "physical_model" / "ccfull_interface" / "ccfull"


# ----------------------------------------------------------------- the emulator's input files

def data_lines(path):
    """The lines of an input file that hold data: blank lines and lines starting with # are
    skipped, and nothing after the dashed line (the explanation) is read."""
    lines = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith("-----"):
            break
        if line and not line.startswith("#"):
            lines.append(line)
    return lines


def numbers(line):
    """The numbers of an input line, separated by commas or blanks."""
    return [float(item) for item in line.replace(",", " ").split()]


def project_path(path):
    """Paths in the input files are relative to the project root."""
    return (PROJECT / Path(path).expanduser()).resolve()


def read_experimental(path):
    """A file of measured fusion cross sections: E_cm (MeV), sigma (mb), the lower and upper
    error (mb), and E_lab, one energy per line.  Returns the first four columns."""
    rows = [numbers(line) for line in data_lines(path)]
    energy, sigma, lower, upper = (np.array([row[k] for row in rows]) for k in range(4))
    return energy, sigma, lower, upper


# ----------------------------------------------------------------- CCFULL input and output

def ccfull_input(records, theta, energy, r_min, last_j=None):
    """ccfull.inp for one parameter point and one energy.

    records are the nine CCFULL records (lines 1-9 of emulator_train.in) of a rotational
    target; line 3 is E2T, beta2, beta4, beta6, NrotT and line 7 is V0, r0, a0.  With
    r_min None there is no line 10, and CCFULL puts r_min at the bottom of the pocket.
    last_j (line 11, read by the interface program only): after the J sum has ended, the
    partial waves up to last_j are solved as well, for their wave functions.
    """
    beta2, beta4, beta6, v0, r0, a0 = (float(value) for value in theta)
    lines = list(records)
    e2, nrot = lines[2].replace(",", " ").split()[0], lines[2].replace(",", " ").split()[-1]
    lines[2] = f"{e2},{beta2!r},{beta4!r},{beta6!r},{nrot}"
    lines[6] = f"{v0!r},{r0!r},{a0!r}"
    lines[7] = f"{float(energy)!r},{float(energy)!r},1.0"
    if r_min is not None:
        lines.append(repr(float(r_min)))              # line 10: the fixed r_min
        if last_j is not None:
            lines.append(str(int(last_j)))            # line 11: the last J to solve and write
    return "\n".join(lines) + "\n"


def run_ccfull(program, text, work):
    """Run one CCFULL program on the input text, in the directory work."""
    if not program.exists():
        raise FileNotFoundError(f"{program} is not built: run make in {program.parent}")
    (work / "ccfull.inp").write_text(text)
    # CCFULL asks a few questions while it starts (generalised couplings?); every answer is "n".
    done = subprocess.run([str(program)], cwd=work, input="n\n" * 64, text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if done.returncode != 0:
        raise RuntimeError(f"{program} stopped:\n{done.stdout[-3000:]}")


def exact_cross_sections(records, theta, energies, r_min):
    """sigma (mb) of the clean CCFULL program at one parameter point: one run per energy,
    because the measured energies are not evenly spaced."""
    sigma = []
    with tempfile.TemporaryDirectory(prefix="ccfull-") as folder:
        for energy in energies:
            run_ccfull(CLEAN, ccfull_input(records, theta, energy, r_min), Path(folder))
            sigma.append(np.loadtxt(Path(folder) / "sigma.dat", skiprows=1)[1])
    return np.array(sigma)


def wave_functions(records, theta, energy, r_min, last_j=None):
    """Run the interface program at one parameter point and energy.

    Returns a dict: J and P_J of every partial wave solved (those of the J sum, and up to
    last_j if given), their wave functions u (J, radius, channel), the matching table (one row
    per J and channel: E, J, channel, the channel energies at r_max - dr and r_max + dr, F and
    G at both), and the grid.
    """
    with tempfile.TemporaryDirectory(prefix="ccfull-") as folder:
        work = Path(folder)
        run_ccfull(INTERFACE, ccfull_input(records, theta, energy, r_min, last_j), work)
        lines = [line for line in (work / "waves_info.txt").read_text().splitlines()
                 if not line.startswith("#")]
        channels, points = (int(value) for value in lines[0].split()[:2])
        r_min, dr, reduced_mass = (float(value) for value in lines[0].split()[2:5])
        table = np.array([numbers(line) for line in lines[1:]])
        waves = np.fromfile(work / "waves.bin", dtype=np.complex128)
        matching = np.loadtxt(work / "matching.txt", ndmin=2)
    return {"J": table[:, 1].astype(int), "P": table[:, 2],
            "waves": waves.reshape(len(table), points, channels), "matching": matching,
            "radius": r_min + dr * np.arange(points), "reduced_mass": reduced_mass}
