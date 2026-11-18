"""Bayesian calibration of the six CCFULL parameters with a trained reduced model.

    python -m emulator.mcmc --config emulator/emulator_mcmc_34S_168Er.in      # a new run
    python -m emulator.mcmc --config <file> --check                           # only prints the input
    python -m emulator.mcmc --config <file> --until best_fit                  # stops after a stage
    python -m emulator.mcmc --resume emulator/outputs/<run folder>            # continues a run

The lines of the input file are explained at its end.  The stages run in order, into one run
folder; a stage that is already finished there is skipped, so --resume continues a run:
  1. best_fit        a Latin-hypercube scan of the training box, then bounded least squares
                     from the best scan points and from the nominal point of the model;
  2. pilot           a short ensemble chain around the best fit, with the experimental errors;
  3. emulator_error  exact CCFULL at points of the pilot chain: the error of the emulator at
                     every measured energy, added in quadrature to the experimental error;
  4. production      the ensemble chain with both errors, saved every 50 steps;
  5. check           exact CCFULL at points of the posterior: at the measured energies, and on an
                     energy grid (the excitation function and the barrier distribution).
The prior is uniform in the training box and the likelihood is Gaussian in the measured fusion
cross sections.  The chains use the affine-invariant ensemble sampler of Goodman and Weare, as
implemented in emcee.
"""

import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import shutil
import time

for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(variable, "1")

import emcee
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from emulator import rom
from emulator.ccfull_run import data_lines, exact_cross_sections, numbers, project_path, read_experimental


STAGES = ("best_fit", "pilot", "emulator_error", "production", "check")
SAVE_EVERY = 50          # chain steps between two saves of a chain file
EDGE = 0.02              # a best-fit parameter this close to a bound, in units of its range, is at the edge
DIFF_STEP = 1.0e-3       # finite-difference step of the least squares, in units of the box width
FIT_EVALUATIONS = 100    # cross-section evaluations of one least-squares fit, besides its Jacobians
SPREAD = 0.05            # largest standard deviation of the starting walkers, in units of the box width


# ----------------------------------------------------------------- small helpers

def log(folder, message):
    """Print a line and append it to progress.log of the run."""
    line = f"{datetime.now().strftime('%H:%M:%S')} {message}"
    print(line, flush=True)
    with open(folder / "progress.log", "a") as stream:
        stream.write(line + "\n")


def update_summary(folder, **items):
    """Add items to mcmc_summary.json of the run."""
    file = folder / "mcmc_summary.json"
    summary = json.loads(file.read_text()) if file.exists() else {}
    summary.update(items)
    file.write_text(json.dumps(summary, indent=2) + "\n")


def named(theta):
    """The six values by name; an undefined value (a chain too short for its autocorrelation) is None."""
    return {name: float(value) if np.isfinite(value) else None for name, value in zip(rom.PARAMETER_NAMES, theta)}


def log_posterior(theta, model, energies, data, inverse_error, box):
    """emcee's function: the log posterior and, kept as the blob of the sample, the cross sections."""
    theta = np.asarray(theta, dtype=float)
    nothing = np.full(len(energies), np.nan)
    if np.any(theta < box[:, 0]) or np.any(theta > box[:, 1]):
        return -np.inf, nothing
    try:
        sigma = rom.cross_sections(model, theta, energies)
    except np.linalg.LinAlgError:
        return -np.inf, nothing
    chi2 = float(np.sum(((sigma - data) * inverse_error) ** 2))
    return (-0.5 * chi2 if np.isfinite(chi2) else -np.inf), sigma


def autocorrelation_time(chain):
    """Integrated autocorrelation time of each parameter, in steps (emcee's estimate, never refused)."""
    with np.errstate(invalid="ignore", divide="ignore"):       # a chain of a few steps has no estimate
        return emcee.autocorr.integrated_time(chain, tol=0, quiet=True)


def run_chain(start, steps, seed, file, folder, stage, posterior_args):
    """An ensemble chain of steps steps, saved to file every SAVE_EVERY steps.  A chain already in
    file is continued from its last state and random state, with the same result."""
    previous = None
    if file.exists():
        with np.load(file) as saved:
            previous = {key: saved[key] for key in saved.files}
        done = previous["chain"].shape[0]
        if done >= steps:
            return previous
        state = emcee.State(previous["chain"][-1], log_prob=previous["log_prob"][-1], blobs=previous["blobs"][-1],
                            random_state=("MT19937", previous["rng_keys"], int(previous["rng_position"]),
                                          int(previous["rng_has_gauss"]), float(previous["rng_cached_gaussian"])))
        log(folder, f"[{stage}] continuing the chain at step {done}")
    else:
        done = 0
        state = emcee.State(np.asarray(start, dtype=float), random_state=np.random.RandomState(seed).get_state())
    sampler = emcee.EnsembleSampler(len(state.coords), len(rom.PARAMETER_NAMES), log_posterior, args=posterior_args)
    started, first = time.perf_counter(), done
    while done < steps:
        count = min(SAVE_EVERY, steps - done)
        state = sampler.run_mcmc(state, count)
        done += count
        arrays = {"chain": sampler.get_chain(), "log_prob": sampler.get_log_prob(), "blobs": sampler.get_blobs(),
                  "accepted": sampler.backend.accepted.astype(float)}
        if previous is not None:
            for key in ("chain", "log_prob", "blobs"):
                arrays[key] = np.concatenate([previous[key], arrays[key]])
            arrays["accepted"] = arrays["accepted"] + previous["accepted"]
        name, keys, position, has_gauss, cached = sampler.random_state
        arrays.update(rng_keys=np.asarray(keys), rng_position=position, rng_has_gauss=has_gauss,
                      rng_cached_gaussian=cached, energies=posterior_args[1])
        np.savez(file.with_suffix(".tmp.npz"), **arrays)
        file.with_suffix(".tmp.npz").replace(file)
        tau = autocorrelation_time(arrays["chain"])
        log(folder, f"[{stage}] step {done}/{steps}, {(time.perf_counter() - started) / (done - first):.1f} s per step, "
                    f"acceptance {np.mean(arrays['accepted']) / done:.2f}, autocorrelation time "
                    f"{tau.min():.0f}-{tau.max():.0f} steps, lowest chi2 {-2.0 * float(np.max(arrays['log_prob'])):.1f}")
    return arrays


# ----------------------------------------------------------------- the stages

def best_fit(folder, model, setup):
    """Scan the box, then least squares from the nominal point and from the best scan points."""
    box, energies, data, errors = setup["box"], setup["energies"], setup["data"], setup["errors"]
    lower, width = box[:, 0], box[:, 1] - box[:, 0]
    dof = len(energies) - 6
    evaluations = [0]
    inverse_error = 1.0 / np.sqrt(errors ** 2)               # the experimental errors only

    def residuals(theta):
        evaluations[0] += 1
        return (rom.cross_sections(model, theta, energies) - data) * inverse_error

    started = time.perf_counter()
    scan = rom.latin_hypercube(box, setup["scan_points"], setup["scan_seed"])
    chi2 = np.empty(len(scan))
    for k, theta in enumerate(scan):
        chi2[k] = float(np.sum(residuals(theta) ** 2))
        if (k + 1) % max(1, len(scan) // 10) == 0 or k + 1 == len(scan):
            log(folder, f"[best_fit] scan {k + 1}/{len(scan)}, lowest chi2/dof {chi2[:k + 1].min() / dof:.2f}")
    pd.DataFrame([{"scan_index": k + 1, **named(theta), "chi2": value} for k, (theta, value) in enumerate(zip(scan, chi2))]
                 ).to_csv(folder / "scan.csv", index=False)

    starts = [("CCFULL input", setup["nominal"])] + [(f"scan point {k + 1}", scan[k])
                                                      for k in np.argsort(chi2)[:setup["local_fits"]]]
    fits = []
    for label, start in starts:
        calls, clock = evaluations[0], time.perf_counter()
        result = least_squares(lambda unit: residuals(lower + width * np.clip(unit, 0.0, 1.0)),
                               np.clip((start - lower) / width, 0.0, 1.0), bounds=(0.0, 1.0), diff_step=DIFF_STEP,
                               max_nfev=FIT_EVALUATIONS, ftol=1.0e-6, xtol=1.0e-6, gtol=1.0e-6)
        covariance = np.linalg.pinv(result.jac.T @ result.jac)    # Laplace approximation, in units of the box
        edge = [name for name, value in zip(rom.PARAMETER_NAMES, result.x) if value < EDGE or value > 1.0 - EDGE]
        fit = {"start": label, "start_point": named(start), "point": named(lower + width * np.clip(result.x, 0.0, 1.0)),
               "unit_point": [float(value) for value in result.x], "chi2": 2.0 * float(result.cost),
               "chi2_per_dof": 2.0 * float(result.cost) / dof,
               "laplace_sd": named(np.sqrt(np.clip(np.diag(covariance), 0.0, None)) * width),
               "at_edge": edge, "evaluations": evaluations[0] - calls, "seconds": time.perf_counter() - clock,
               "status": int(result.status), "message": str(result.message), "covariance_unit": covariance.tolist()}
        fits.append(fit)
        log(folder, f"[best_fit] from {label}: chi2/dof {fit['chi2_per_dof']:.2f} after {fit['evaluations']} evaluations "
                    f"({fit['seconds']:.0f} s); " + ", ".join(f"{n} {v:.5g}" for n, v in fit["point"].items())
            + (f"; at the edge: {', '.join(edge)}" if edge else ""))
    pd.DataFrame([{"start": f["start"], **{f"start_{k}": v for k, v in f["start_point"].items()}, **f["point"],
                   "chi2": f["chi2"], "chi2_per_dof": f["chi2_per_dof"],
                   **{f"laplace_sd_{k}": v for k, v in f["laplace_sd"].items()}, "at_edge": " ".join(f["at_edge"]),
                   "evaluations": f["evaluations"], "status": f["status"]} for f in fits]
                 ).to_csv(folder / "best_fit.csv", index=False)
    best = min(fits, key=lambda f: f["chi2"])
    np.savez(folder / "best_fit.npz", unit=np.array(best["unit_point"]), point=np.array(list(best["point"].values())),
             covariance=np.array(best["covariance_unit"]))
    update_summary(folder, best_fit={
        "degrees_of_freedom": dof, "seconds": time.perf_counter() - started,
        "scan": {"points": len(scan), "lowest_chi2": float(chi2.min())},
        "fits": [{k: v for k, v in f.items() if k != "covariance_unit"} for f in fits],
        "best": {k: best[k] for k in ("start", "point", "chi2", "chi2_per_dof", "laplace_sd", "at_edge")}})
    if best["at_edge"]:
        log(folder, f"[best_fit] the best fit is at the edge of the training box in {', '.join(best['at_edge'])}")


def pilot(folder, model, setup):
    """Walkers drawn around the best fit (its Laplace approximation, each standard deviation at most
    SPREAD of the box and reflected into it), then a short chain with the experimental errors."""
    box = setup["box"]
    with np.load(folder / "best_fit.npz") as saved:
        center, covariance = saved["unit"], saved["covariance"]
    values, vectors = np.linalg.eigh(0.5 * (covariance + covariance.T))
    scale = vectors * np.sqrt(np.clip(values, 1.0e-14, SPREAD ** 2))
    draws = np.abs(center + np.random.default_rng(setup["seed"]).standard_normal((setup["walkers"], 6)) @ scale.T)
    draws = np.where(draws > 1.0, 2.0 - draws, draws)
    start = box[:, 0] + (box[:, 1] - box[:, 0]) * draws
    started = time.perf_counter()
    args = (model, setup["energies"], setup["data"], 1.0 / np.sqrt(setup["errors"] ** 2), box)
    arrays = run_chain(start, setup["pilot_steps"], setup["seed"], folder / "pilot_chain.npz", folder, "pilot", args)
    update_summary(folder, pilot={"walkers": setup["walkers"], "steps": int(arrays["chain"].shape[0]),
                                  "acceptance": float(np.mean(arrays["accepted"]) / arrays["chain"].shape[0]),
                                  "autocorrelation_time": named(autocorrelation_time(arrays["chain"])),
                                  "seconds": time.perf_counter() - started})


def emulator_error(folder, model, setup):
    """RMS of sigma_ROM - sigma_exact at each measured energy, over points of the second half of the
    pilot chain."""
    energies, data, errors = setup["energies"], setup["data"], setup["errors"]
    started = time.perf_counter()
    with np.load(folder / "pilot_chain.npz") as saved:
        half = saved["chain"].shape[0] // 2
        points = saved["chain"][half:].reshape(-1, 6)
        blobs = saved["blobs"][half:].reshape(-1, len(energies))
    points, first = np.unique(points, axis=0, return_index=True)        # rejected moves repeat points
    blobs = blobs[first]
    picks = np.random.default_rng(setup["seed"] + 1).choice(len(points), size=min(setup["error_points"], len(points)),
                                                            replace=False)
    exact, rows = [], []
    for done, pick in enumerate(picks, start=1):
        exact.append(exact_cross_sections(model["records"], points[pick], energies, model["radius"][0]))
        rows.extend({"point_index": done, **named(points[pick]), "energy_cm_mev": float(e),
                     "emulator_cross_section_mb": float(a), "exact_cross_section_mb": float(b)}
                    for e, a, b in zip(energies, blobs[pick], exact[-1]))
        if done % 5 == 0 or done == len(picks):
            log(folder, f"[emulator_error] exact CCFULL at {done}/{len(picks)} pilot points")
    if len(picks):
        difference = blobs[picks] - np.array(exact)
        error = np.sqrt(np.mean(difference ** 2, axis=0))
        relative = np.median(np.abs(difference) / np.abs(np.array(exact)), axis=0)
        worst = np.max(np.abs(difference), axis=0) / errors
        pd.DataFrame(rows).to_csv(folder / "emulator_error_points.csv", index=False)
    else:
        error = relative = worst = np.zeros(len(energies))
    pd.DataFrame({"energy_cm_mev": energies, "experimental_cross_section_mb": data, "experimental_error_mb": errors,
                  "emulator_error_mb": error, "emulator_over_experimental_error": error / errors,
                  "median_relative_error": relative, "max_error_over_experimental_error": worst}
                 ).to_csv(folder / "emulator_error.csv", index=False)
    ratio = error / errors
    log(folder, f"[emulator_error] emulator / experimental error: median {np.median(ratio):.2f}, max {ratio.max():.2f}")
    update_summary(folder, emulator_error={"points": int(len(picks)), "emulator_over_experimental_error": {
        "median": float(np.median(ratio)), "max": float(ratio.max())}, "seconds": time.perf_counter() - started})


def production(folder, model, setup):
    """The chain with both errors, continued from the last state of the pilot chain."""
    with np.load(folder / "pilot_chain.npz") as saved:
        start = saved["chain"][-1]
    emulator = pd.read_csv(folder / "emulator_error.csv")["emulator_error_mb"].to_numpy()
    inverse_error = 1.0 / np.sqrt(setup["errors"] ** 2 + emulator ** 2)
    started = time.perf_counter()
    args = (model, setup["energies"], setup["data"], inverse_error, setup["box"])
    arrays = run_chain(start, setup["production_steps"], setup["seed"] + 2, folder / "chain.npz", folder,
                       "production", args)
    update_summary(folder, production={"walkers": setup["walkers"], "steps": int(arrays["chain"].shape[0]),
                                       "acceptance": float(np.mean(arrays["accepted"]) / arrays["chain"].shape[0]),
                                       "autocorrelation_time": named(autocorrelation_time(arrays["chain"])),
                                       "seconds_this_session": time.perf_counter() - started})


def check(folder, model, setup):
    """Posterior samples, and exact CCFULL at the best sample and at points drawn from the posterior."""
    energies, data, errors, grid = setup["energies"], setup["data"], setup["errors"], setup["grid"]
    started = time.perf_counter()
    with np.load(folder / "chain.npz") as saved:
        chain, log_prob = saved["chain"], saved["log_prob"]
    tau = autocorrelation_time(chain)
    tau = np.where(np.isfinite(tau) & (tau > 0.0), tau, 1.0)
    burn = min(int(max(SAVE_EVERY, np.ceil(2.0 * tau.max()))), chain.shape[0] // 2)
    thin = max(1, int(tau.min() // 2))
    samples = chain[burn::thin].reshape(-1, 6)
    pd.DataFrame({**{name: samples[:, k] for k, name in enumerate(rom.PARAMETER_NAMES)},
                  "log_posterior": log_prob[burn::thin].reshape(-1)}).to_csv(folder / "posterior_samples.csv", index=False)
    best = np.unravel_index(int(np.argmax(log_prob)), log_prob.shape)
    picks = np.random.default_rng(setup["seed"] + 3).choice(len(samples), size=min(setup["check_points"], len(samples)),
                                                            replace=False)
    points = [("best sample", chain[best])] + [("posterior", samples[pick]) for pick in picks]
    at_data, on_grid = [], []
    for k, (label, theta) in enumerate(points):
        sigma_rom = rom.cross_sections(model, theta, energies)
        sigma_exact = exact_cross_sections(model["records"], theta, energies, model["radius"][0])
        sigma_grid = exact_cross_sections(model["records"], theta, grid, model["radius"][0])
        at_data.extend({"draw": k, "kind": label, **named(theta), "energy_cm_mev": float(e),
                        "emulator_cross_section_mb": float(a), "exact_cross_section_mb": float(b),
                        "experimental_cross_section_mb": float(c), "experimental_error_mb": float(d)}
                       for e, a, b, c, d in zip(energies, sigma_rom, sigma_exact, data, errors))
        on_grid.extend({"draw": k, "kind": label, **named(theta), "energy_cm_mev": float(e),
                        "exact_cross_section_mb": float(value)} for e, value in zip(grid, sigma_grid))
        if k % 5 == 0 or k + 1 == len(points):
            log(folder, f"[check] exact CCFULL at {k + 1}/{len(points)} points")
    at_data, on_grid = pd.DataFrame(at_data), pd.DataFrame(on_grid)
    at_data.to_csv(folder / "posterior_check.csv", index=False)
    on_grid.to_csv(folder / "posterior_predictive.csv", index=False)
    deviation = (at_data["emulator_cross_section_mb"] - at_data["exact_cross_section_mb"]).abs() / at_data["experimental_error_mb"]
    plain = (data - at_data[at_data["draw"] == 0]["exact_cross_section_mb"].to_numpy()) / errors
    quantiles = np.percentile(samples, [2.5, 16.0, 50.0, 84.0, 97.5], axis=0)
    update_summary(folder, posterior={
        "burn_in_steps": burn, "thin": thin, "samples": int(len(samples)),
        "effective_samples": float(chain.shape[1] * (chain.shape[0] - burn) / tau.max()),
        "autocorrelation_time": named(tau),
        "quantiles_2.5_16_50_84_97.5": {name: quantiles[:, k].tolist() for k, name in enumerate(rom.PARAMETER_NAMES)},
        "correlation": np.corrcoef(samples.T).tolist(),
        "best_sample": {"point": named(chain[best]), "log_posterior": float(log_prob[best]),
                        "chi2_emulator_with_both_errors": -2.0 * float(log_prob[best]),
                        "chi2_exact_with_experimental_errors": float(np.sum(plain ** 2))},
    }, check={"points": len(points), "grid_mev": [float(grid[0]), float(grid[-1]), int(grid.size)],
              "emulator_minus_exact_over_experimental_error": {"median": float(deviation.median()),
                                                               "p95": float(deviation.quantile(0.95)),
                                                               "max": float(deviation.max())},
              "seconds": time.perf_counter() - started})
    log(folder, f"[check] emulator - exact at the posterior points, in experimental errors: median "
                f"{deviation.median():.2f}, 95th percentile {deviation.quantile(0.95):.2f}, max {deviation.max():.2f}")


# ----------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(description="Bayesian calibration with a trained reduced model of CCFULL.")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--config", help="input file; the run gets a new folder")
    source.add_argument("--resume", help="run folder to continue, with its copied input")
    parser.add_argument("--until", choices=STAGES, default=STAGES[-1], help="the last stage to run")
    parser.add_argument("--check", action="store_true", help="print how each line was read, and stop")
    args = parser.parse_args()

    # 1. Read the input (in a resumed run, the copy in its folder).
    if args.resume:
        folder = project_path(args.resume)
        inputs = sorted(folder.glob("*.in"))
        if len(inputs) != 1:
            raise SystemExit(f"{folder} should hold exactly one copied input file")
        config = inputs[0]
    else:
        config = Path(args.config).resolve()
    lines = data_lines(config)
    model = rom.load_model(project_path(lines[1]) / lines[0], int(numbers(lines[2])[0]))
    lowest, highest = numbers(lines[3])
    measured_energy, measured_sigma, lower_error, upper_error = read_experimental(project_path(model["experimental_file"]))
    keep = (measured_energy >= lowest) & (measured_energy <= highest)
    scan_points, scan_seed, local_fits = (int(value) for value in numbers(lines[4]))
    walkers, pilot_steps, production_steps, seed = (int(value) for value in numbers(lines[5]))
    check_values = numbers(lines[7])
    first, last, step = check_values[1:4]
    setup = {"box": model["box"], "nominal": rom.nominal_point(model), "energies": measured_energy[keep],
             "data": measured_sigma[keep],
             "errors": 0.5 * (lower_error[keep] + upper_error[keep]),     # the files hold equal lower and upper errors
             "scan_points": scan_points, "scan_seed": scan_seed, "local_fits": local_fits,
             "walkers": walkers, "pilot_steps": pilot_steps, "production_steps": production_steps, "seed": seed,
             "error_points": int(numbers(lines[6])[0]), "check_points": int(check_values[0]),
             "grid": np.round(first + step * np.arange(int(np.floor((last - first) / step + 1.0e-9)) + 1), 10)}
    if len(setup["energies"]) <= 6:
        raise SystemExit(f"line 4 keeps {len(setup['energies'])} measured cross sections; more than six are needed")
    if walkers < 12 or pilot_steps < 2:
        raise SystemExit("line 6: emcee needs at least 12 walkers, and the pilot chain at least two steps")
    system = model["system"]
    reaction = f"A={system['apro']}, Z={system['zpro']} + A={system['atar']}, Z={system['ztar']}"

    # 2. Say how each line was read.
    print(f"MCMC input {config}")
    print(f"  lines 1-3  {model['type']} model {model['name']} ({reaction}), rank {model['rank']} of {model['trained_rank']}")
    print(f"  line 4     {len(setup['energies'])} measured cross sections from {setup['energies'].min():g} to "
          f"{setup['energies'].max():g} MeV")
    print(f"  line 5     best fit: {scan_points} scan points (seed {scan_seed}), least squares from the best "
          f"{local_fits} and from the nominal point " + ", ".join(f"{n}={v:g}" for n, v in zip(rom.PARAMETER_NAMES, setup["nominal"])))
    print(f"  line 6     {walkers} walkers, {pilot_steps} pilot and {production_steps} production steps, seed {seed}")
    print(f"  line 7     emulator error from exact CCFULL at {setup['error_points']} pilot points")
    print(f"  line 8     check: {setup['check_points']} posterior points and the best sample; excitation functions "
          f"at {setup['grid'].size} energies from {setup['grid'][0]:g} to {setup['grid'][-1]:g} MeV")
    print(f"  line 9     output folder {project_path(lines[8])}")
    if args.check:
        return

    # 3. The run folder: new, with a copy of the input, or the resumed one.
    if not args.resume:
        created = datetime.now().astimezone()
        suffix = "" if model["rank"] == model["trained_rank"] else f"_r{model['rank']}"
        folder = project_path(lines[8]) / f"{created.strftime('%Y%m%d_%H%M%S_%f')}_{model['name']}{suffix}_mcmc"
        folder.mkdir(parents=True)
        shutil.copy2(config, folder / config.name)
        update_summary(folder, format="ccfull_rom_mcmc_v2", created_at=created.isoformat(), configuration=str(config),
                       model_id=model["name"], model_folder=str(model["folder"]), reaction_system=reaction,
                       model_type=model["type"], basis_rank=model["rank"], trained_rank=model["trained_rank"],
                       experimental_file=model["experimental_file"], energies_mev=setup["energies"].tolist(),
                       prior={"kind": "uniform in the training box",
                              "bounds": {n: list(map(float, b)) for n, b in zip(rom.PARAMETER_NAMES, model["box"])}},
                       nominal_point=named(setup["nominal"]))
    log(folder, f"[setup] {folder}")
    started = time.perf_counter()
    if model["type"] == "local":
        rom.load_blocks(model, setup["energies"])
    else:
        rom.prepare_global(model, setup["energies"])
    log(folder, f"[setup] {model['name']} at rank {model['rank']}, prepared in {time.perf_counter() - started:.0f} s")

    # 4. The stages, in order; a finished one is skipped.
    def steps_in(name):
        if not (folder / name).exists():
            return 0
        with np.load(folder / name) as saved:
            return saved["chain"].shape[0]
    finished = {"best_fit": lambda: (folder / "best_fit.npz").exists(),
                "pilot": lambda: steps_in("pilot_chain.npz") >= pilot_steps,
                "emulator_error": lambda: (folder / "emulator_error.csv").exists(),
                "production": lambda: steps_in("chain.npz") >= production_steps,
                "check": lambda: (folder / "posterior_predictive.csv").exists()}
    stage_function = {"best_fit": best_fit, "pilot": pilot, "emulator_error": emulator_error,
                      "production": production, "check": check}
    for stage in STAGES[:STAGES.index(args.until) + 1]:
        if finished[stage]():
            log(folder, f"[{stage}] already done")
            continue
        log(folder, f"[{stage}] starting")
        stage_function[stage](folder, model, setup)
        log(folder, f"[{stage}] done")
    log(folder, f"[complete] {folder}")


if __name__ == "__main__":
    main()
