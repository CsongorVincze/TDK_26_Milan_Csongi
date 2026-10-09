#!/usr/bin/env python3
"""Measure late-collapse scales and PDE-term balance for the Module 5 claim.

The solver is imported from Fraction_of_central/half_radius/study.py so the
experiment uses the established Module 4 equation and radial grid. Example:

    python3 run_audit.py --case 500:0.01 --case 1000:0.01 --case 2000:0.01

Each case is ``N:offset_percent[:CFL[:field_step]]``. The amplitude is the
same continuum Gaussian for every resolution, referenced to the existing
N=500 finite-horizon threshold; this is not a per-resolution critical search.
"""

import argparse
import csv
import importlib.util
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
SOLVER = HERE.parents[1] / "Fraction_of_central" / "half_radius" / "study.py"
A_REFERENCE = 1.202609700895846
LEVELS = (10, 15, 20, 30, 40, 60, 80, 120, 160, 240, 320, 480, 640, 960, 1280)
FIT_BANDS = ((20, 80), (20, 160), (40, 320), (80, 640),
             (160, 1280), (20, 320), (40, 640), (80, 1280), (20, 1280))
PROFILE_LEVELS = (20, 40, 80, 160, 320, 640, 1280)
PROFILE_RHO = np.linspace(0.0, 2.0, 101)


def load_solver():
    spec = importlib.util.spec_from_file_location("module5_audit_solver", SOLVER)
    solver = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = solver
    spec.loader.exec_module(solver)
    return solver


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def radial_norm(values, radius, limit):
    """Volume-weighted L2 norm on [0, limit] in spherical symmetry."""
    end = int(np.searchsorted(radius, limit, side="right"))
    r = radius[:end]
    v = values[:end]
    if len(r) < 2:
        return float("nan")
    if r[-1] < limit:
        r = np.append(r, limit)
        v = np.append(v, np.interp(limit, radius, values))
    return float(np.sqrt(np.trapezoid(4*np.pi*r*r*v*v, r)))


def fractional_radius(phi, radius, fraction):
    """First interpolated radius where |phi| falls to fraction of its center."""
    values = np.abs(phi)
    if values[0] <= 0.0 or not 0.0 < fraction < 1.0:
        return float("nan")
    target = fraction*values[0]
    crossings = np.flatnonzero(values[1:] <= target)
    if not len(crossings):
        return float("nan")
    i = int(crossings[0]+1)
    y0, y1 = values[i-1], values[i]
    if y0 == y1:
        return float(radius[i])
    frac = np.clip((target-y0)/(y1-y0), 0.0, 1.0)
    return float(radius[i-1]+frac*(radius[i]-radius[i-1]))


def profile_metrics(solver, phi, pi, grid, time, level, case):
    radius = grid.r
    r50 = solver.half_radius(phi, radius)
    if not np.isfinite(r50) or r50 <= 0:
        return None

    lap = solver.laplacian(phi, grid)
    mass = -phi
    cubic = phi**3
    acceleration = lap + mass + cubic - grid.damping*pi
    cubic_norm = radial_norm(cubic, radius, r50)
    grad = np.gradient(phi, radius, edge_order=2)
    gradient_at_r50 = float(np.interp(r50, radius, np.abs(grad)))
    center = abs(float(phi[0]))
    r25_fraction = fractional_radius(phi, radius, 0.25)
    r75_fraction = fractional_radius(phi, radius, 0.75)
    return {
        "case": case["name"], "N": case["n"], "offset_percent": case["offset"],
        "amplitude": case["amplitude"], "cfl": case["cfl"],
        "field_step": case["field_step"], "level": level,
        "time": time, "phi_center": float(phi[0]),
        "abs_phi_center": center, "r50": r50,
        "r25_fraction": r25_fraction, "r75_fraction": r75_fraction,
        "pi_center": float(pi[0]), "pi_over_phi2_center": float(pi[0])/center**2,
        "cells_per_r50": r50/grid.dr_min,
        "gradient_r50": gradient_at_r50,
        "shape_gradient_ratio": gradient_at_r50*r50/center,
        "lap_norm_over_cubic_norm": radial_norm(lap, radius, r50)/cubic_norm,
        "mass_norm_over_cubic_norm": radial_norm(mass, radius, r50)/cubic_norm,
        "rhs_norm_over_cubic_norm": radial_norm(acceleration, radius, r50)/cubic_norm,
        "lap_over_cubic_center": abs(float(lap[0]))/center**3,
        "rhs_over_cubic_center": abs(float(acceleration[0]))/center**3,
    }


def simulate(solver, case, max_level, capture_profile_levels=(), verbose=True):
    config = solver.Config(n=case["n"], cfl=case["cfl"])
    grid = solver.make_grid(config)
    phi = case["amplitude"]*np.exp(-grid.r**2/(2*config.r0**2))
    pi = np.zeros_like(phi)
    phi[-1] = pi[-1] = 0.0
    targets = [value for value in LEVELS if value <= max_level]
    next_target = 0
    time = 0.0
    steps = 0
    events = []
    profiles = []
    if verbose:
        print(f"{case['name']}: N={case['n']}, A={case['amplitude']:.12f}, "
              f"CFL={case['cfl']:g}, field-step={case['field_step']:g}; "
              f"target |phi(0)|={max(targets):g}", flush=True)

    while next_target < len(targets):
        center = abs(float(phi[0]))
        peak = max(float(np.max(np.abs(phi))), 1.0)
        dt = min(grid.dt, case["field_step"]/peak)
        new_phi, new_pi = solver.rk4(phi, pi, dt, grid)
        if not np.all(np.isfinite(new_phi)) or not np.all(np.isfinite(new_pi)):
            raise RuntimeError(f"Non-finite RK4 state in {case['name']} at t={time:g}.")
        next_center = abs(float(new_phi[0]))
        if next_center > center:
            while next_target < len(targets) and center < targets[next_target] <= next_center:
                level = targets[next_target]
                fraction = (level-center)/(next_center-center)
                event_phi = phi + fraction*(new_phi-phi)
                event_pi = pi + fraction*(new_pi-pi)
                event_phi[0] = np.copysign(level, new_phi[0])
                event_time = time + fraction*dt
                row = profile_metrics(solver, event_phi, event_pi, grid,
                                     event_time, level, case)
                if row:
                    events.append(row)
                    if level in capture_profile_levels:
                        normalized = np.interp(PROFILE_RHO*row["r50"], grid.r,
                                               np.abs(event_phi))/level
                        profiles.extend({
                            "case": case["name"], "N": case["n"],
                            "offset_percent": case["offset"],
                            "cfl": case["cfl"], "field_step": case["field_step"],
                            "level": level,
                            "rho_r_over_r50": float(rho),
                            "normalized_abs_phi": float(value),
                        } for rho, value in zip(PROFILE_RHO, normalized))
                    if verbose:
                        print(f"  |phi0|={level:g} at t={event_time:.9f}; "
                              f"R50/dr_min={row['cells_per_r50']:.1f}", flush=True)
                next_target += 1
        time += dt
        phi, pi = new_phi, new_pi
        steps += 1
        if verbose and steps % 5000 == 0:
            print(f"  progress t={time:.6f}, |phi0|={abs(phi[0]):.3g}, "
                  f"steps={steps}", flush=True)
        if time > 100.0:
            break

    return events, steps, time, profiles


def fit_profile_shapes(profile_rows):
    grouped = {}
    for row in profile_rows:
        grouped.setdefault((row["case"], row["level"]), []).append(row)
    output = []
    exponents = np.linspace(0.5, 4.0, 701)
    for (case, level), rows in grouped.items():
        rows.sort(key=lambda row: row["rho_r_over_r50"])
        rho = np.asarray([row["rho_r_over_r50"] for row in rows])
        profile = np.asarray([row["normalized_abs_phi"] for row in rows])
        rational = 1/(1+rho**2)
        gaussian = np.exp(-np.log(2)*rho**2)
        mask = rho <= 2.0
        family = 1/(1+rho[mask][None, :]**exponents[:, None])
        errors = np.mean((family-profile[mask][None, :])**2, axis=1)
        best = int(np.argmin(errors))
        derivative = np.gradient(profile, rho, edge_order=2)
        output.append({
            "case": case, "N": rows[0]["N"], "level": level,
            "field_step": rows[0]["field_step"],
            "rational_rmse_rho_0_2": float(np.sqrt(np.mean((rational[mask]-profile[mask])**2))),
            "gaussian_rmse_rho_0_2": float(np.sqrt(np.mean((gaussian[mask]-profile[mask])**2))),
            "best_rational_power_m": float(exponents[best]),
            "best_rational_family_rmse": float(np.sqrt(errors[best])),
            "measured_shape_gradient_at_r50": float(abs(np.interp(1.0, rho, derivative))),
        })
    return output


def plot_profile_shapes(profile_rows, fit_rows):
    grouped = {}
    for row in profile_rows:
        grouped.setdefault((row["case"], row["level"]), []).append(row)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    case4000 = next((case for case, level in grouped if "N4000_plus_0.01pct" in case), None)
    if case4000:
        for level in (20, 80, 320, 1280):
            rows = sorted(grouped.get((case4000, level), []),
                          key=lambda row: row["rho_r_over_r50"])
            if rows:
                axes[0].plot([row["rho_r_over_r50"] for row in rows],
                             [row["normalized_abs_phi"] for row in rows],
                             label=fr"$|\phi_0|={level}$")
    rho = PROFILE_RHO
    axes[0].plot(rho, 1/(1+rho**2), "k--", label=r"local-ODE $1/(1+\rho^2)$")
    axes[0].plot(rho, np.exp(-np.log(2)*rho**2), "k:", label="Gaussian")
    axes[0].set(xlabel=r"$\rho=r/R_{50}$", ylabel=r"$|\phi(r)|/|\phi(0)|$",
                title="Normalized core profile (N=4000)", xlim=(0, 2), ylim=(0, 1.03))
    axes[0].grid(True, alpha=.25)
    axes[0].legend(fontsize=8)

    fit_groups = {}
    for row in fit_rows:
        fit_groups.setdefault((row["case"], row["N"]), []).append(row)
    for (case, n), rows in fit_groups.items():
        rows.sort(key=lambda row: row["level"])
        label = f"N={n}, step={rows[0]['field_step']:g} rational"
        axes[1].loglog([row["level"] for row in rows],
                       [row["rational_rmse_rho_0_2"] for row in rows], "o-", label=label)
        axes[1].loglog([row["level"] for row in rows],
                       [row["gaussian_rmse_rho_0_2"] for row in rows], "s--",
                       label=f"N={n}, step={rows[0]['field_step']:g} Gaussian")
    axes[1].set(xlabel=r"central field $|\phi(0)|$", ylabel="profile RMSE over 0 ≤ r/R50 ≤ 2",
                title="Which profile shape fits better?")
    axes[1].grid(True, which="both", alpha=.25)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS/"profile_shape_comparison.png", dpi=180)
    plt.close(fig)


def fit_powerlaw(rows, lower, upper, y_key):
    selected = [row for row in rows if lower <= row["level"] <= upper
                and row[y_key] > 0 and np.isfinite(row[y_key])]
    if len(selected) < 5:
        return None
    times = np.asarray([row["time"] for row in selected])
    yvals = np.log(np.asarray([row[y_key] for row in selected]))
    # Search tc on a logarithmic offset grid, then refine around its best value.
    span = max(float(np.ptp(times)), 1e-5)
    left_delta, right_delta = max(1e-12, span*1e-7), max(0.1, span*20)
    best = None
    for _ in range(4):
        deltas = np.geomspace(left_delta, right_delta, 3000)
        trials = times[-1] + deltas
        log_tau = np.log(trials[:, None]-times[None, :])
        centered = log_tau-log_tau.mean(axis=1, keepdims=True)
        centered_y = yvals-yvals.mean()
        den = np.sum(centered**2, axis=1)
        slopes = np.sum(centered*centered_y, axis=1)/den
        intercepts = yvals.mean()-slopes*log_tau.mean(axis=1)
        residual = yvals[None, :]-(slopes[:, None]*log_tau+intercepts[:, None])
        rss = np.sum(residual**2, axis=1)
        rss[slopes >= 0] = np.inf
        index = int(np.argmin(rss))
        best = (trials[index], slopes[index], rss[index], index, deltas)
        lo = max(1e-12, deltas[max(0, index-3)])
        hi = deltas[min(len(deltas)-1, index+3)]
        left_delta, right_delta = lo, hi
    tc, slope, rss, _, _ = best
    total = float(np.sum((yvals-yvals.mean())**2))
    return float(tc), float(-slope), float(1-rss/total if total else np.nan), len(selected)


def fit_rows(events, case):
    output = []
    for lower, upper in FIT_BANDS:
        alpha_fit = fit_powerlaw(events, lower, upper, "abs_phi_center")
        if not alpha_fit:
            continue
        tc, alpha, alpha_r2, points = alpha_fit
        selected = [row for row in events if lower <= row["level"] <= upper]
        tau = tc-np.asarray([row["time"] for row in selected])
        if np.any(tau <= 0):
            continue
        log_tau = np.log(tau)
        log_r = np.log([row["r50"] for row in selected])
        beta_coeff = np.polyfit(log_tau, log_r, 1)
        beta_slope = float(beta_coeff[0])
        q_coeff = np.polyfit(np.log([row["abs_phi_center"] for row in selected]), log_r, 1)
        q_field = float(-q_coeff[0])
        beta_rss = float(np.sum((log_r-np.polyval(beta_coeff, log_tau))**2))
        r_beta = float(1-beta_rss/np.sum((log_r-log_r.mean())**2))
        log_g = np.log([row["gradient_r50"] for row in selected])
        grad_coeff = np.polyfit(log_tau, log_g, 1)
        grad_slope = float(grad_coeff[0])
        grad_rss = float(np.sum((log_g-np.polyval(grad_coeff, log_tau))**2))
        r_grad = float(1-grad_rss/np.sum((log_g-log_g.mean())**2))
        output.append({
            "case": case["name"], "N": case["n"], "offset_percent": case["offset"],
            "cfl": case["cfl"], "field_step": case["field_step"],
            "level_min": lower, "level_max": upper, "points": points,
            "tc_fit": tc, "alpha_center": alpha, "alpha_r2": alpha_r2,
            "beta_r50": beta_slope, "beta_r2": r_beta,
            "q_r50_vs_center": q_field,
            "alpha_plus_beta": alpha+beta_slope,
            # Report gamma in |d_r phi| ~ tau^(-gamma), hence negate the
            # ordinary log(y)-versus-log(tau) regression slope.
            "gradient_tau_exponent": -grad_slope,
            "gradient_r2": r_grad,
        })
    return output


def term_ratio_fits(events):
    grouped = {}
    for row in events:
        grouped.setdefault(row["case"], []).append(row)
    metrics = ("lap_norm_over_cubic_norm", "mass_norm_over_cubic_norm",
               "rhs_norm_over_cubic_norm")
    output = []
    for case, rows in grouped.items():
        selected = [row for row in rows if row["level"] >= 20]
        x = np.log([row["abs_phi_center"] for row in selected])
        for key in metrics:
            y = np.log([row[key] for row in selected])
            coeff = np.polyfit(x, y, 1)
            residual = y-np.polyval(coeff, x)
            total = np.sum((y-y.mean())**2)
            output.append({
                "case": case, "metric": key, "field_min": selected[0]["level"],
                "field_max": selected[-1]["level"], "points": len(selected),
                "log_ratio_slope_vs_log_center": float(coeff[0]),
                "r2": float(1-np.sum(residual**2)/total if total else np.nan),
            })
    return output


def central_ode_fits(events, scaling_fits):
    """Linearized check of phi(0) ~ C/(tc-t): fit 1/|phi0| against t."""
    grouped = {}
    for row in events:
        grouped.setdefault(row["case"], []).append(row)
    output = []
    for case, rows in grouped.items():
        selected = sorted((row for row in rows if 160 <= row["level"] <= 1280),
                          key=lambda row: row["time"])
        if len(selected) < 5:
            continue
        time = np.asarray([row["time"] for row in selected])
        inverse_center = 1/np.asarray([row["abs_phi_center"] for row in selected])
        slope, intercept = np.polyfit(time, inverse_center, 1)
        residual = inverse_center-(slope*time+intercept)
        total = np.sum((inverse_center-inverse_center.mean())**2)
        if slope >= 0:
            continue
        tc_linear = float(-intercept/slope)
        coefficient = float(-1/slope)
        exponent_fit = next((row for row in scaling_fits
                             if row["case"] == case and row["level_min"] == 160
                             and row["level_max"] == 1280), None)
        tc_power = exponent_fit["tc_fit"] if exponent_fit else float("nan")
        output.append({
            "case": case, "N": selected[0]["N"],
            "offset_percent": selected[0]["offset_percent"],
            "level_min": selected[0]["level"], "level_max": selected[-1]["level"],
            "points": len(selected), "linear_tc": tc_linear,
            "linear_ode_prefactor_C": coefficient,
            "relative_difference_C_from_sqrt2": coefficient/np.sqrt(2)-1,
            "linear_fit_r2": float(1-np.sum(residual**2)/total if total else np.nan),
            "mean_pi_over_phi2": float(np.mean([row["pi_over_phi2_center"] for row in selected])),
            "mean_abs_pi_over_phi2": float(np.mean(
                [abs(row["pi_over_phi2_center"]) for row in selected])),
            "relative_difference_abs_pi_ratio_from_inv_sqrt2": float(
                np.mean([abs(row["pi_over_phi2_center"]) for row in selected])*np.sqrt(2)-1),
            "power_fit_tc": tc_power,
            "tc_difference_linear_minus_power": tc_linear-tc_power,
        })
    return output


def plot_results(events, fits):
    RESULTS.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    grouped = {}
    for row in events:
        grouped.setdefault(row["case"], []).append(row)
    baseline = [name for name in grouped if "cfl0.4_step0.15" in name
                and ("plus_0.01pct" in name or "plus_0.05pct" in name
                     or "plus_0.1pct" in name)]
    resolution_cases = [name for name in baseline if "plus_0.01pct" in name]
    for name in resolution_cases:
        rows = grouped[name]
        rows.sort(key=lambda row: row["level"])
        label = f"N={rows[0]['N']}"
        center = np.asarray([row["abs_phi_center"] for row in rows])
        axes[1].loglog(center, [row["r50"] for row in rows], "o-", label=label)
    representative = next((rows for name, rows in grouped.items()
                           if "N2000_plus_0.01pct_cfl0.4_step0.15" in name), None)
    if representative:
        representative.sort(key=lambda row: row["level"])
        center = np.asarray([row["abs_phi_center"] for row in representative])
        for key, label, marker in (
            ("lap_norm_over_cubic_norm", r"$\|\Delta\phi\|/\|\phi^3\|$", "o-"),
            ("mass_norm_over_cubic_norm", r"$\|\phi\|/\|\phi^3\|$", "s-"),
            ("rhs_norm_over_cubic_norm", r"$\|\mathrm{RHS}\|/\|\phi^3\|$", "^-"),
        ):
            axes[0].loglog(center, [row[key] for row in representative], marker, label=label)
        xref = np.asarray([20.0, 1280.0])
        axes[0].loglog(xref, 0.12*(xref/20.0)**-1, "k--", alpha=.55,
                       label=r"reference $\propto |\phi_0|^{-1}$")
        axes[0].loglog(xref, 4.8e-3*(xref/20.0)**-2, "k:", alpha=.65,
                       label=r"reference $\propto |\phi_0|^{-2}$")
    axes[0].set(xlabel=r"central field $|\phi(0)|$",
                ylabel="volume-weighted norm ratio inside R50",
                title="Core PDE balance (N=2000, +0.01%)")
    axes[1].set(xlabel=r"central field $|\phi(0)|$", ylabel=r"half-height radius $R_{50}$",
                title="Core width versus amplitude")
    for ax in axes:
        ax.grid(True, which="both", alpha=.25)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(RESULTS/"pde_balance_and_width.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    for case in baseline:
        rows = grouped[case]
        fit = next((row for row in fits if row["case"] == case
                    and row["level_min"] == 160 and row["level_max"] == 1280), None)
        if not fit:
            continue
        rows.sort(key=lambda row: row["level"])
        rows = [row for row in rows if 160 <= row["level"] <= 1280]
        tau = fit["tc_fit"]-np.asarray([row["time"] for row in rows])
        good = tau > 0
        gradient = np.asarray([row["gradient_r50"] for row in rows])
        ax.loglog(tau[good], gradient[good], "o-",
                  label=f"N={grouped[case][0]['N']}, +{grouped[case][0]['offset_percent']:g}% "
                        f"(gamma={fit['gradient_tau_exponent']:.3f})")
    ax.set(xlabel=r"remaining time $\tau=t_c-t$", ylabel=r"$|\partial_r\phi|_{R_{50}}$",
           title="Direct test of gradient scaling")
    ax.grid(True, which="both", alpha=.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS/"gradient_scaling.png", dpi=180)
    plt.close(fig)


def parse_case(raw):
    parts = raw.split(":")
    if len(parts) not in (2, 3, 4):
        raise argparse.ArgumentTypeError("case format is N:offset[:CFL[:field_step]]")
    try:
        n, offset = int(parts[0]), float(parts[1])
        cfl = float(parts[2]) if len(parts) >= 3 else 0.4
        field_step = float(parts[3]) if len(parts) >= 4 else 0.15
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error
    if n < 100 or offset <= 0 or cfl <= 0 or field_step <= 0:
        raise argparse.ArgumentTypeError("N, offset, CFL and field_step must be positive (N>=100)")
    amplitude = A_REFERENCE*(1+offset/100)
    return {"n": n, "offset": offset, "cfl": cfl, "field_step": field_step,
            "amplitude": amplitude,
            "name": f"N{n}_plus_{offset:g}pct_cfl{cfl:g}_step{field_step:g}"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=parse_case,
                        help="N:offset_percent[:CFL[:field_step]], repeatable")
    parser.add_argument("--max-level", type=float, default=1280.0)
    parser.add_argument("--profiles-only", action="store_true",
                        help="capture normalized radial profiles and compare ODE/Gaussian shapes")
    args = parser.parse_args()
    if args.max_level < 80:
        parser.error("max-level must be at least 80")
    if args.case:
        cases = args.case
    elif args.profiles_only:
        cases = [parse_case(value) for value in (
            "1000:0.01", "1000:0.01:0.4:0.05", "2000:0.01", "4000:0.01")]
    else:
        cases = [parse_case(value) for value in (
            "500:0.01", "1000:0.01", "2000:0.01", "4000:0.01",
            "1000:0.01:0.2", "1000:0.01:0.4:0.05", "1000:0.05", "1000:0.1")]
    solver = load_solver()
    if args.profiles_only:
        profiles = []
        for case in cases:
            _, _, _, case_profiles = simulate(
                solver, case, args.max_level, capture_profile_levels=PROFILE_LEVELS,
            )
            profiles.extend(case_profiles)
        shape_fits = fit_profile_shapes(profiles)
        write_csv(RESULTS/"profile_shape_snapshots.csv", profiles)
        write_csv(RESULTS/"profile_shape_fits.csv", shape_fits)
        plot_profile_shapes(profiles, shape_fits)
        print(f"Wrote profile-shape audit under {RESULTS}", flush=True)
        return
    all_events, all_fits, settings = [], [], []
    for case in cases:
        events, steps, final_time, _ = simulate(solver, case, args.max_level)
        all_events.extend(events)
        all_fits.extend(fit_rows(events, case))
        settings.append({**case, "steps": steps, "final_time": final_time,
                         "highest_level_reached": max((row["level"] for row in events), default="")})
    write_csv(RESULTS/"event_diagnostics.csv", all_events)
    write_csv(RESULTS/"scaling_fits.csv", all_fits)
    write_csv(RESULTS/"term_ratio_fits.csv", term_ratio_fits(all_events))
    write_csv(RESULTS/"central_ode_fits.csv", central_ode_fits(all_events, all_fits))
    write_csv(RESULTS/"run_settings.csv", settings)
    plot_results(all_events, all_fits)
    print(f"Wrote diagnostics under {RESULTS}", flush=True)


if __name__ == "__main__":
    main()
