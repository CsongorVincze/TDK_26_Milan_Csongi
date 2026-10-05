#!/usr/bin/env python3
"""Standalone Module 4 collapse study using the field half-radius R50."""

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
SAMPLE_INTERVAL = 0.01
CRITICAL_LEVEL = 12.0
BASE_LEVEL = 20.0
OFFSETS = (1e-4, 5e-4, 1e-3)
FIT_WINDOWS = ((0.10, 0.30), (0.20, 0.60), (0.40, 1.00), (0.80, 1.60))


@dataclass
class Config:
    n: int = 500
    r_max: float = 25.0
    stretch: float = 5.0
    sponge_fraction: float = 0.8
    sponge_strength: float = 7.0
    sponge_power: float = 3.0
    cfl: float = 0.4
    r0: float = 2.0
    sample_interval: float = SAMPLE_INTERVAL


@dataclass
class Grid:
    r: np.ndarray
    rx: np.ndarray
    rxx: np.ndarray
    damping: np.ndarray
    dx: float
    dr_min: float
    dt: float


def make_grid(c):
    x = np.linspace(0.0, 1.0, c.n)
    dx = x[1] - x[0]
    r = c.r_max * np.sinh(c.stretch*x) / np.sinh(c.stretch)
    rx = c.r_max*c.stretch*np.cosh(c.stretch*x) / np.sinh(c.stretch)
    rxx = c.stretch**2*r
    sponge_start = c.sponge_fraction*c.r_max
    damping = c.sponge_strength*np.maximum(
        (r-sponge_start)/(c.r_max-sponge_start), 0.0
    )**c.sponge_power
    dr_min = rx[0]*dx
    return Grid(r, rx, rxx, damping, dx, dr_min, c.cfl*dr_min)


def laplacian(phi, grid):
    """Spherically symmetric 3D Laplacian on the stretched radial grid."""
    result = np.zeros_like(phi)
    result[0] = 6.0*(phi[1]-phi[0])/(grid.rx[0]*grid.dx)**2
    first = (phi[2:]-phi[:-2])/(2.0*grid.dx)
    second = (phi[2:]-2.0*phi[1:-1]+phi[:-2])/grid.dx**2
    result[1:-1] = second/grid.rx[1:-1]**2 + (
        2.0/(grid.r[1:-1]*grid.rx[1:-1])
        - grid.rxx[1:-1]/grid.rx[1:-1]**3
    )*first
    return result


def rhs(phi, pi, grid):
    dphi = pi.copy()
    dpi = laplacian(phi, grid) - (phi-phi**3) - grid.damping*pi
    dphi[-1] = dpi[-1] = 0.0
    return dphi, dpi


def rk4(phi, pi, dt, grid):
    k1 = rhs(phi, pi, grid)
    k2 = rhs(phi+dt*k1[0]/2, pi+dt*k1[1]/2, grid)
    k3 = rhs(phi+dt*k2[0]/2, pi+dt*k2[1]/2, grid)
    k4 = rhs(phi+dt*k3[0], pi+dt*k3[1], grid)
    new_phi = phi+dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
    new_pi = pi+dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
    new_phi[-1] = new_pi[-1] = 0.0
    return new_phi, new_pi


def half_radius(phi, radius):
    """Interpolate the first radius with |phi(r)| <= |phi(0)|/2."""
    values = np.abs(phi)
    if values[0] <= 0.0:
        return float("nan")
    crossings = np.flatnonzero(values[1:] <= 0.5*values[0])
    if not len(crossings):
        return float("nan")
    i = int(crossings[0]+1)
    y0, y1 = values[i-1], values[i]
    if y0 == y1:
        return float(radius[i])
    frac = np.clip((0.5*values[0]-y0)/(y1-y0), 0.0, 1.0)
    return float(radius[i-1]+frac*(radius[i]-radius[i-1]))


def sample(t, phi, pi, grid):
    r50 = half_radius(phi, grid.r)
    return {
        "time": float(t),
        "phi_center": float(phi[0]),
        "abs_phi_center": float(abs(phi[0])),
        "pi_center": float(pi[0]),
        "r50": r50,
        "r50_squared": r50*r50 if np.isfinite(r50) else float("nan"),
    }


def evolve(amplitude, n, t_max, sample_interval, stop_level=None, record=True):
    c = Config(n=n, sample_interval=sample_interval)
    grid = make_grid(c)
    phi = amplitude*np.exp(-grid.r**2/(2*c.r0**2))
    pi = np.zeros_like(phi)
    phi[-1] = pi[-1] = 0.0
    t, steps, next_sample = 0.0, 0, sample_interval
    rows = [sample(t, phi, pi, grid)] if record else []
    crossed = stop_level is not None and abs(phi[0]) >= stop_level
    failed = False

    while t < t_max and not crossed:
        peak = max(float(np.max(np.abs(phi))), 1.0)
        dt = min(grid.dt, 0.15/peak, t_max-t)
        if record:
            dt = min(dt, max(next_sample-t, 1e-12))
        next_phi, next_pi = rk4(phi, pi, dt, grid)
        if not np.all(np.isfinite(next_phi)) or not np.all(np.isfinite(next_pi)):
            failed = True
            break
        crossed = stop_level is not None and abs(next_phi[0]) >= stop_level
        if crossed:
            change = abs(next_phi[0])-abs(phi[0])
            frac = np.clip((stop_level-abs(phi[0]))/change, 0.0, 1.0)
            t += dt*frac
            phi += frac*(next_phi-phi)
            pi += frac*(next_pi-pi)
            phi[0] = np.copysign(stop_level, next_phi[0])
        else:
            t += dt
            phi, pi = next_phi, next_pi
        steps += 1
        if record and (t >= next_sample-1e-11 or crossed):
            rows.append(sample(t, phi, pi, grid))
            while next_sample <= t+1e-11:
                next_sample += sample_interval

    if record and rows[-1]["time"] < t-1e-11:
        rows.append(sample(t, phi, pi, grid))
    return rows, crossed, failed, (t if crossed else float("nan")), steps


def search_critical(args):
    history = []

    def classify(amplitude):
        _, crossed, failed, event_time, steps = evolve(
            amplitude, args.n, args.critical_tmax, args.sample_interval,
            CRITICAL_LEVEL, record=False,
        )
        if failed:
            raise RuntimeError(f"Numerical failure during threshold search at A={amplitude:.10f}.")
        history.append({"amplitude": amplitude, "collapsed": int(crossed),
                        "event_time": event_time if crossed else "", "steps": steps})
        print(f"threshold trial A={amplitude:.10f}: "
              f"{'collapse' if crossed else 'no cutoff by horizon'}", flush=True)
        return crossed

    lo, hi = args.a_min, args.a_max
    while classify(lo):
        lo *= 0.5
    while not classify(hi):
        hi *= 1.5
        if hi > 10:
            raise RuntimeError("Could not find a collapsing upper amplitude.")
    while hi-lo > args.critical_tolerance:
        mid = 0.5*(lo+hi)
        if classify(mid):
            hi = mid
        else:
            lo = mid
    return lo, hi, history


def write_rows(path, rows):
    if not rows:
        return
    columns = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def fit_line(x, y):
    if len(x) < 5 or np.ptp(x) == 0:
        return float("nan"), float("nan"), float("nan")
    slope, intercept = np.polyfit(x, y, 1)
    rss = np.sum((y-(slope*x+intercept))**2)
    total = np.sum((y-y.mean())**2)
    return float(slope), float(intercept), float(1-rss/total if total else float("nan"))


def estimate_tc_joint(rows, count=40):
    """Jointly fit tc and alpha on the late, high-field central-field tail."""
    time = np.asarray([row["time"] for row in rows])
    center = np.asarray([row["abs_phi_center"] for row in rows])
    ids = np.flatnonzero(center >= 4.0)[-count:]
    if len(ids) < 8:
        return None
    fit_time, log_center = time[ids], np.log(center[ids])
    left, right = time[-1]+1e-5, time[-1]+1.0
    best = None
    for _ in range(4):
        trial_tc = np.linspace(left, right, 1201)
        log_tau = np.log(trial_tc[:, None]-fit_time[None, :])
        centered_x = log_tau-log_tau.mean(axis=1, keepdims=True)
        centered_y = log_center-log_center.mean()
        denom = np.sum(centered_x**2, axis=1)
        slopes = np.sum(centered_x*centered_y, axis=1)/denom
        intercepts = log_center.mean()-slopes*log_tau.mean(axis=1)
        residual = log_center[None, :]-(slopes[:, None]*log_tau+intercepts[:, None])
        rss = np.sum(residual**2, axis=1)
        rss[slopes >= 0] = np.inf
        index = int(np.argmin(rss))
        step = (right-left)/(len(trial_tc)-1)
        best = (float(trial_tc[index]), float(slopes[index]), float(rss[index]), step)
        left = max(time[-1]+1e-6, best[0]-5*step)
        right = best[0]+5*step
    tc, slope, rss, _ = best
    total = np.sum((log_center-log_center.mean())**2)
    r2 = 1-rss/total if total > 0 else float("nan")
    return tc, -slope, r2, len(ids)


def fit_results(cases):
    fits, scaling = [], []
    for case in cases:
        rows = case["rows"]
        time = np.asarray([row["time"] for row in rows])
        center = np.asarray([row["abs_phi_center"] for row in rows])
        radius = np.asarray([row["r50"] for row in rows])
        joint = estimate_tc_joint(rows)
        tc_methods = []
        if joint:
            tc_methods.append(("joint_tc_alpha", joint[0], joint[1], joint[2], joint[3]))
        tc_estimate = joint[0] if joint else float("nan")
        for method, tc, alpha, tc_r2, tc_n in tc_methods:
            tau = tc-time
            for lower, upper in FIT_WINDOWS:
                use = ((tau >= lower) & (tau <= upper) & (center >= 4.0)
                       & np.isfinite(radius)
                       & (radius > 0))
                x, y = np.log10(tau[use]), np.log10(radius[use])
                beta, intercept, r2 = fit_line(x, y)
                fits.append({
                    "case": case["case"], "N": case.get("n", 500),
                    "amplitude": case["amplitude"],
                    "offset_percent": 100*case["offset"],
                    "tc_method": method, "tc": tc,
                    "alpha_central": alpha,
                    "tc_fit_r2": tc_r2, "tc_fit_points": tc_n,
                    "tau_min": lower, "tau_max": upper, "n": int(use.sum()),
                    "beta_r50": beta, "fit_r2": r2,
                })
        ids = np.flatnonzero(center >= 4.0)[-40:]
        use = ids[np.isfinite(radius[ids]) & (radius[ids] > 0)]
        q, intercept, r2 = fit_line(np.log10(center[use]), np.log10(radius[use])) \
            if len(use) >= 5 else (float("nan"), float("nan"), float("nan"))
        scaling.append({
            "case": case["case"], "N": case.get("n", 500),
            "amplitude": case["amplitude"],
            "offset_percent": 100*case["offset"], "tc_estimate": tc_estimate,
            "central_field_min": float(center[use[0]]) if len(use) else "",
            "central_field_max": float(center[use[-1]]) if len(use) else "",
            "n": len(use), "beta_over_alpha_q": -q, "fit_r2": r2,
        })
    return fits, scaling


def plot_critical_search(history, acrit, out, tmax):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for collapsed, marker, label, color in (
        (0, "x", "no cutoff by horizon", "#2563eb"),
        (1, "o", "reached |phi(0)|=12", "#dc2626"),
    ):
        rows = [row for row in history if row["collapsed"] == collapsed]
        if rows:
            y = [row["event_time"] if collapsed else tmax for row in rows]
            ax.scatter([row["amplitude"] for row in rows], y,
                       marker=marker, color=color, label=label)
    ax.axvline(acrit, color="black", ls="--", label=f"bracket midpoint {acrit:.9f}")
    ax.set(xlabel="initial amplitude A", ylabel="event time (or horizon)",
           title="Finite-horizon critical-amplitude search")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_central_fields(cases, out, level):
    fig, ax = plt.subplots(figsize=(9, 5))
    for case in cases:
        rows = case["rows"]
        ax.plot([r["time"] for r in rows], [r["phi_center"] for r in rows],
                label=case["case"])
    ax.axhline(CRITICAL_LEVEL, color="firebrick", ls=":", lw=1,
               label=rf"critical-search cutoff $|\phi_0|={CRITICAL_LEVEL:g}$")
    ax.axhline(-CRITICAL_LEVEL, color="firebrick", ls=":", lw=1)
    ax.axhline(level, color="black", ls="--", lw=0.9,
               label=rf"supercritical cutoff $|\phi_0|={level:g}$")
    ax.axhline(-level, color="black", ls="--", lw=0.9)
    ax.set(xlabel="time t", ylabel=r"central field $\phi(0,t)$",
           title="Critical-bracket and supercritical central-field histories")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_r50_time(cases, out):
    fig, ax = plt.subplots(figsize=(9, 5))
    for case in cases:
        rows = case["rows"]
        center = np.asarray([r["abs_phi_center"] for r in rows])
        radius = np.asarray([r["r50"] for r in rows])
        time = np.asarray([r["time"] for r in rows])
        use = (center >= 4.0) & np.isfinite(radius) & (radius > 0)
        ax.plot(time[use], radius[use],
                label=f"{case['case']} ({case['amplitude']:.9f})")
    ax.set(xlabel="time t", ylabel=r"half-central-field radius $R_{50}(t)$",
           title=r"Half-radius during collapse ($|\phi(0,t)|\geq4$)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_loglog(cases, fits, out, tc_method="joint_tc_alpha"):
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.get_cmap("tab10").colors
    for color, case in zip(colors, cases):
        tc_fit = estimate_tc_joint(case["rows"])
        if not tc_fit:
            continue
        tc = tc_fit[0]
        rows = case["rows"]
        tau = tc-np.asarray([r["time"] for r in rows])
        center = np.asarray([r["abs_phi_center"] for r in rows])
        radius = np.asarray([r["r50"] for r in rows])
        good = ((tau >= 0.1) & (tau <= 1.6) & (center >= 4.0)
                & np.isfinite(radius) & (radius > 0))
        ax.plot(np.log10(tau[good]), np.log10(radius[good]), ".", ms=2.7,
                color=color, alpha=0.55, label=f"{case['case']} data")
        fit = next((row for row in fits if row["case"] == case["case"]
                    and row["tc_method"] == tc_method
                    and row["tau_min"] == 0.2 and row["tau_max"] == 0.6), None)
        if fit and np.isfinite(fit["beta_r50"]):
            x = np.linspace(np.log10(0.2), np.log10(0.6), 60)
            ids = ((tau >= 0.2) & (tau <= 0.6) & (center >= 4.0)
                   & np.isfinite(radius) & (radius > 0))
            beta, intercept, _ = fit_line(np.log10(tau[ids]), np.log10(radius[ids]))
            ax.plot(np.log10(tau[ids]), np.log10(radius[ids]), "o", ms=3.2,
                    mfc="none", mec=color, label=f"{case['case']} fit points")
            ax.plot(x, beta*x+intercept, "--", color=color,
                    label=rf"$\beta={fit['beta_r50']:.3f}$, $R^2={fit['fit_r2']:.3f}$")
    method_title = r"joint high-field $t_c,\alpha$ fit"
    ax.set(xlabel=r"$\log_{10}(t_c-t)$", ylabel=r"$\log_{10}R_{50}$",
           title=rf"Half-radius collapse fits ({method_title}): $R_{{50}}\sim(t_c-t)^\beta$")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(out, dpi=170)
    plt.close(fig)


def plot_width_field(cases, scaling, out):
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.get_cmap("tab10").colors
    for color, case in zip(colors, cases):
        rows = case["rows"]
        center = np.asarray([r["abs_phi_center"] for r in rows])
        radius = np.asarray([r["r50"] for r in rows])
        use = (center >= 4) & np.isfinite(radius) & (radius > 0)
        summary = next(row for row in scaling if row["case"] == case["case"])
        ax.plot(np.log10(center[use]), np.log10(radius[use]), ".", ms=3,
                color=color, label=rf"{case['case']}: $q={summary['beta_over_alpha_q']:.3f}$")
    ax.set(xlabel=r"$\log_{10}|\phi(0,t)|$", ylabel=r"$\log_{10}R_{50}$",
           title=r"Direct width scaling, $R_{50}\sim|\phi(0,t)|^{-q}$")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=170)
    plt.close(fig)


def plot_continuation(continuation, out):
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.get_cmap("viridis")(np.linspace(0.12, 0.88, len(continuation)))
    for color, case in zip(colors, continuation):
        rows = case["rows"]
        center = np.asarray([r["abs_phi_center"] for r in rows])
        radius = np.asarray([r["r50"] for r in rows])
        use = (center >= 4.0) & np.isfinite(radius) & (radius > 0)
        ax.plot(np.log10(center[use]), np.log10(radius[use]), ".-", ms=2.5,
                lw=1, color=color,
                label=f"N={case['n']}, stop at |phi0|={case['level']:g}")
    ax.set(xlabel=r"$\log_{10}|\phi(0,t)|$", ylabel=r"$\log_{10}R_{50}$",
           title="Half-radius scaling continued toward stronger collapse")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=170)
    plt.close(fig)


def plot_resolution(cases, out):
    fig, ax = plt.subplots(figsize=(8, 6))
    for case in cases:
        rows = case["rows"]
        center = np.asarray([r["abs_phi_center"] for r in rows])
        radius = np.asarray([r["r50"] for r in rows])
        use = (center >= 4.0) & np.isfinite(radius) & (radius > 0)
        ax.plot(np.log10(center[use]), np.log10(radius[use]), ".-", ms=2.5,
                lw=1, label=f"N={case.get('n', 500)}")
    ax.set(xlabel=r"$\log_{10}|\phi(0,t)|$", ylabel=r"$\log_{10}R_{50}$",
           title="Resolution check at fixed supercritical amplitude")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=170)
    plt.close(fig)


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        parsed = []
        for row in csv.DictReader(handle):
            values = {}
            for key, value in row.items():
                try:
                    values[key] = float(value)
                except ValueError:
                    values[key] = value
            parsed.append(values)
        return parsed


def reanalyze_saved_runs():
    """Recreate fits and plots from this folder's existing trajectory CSVs."""
    summary = read_rows(OUT/"supercritical_sweep_summary.csv")
    cases = [{
        "case": row["case"], "amplitude": row["amplitude"],
        "offset": row["relative_offset_percent"]/100.0,
        "rows": read_rows(OUT/row["csv_file"]),
    } for row in summary]
    critical = [{
        "case": name, "amplitude": rows[0]["abs_phi_center"], "rows": rows,
    } for name, rows in (
        ("critical_lower_boundary", read_rows(OUT/"critical_lower_boundary.csv")),
        ("critical_upper_boundary", read_rows(OUT/"critical_upper_boundary.csv")),
    )]
    trials = read_rows(OUT/"critical_threshold_search.csv")
    lo = max(row["amplitude"] for row in trials if row["collapsed"] == 0.0)
    hi = min(row["amplitude"] for row in trials if row["collapsed"] == 1.0)
    settings_file = OUT/"run_settings.csv"
    if settings_file.exists():
        settings = read_rows(settings_file)[0]
    else:
        settings = {
            "grid_N": 500, "a_min_start": 0.1, "a_max_start": 1.2,
            "critical_horizon": 100.0, "critical_tolerance": 1e-8,
            "supercritical_horizon": 100.0,
            "sample_interval": SAMPLE_INTERVAL,
            "critical_cutoff": CRITICAL_LEVEL,
            "supercritical_cutoff": BASE_LEVEL,
        }
        write_rows(settings_file, [settings])
    write_rows(OUT/"critical_summary.csv", [{
        "grid_N": settings["grid_N"], "horizon": settings["critical_horizon"],
        "collapse_cutoff_abs_phi0": CRITICAL_LEVEL,
        "amplitude_lower_noncollapsing": lo,
        "amplitude_upper_collapsing": hi,
        "Acrit_midpoint_estimate": 0.5*(lo+hi),
        "bracket_width": hi-lo, "tolerance": settings["critical_tolerance"],
    }])
    continuation_specs = ((500, 30), (500, 40), (500, 60), (500, 80), (1000, 80))
    continuation = [{
        "case": f"continuation_N{n}_phi{level}", "n": n, "level": level,
        "amplitude": cases[0]["amplitude"], "offset": cases[0]["offset"],
        "rows": read_rows(OUT/f"continuation_N{n}_phi{level}.csv"),
    } for n, level in continuation_specs]
    resolution = [{
        "case": "supercritical_plus_0.01pct", "n": 500,
        "amplitude": cases[0]["amplitude"], "offset": cases[0]["offset"],
        "rows": cases[0]["rows"],
    }]
    for n in (250, 1000):
        resolution.append({
            "case": f"resolution_N{n}_phi20", "n": n,
            "amplitude": cases[0]["amplitude"], "offset": cases[0]["offset"],
            "rows": read_rows(OUT/f"resolution_N{n}_phi20.csv"),
        })

    fits, scaling = fit_results(cases+continuation+resolution[1:])
    write_rows(OUT/"beta_fit_sensitivity.csv", fits)
    write_rows(OUT/"width_vs_center_scaling.csv", scaling)
    plot_central_fields(critical+cases, OUT/"critical_and_supercritical_phi.png", BASE_LEVEL)
    plot_r50_time(cases, OUT/"supercritical_r50_vs_time.png")
    plot_continuation(continuation, OUT/"half_radius_continuation.png")
    plot_resolution(resolution, OUT/"half_radius_resolution.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=500, help="grid size for Acrit search")
    parser.add_argument("--a-min", type=float, default=0.1)
    parser.add_argument("--a-max", type=float, default=1.2)
    parser.add_argument("--critical-tmax", type=float, default=100.0)
    parser.add_argument("--critical-tolerance", type=float, default=1e-8)
    parser.add_argument("--t-max", type=float, default=100.0)
    parser.add_argument("--sample-interval", type=float, default=SAMPLE_INTERVAL)
    parser.add_argument("--analyze-existing", action="store_true",
                        help="redo fits/plots from saved CSVs without rerunning PDEs")
    args = parser.parse_args()
    if args.n < 20 or args.t_max <= 0 or args.critical_tmax <= 0:
        parser.error("Grid size must be >=20 and both time horizons must be positive.")
    if args.sample_interval <= 0 or args.critical_tolerance <= 0:
        parser.error("Sample interval and critical tolerance must be positive.")
    OUT.mkdir(parents=True, exist_ok=True)
    if args.analyze_existing:
        try:
            reanalyze_saved_runs()
        except FileNotFoundError as error:
            parser.error(f"Cannot analyze existing output; missing {error.filename or error}")
        print(f"Re-analyzed saved CSVs under {OUT}")
        return

    print("Searching the finite-horizon critical bracket...", flush=True)
    lo, hi, history = search_critical(args)
    acrit = 0.5*(lo+hi)
    write_rows(OUT/"critical_threshold_search.csv", history)
    write_rows(OUT/"run_settings.csv", [{
        "grid_N": args.n, "a_min_start": args.a_min, "a_max_start": args.a_max,
        "critical_horizon": args.critical_tmax,
        "critical_tolerance": args.critical_tolerance,
        "supercritical_horizon": args.t_max,
        "sample_interval": args.sample_interval,
        "critical_cutoff": CRITICAL_LEVEL, "supercritical_cutoff": BASE_LEVEL,
    }])
    write_rows(OUT/"critical_summary.csv", [{
        "grid_N": args.n, "horizon": args.critical_tmax,
        "collapse_cutoff_abs_phi0": CRITICAL_LEVEL,
        "amplitude_lower_noncollapsing": lo,
        "amplitude_upper_collapsing": hi,
        "Acrit_midpoint_estimate": acrit,
        "bracket_width": hi-lo,
        "tolerance": args.critical_tolerance,
    }])
    critical_cases = []
    for name, amplitude in (("critical_lower_boundary", lo),
                            ("critical_upper_boundary", hi)):
        rows, _, _, _, _ = evolve(amplitude, args.n, args.critical_tmax,
                                  args.sample_interval, CRITICAL_LEVEL)
        write_rows(OUT/f"{name}.csv", rows)
        critical_cases.append({"case": name, "amplitude": amplitude, "rows": rows})
    plot_critical_search(history, acrit, OUT/"critical_search.png", args.critical_tmax)

    cases = []
    for offset in OFFSETS:
        amplitude = acrit*(1+offset)
        name = f"supercritical_plus_{100*offset:g}pct"
        print(f"Running {name}: A={amplitude:.12f}", flush=True)
        rows, crossed, failed, event_time, steps = evolve(
            amplitude, 500, args.t_max, args.sample_interval, BASE_LEVEL,
        )
        case = {"case": name, "amplitude": amplitude, "offset": offset,
                "n": 500,
                "rows": rows, "crossed": crossed, "failed": failed,
                "event_time": event_time, "steps": steps}
        write_rows(OUT/f"{name}.csv", rows)
        cases.append(case)

    baseline = cases[0]
    central_cases = critical_cases + cases
    plot_central_fields(central_cases, OUT/"critical_and_supercritical_phi.png",
                        BASE_LEVEL)
    plot_r50_time(cases, OUT/"supercritical_r50_vs_time.png")

    continuation = []
    for n, level in ((500, 30.0), (500, 40.0), (500, 60.0),
                     (500, 80.0), (1000, 80.0)):
        label = f"continuation_N{n}_phi{int(level)}"
        print(f"Running {label}: same A={baseline['amplitude']:.12f}", flush=True)
        rows, crossed, failed, event_time, steps = evolve(
            baseline["amplitude"], n, args.t_max, args.sample_interval, level,
        )
        write_rows(OUT/f"{label}.csv", rows)
        continuation.append({"case": label, "n": n, "level": level,
                             "amplitude": baseline["amplitude"],
                             "offset": baseline["offset"],
                             "crossed": crossed, "failed": failed,
                             "event_time": event_time, "steps": steps,
                             "rows": rows})

    resolution = [baseline]
    for n in (250, 1000):
        label = f"resolution_N{n}_phi20"
        print(f"Resolution check {label}", flush=True)
        rows, crossed, failed, event_time, steps = evolve(
            baseline["amplitude"], n, args.t_max, args.sample_interval, BASE_LEVEL,
        )
        write_rows(OUT/f"{label}.csv", rows)
        resolution.append({"case": label, "amplitude": baseline["amplitude"],
                           "offset": baseline["offset"],
                           "n": n, "rows": rows, "crossed": crossed, "failed": failed,
                           "event_time": event_time, "steps": steps})

    plot_continuation(continuation, OUT/"half_radius_continuation.png")
    plot_resolution(resolution, OUT/"half_radius_resolution.png")

    fit_cases = cases+continuation+resolution[1:]
    fits, scaling = fit_results(fit_cases)
    write_rows(OUT/"beta_fit_sensitivity.csv", fits)
    write_rows(OUT/"width_vs_center_scaling.csv", scaling)
    sweep_summary = [{
        "case": case["case"], "amplitude": case["amplitude"],
        "relative_offset_percent": 100*case["offset"],
        "reached_phi20": case["crossed"], "event_time": case["event_time"],
        "steps": case["steps"], "numerical_failure": case["failed"],
        "csv_file": f"{case['case']}.csv",
    } for case in cases]
    write_rows(OUT/"supercritical_sweep_summary.csv", sweep_summary)
    resolution_summary = []
    for case in resolution:
        last = case["rows"][-1]
        resolution_summary.append({
            "case": case["case"], "n": case.get("n", 500),
            "amplitude": case["amplitude"], "reached_phi20": case["crossed"],
            "event_time": case["event_time"], "r50_at_stop": last["r50"],
            "central_field_at_stop": last["abs_phi_center"],
        })
    write_rows(OUT/"resolution_summary.csv", resolution_summary)

    print(f"Acrit bracket [{lo:.10f}, {hi:.10f}], midpoint {acrit:.10f}")
    print(f"Saved CSVs and plots in {OUT}")


if __name__ == "__main__":
    main()
