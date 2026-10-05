#!/usr/bin/env python3
"""Test a time-dependent Gaussian core-width model during Module 4 collapse."""

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
R_MAX, STRETCH, R0 = 25.0, 5.0, 2.0
SPONGE_FRACTION, SPONGE_STRENGTH, SPONGE_POWER = 0.8, 7.0, 3.0
SAMPLE_INTERVAL, A_CRIT_REFERENCE = 0.01, 1.202609700895846
CRITICAL_LEVEL, RUN_LEVEL = 12.0, 20.0
FIT_WINDOWS = ((0.1, 0.3), (0.2, 0.6), (0.4, 1.0), (0.8, 1.6))
OFFSETS = (1e-4, 5e-4, 1e-3)
R50_FACTOR = np.sqrt(2.0*np.log(2.0))
GAUSSIAN_FIT_FRACTIONS = (0.30, 0.90)


@dataclass
class Grid:
    r: np.ndarray
    rx: np.ndarray
    rxx: np.ndarray
    damping: np.ndarray
    dx: float
    dt: float


def make_grid(n, cfl=0.4):
    x = np.linspace(0.0, 1.0, n)
    dx = x[1]-x[0]
    r = R_MAX*np.sinh(STRETCH*x)/np.sinh(STRETCH)
    rx = R_MAX*STRETCH*np.cosh(STRETCH*x)/np.sinh(STRETCH)
    rxx = STRETCH**2*r
    sponge_start = SPONGE_FRACTION*R_MAX
    damping = SPONGE_STRENGTH*np.maximum(
        (r-sponge_start)/(R_MAX-sponge_start), 0.0
    )**SPONGE_POWER
    return Grid(r, rx, rxx, damping, dx, cfl*rx[0]*dx)


def laplacian(phi, grid):
    """Spherically symmetric 3-D Laplacian on the stretched radial grid."""
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
    dpi = laplacian(phi, grid)-(phi-phi**3)-grid.damping*pi
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


def gaussian_width(phi, radius, fit_radius=8.0):
    """Fit log(|phi|/|phi(0)|) = -r^2/(2 sigma^2) over the inner core."""
    center = abs(float(phi[0]))
    if center <= 0.0:
        return (np.nan, np.nan, np.nan, np.nan, 0)
    ratio = np.abs(phi)/center
    lower, upper = GAUSSIAN_FIT_FRACTIONS
    below = np.flatnonzero(ratio[1:] <= lower)
    core_end = int(below[0]+1) if len(below) else len(radius)-1
    use = ((radius > 0.0) & (radius <= fit_radius) & (ratio >= lower)
           & (ratio <= upper) & (np.arange(len(radius)) <= core_end))
    if use.sum() < 5:
        return (np.nan, np.nan, np.nan, np.nan, int(use.sum()))
    x = radius[use]**2
    y = np.log(ratio[use])
    slope = float(np.dot(x, y)/np.dot(x, x))  # normalized Gaussian: intercept is zero
    sigma = np.sqrt(-1.0/(2.0*slope)) if slope < 0.0 else np.nan
    residual = y-slope*x
    total = np.sum((y-y.mean())**2)
    r2 = 1.0-np.sum(residual**2)/total if total > 0 else np.nan
    free_slope, intercept = np.polyfit(x, y, 1)
    return (float(sigma), float(sigma*R50_FACTOR), float(r2),
            float(intercept), int(use.sum()))


def sample(t, phi, grid, fit_radius=8.0):
    sigma, r50, quality, intercept, count = gaussian_width(phi, grid.r, fit_radius)
    return {
        "time": float(t), "phi_center": float(phi[0]),
        "abs_phi_center": float(abs(phi[0])), "sigma_gaussian": sigma,
        "r50_gaussian": r50, "log_shape_r2": quality,
        "free_log_intercept": intercept, "fit_points": count,
    }


def evolve(amplitude, n, t_max, stop_level, sample_interval, case,
           capture_profiles=False, record=True, fit_radius=8.0):
    grid = make_grid(n)
    phi = amplitude*np.exp(-grid.r**2/(2.0*R0**2))
    pi = np.zeros_like(phi)
    phi[-1] = pi[-1] = 0.0
    t, steps, next_sample = 0.0, 0, sample_interval
    rows = [sample(t, phi, grid, fit_radius)] if record else []
    snapshots, targets = [], [4.0, 8.0, 12.0, 20.0, 40.0, 80.0]
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
        old_center, new_center = abs(phi[0]), abs(next_phi[0])
        crossing_fraction = 1.0
        if stop_level is not None and new_center >= stop_level:
            change = new_center-old_center
            crossing_fraction = np.clip((stop_level-old_center)/change, 0.0, 1.0)
            crossed = True
        if capture_profiles and new_center > old_center:
            for target in targets:
                if old_center < target <= new_center:
                    fraction = np.clip((target-old_center)/(new_center-old_center), 0.0, 1.0)
                    snap_phi = phi+fraction*(next_phi-phi)
                    snap_t = t+fraction*dt
                    sigma, _, _, _, _ = gaussian_width(snap_phi, grid.r, fit_radius)
                    for r, value in zip(grid.r, snap_phi):
                        snapshots.append({
                            "case": case, "target_abs_phi_center": target,
                            "time": float(snap_t), "radius": float(r),
                            "normalized_phi": float(value/snap_phi[0]),
                            "gaussian_fit": (float(np.exp(-r*r/(2*sigma*sigma))
                                             if np.isfinite(sigma) else np.nan)),
                            "sigma_gaussian": sigma,
                        })
        if crossed:
            t += dt*crossing_fraction
            phi += crossing_fraction*(next_phi-phi)
            pi += crossing_fraction*(next_pi-pi)
            phi[0] = np.copysign(stop_level, next_phi[0])
        else:
            t += dt
            phi, pi = next_phi, next_pi
        steps += 1
        if record and (t >= next_sample-1e-11 or crossed):
            rows.append(sample(t, phi, grid, fit_radius))
            while next_sample <= t+1e-11:
                next_sample += sample_interval

    if record and rows and rows[-1]["time"] < t-1e-11:
        rows.append(sample(t, phi, grid, fit_radius))
    return rows, snapshots, crossed, failed, (t if crossed else np.nan), steps


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = []
        for row in csv.DictReader(handle):
            for key, value in row.items():
                try:
                    row[key] = float(value)
                except (ValueError, TypeError):
                    pass
            rows.append(row)
    return rows


def fit_line(x, y):
    if len(x) < 5 or np.ptp(x) == 0:
        return np.nan, np.nan, np.nan
    slope, intercept = np.polyfit(x, y, 1)
    rss = np.sum((y-(slope*x+intercept))**2)
    total = np.sum((y-y.mean())**2)
    return float(slope), float(intercept), float(1-rss/total if total else np.nan)


def estimate_tc_joint(rows):
    """Jointly fit tc and alpha on the late, high-field central-field tail."""
    t = np.asarray([row["time"] for row in rows])
    center = np.asarray([row["abs_phi_center"] for row in rows])
    ids = np.flatnonzero(center >= 4.0)[-40:]
    if len(ids) < 8:
        return None
    fit_t, log_center = t[ids], np.log(center[ids])
    left, right = t[-1]+1e-5, t[-1]+1.0
    for _ in range(4):
        trial = np.linspace(left, right, 1201)
        log_tau = np.log(trial[:, None]-fit_t[None, :])
        centered = log_tau-log_tau.mean(axis=1, keepdims=True)
        slopes = np.sum(centered*(log_center-log_center.mean()), axis=1)/np.sum(centered**2, axis=1)
        intercepts = log_center.mean()-slopes*log_tau.mean(axis=1)
        rss = np.sum((log_center[None, :]-(slopes[:, None]*log_tau+intercepts[:, None]))**2, axis=1)
        rss[slopes >= 0] = np.inf
        index = int(np.argmin(rss))
        step = (right-left)/(len(trial)-1)
        tc, slope, best_rss = float(trial[index]), float(slopes[index]), float(rss[index])
        left, right = max(t[-1]+1e-6, tc-5*step), tc+5*step
    total = np.sum((log_center-log_center.mean())**2)
    r2 = 1.0-best_rss/total if total > 0 else np.nan
    return tc, -slope, float(r2), len(ids)


def analyze(cases):
    fits, scaling = [], []
    for case in cases:
        rows = case["rows"]
        t = np.asarray([row["time"] for row in rows])
        center = np.asarray([row["abs_phi_center"] for row in rows])
        radius = np.asarray([row["r50_gaussian"] for row in rows])
        methods = (("joint_tc_alpha", estimate_tc_joint(rows)),)
        for method, estimate in methods:
            if estimate is None:
                continue
            tc, alpha, tc_r2, tc_n = estimate
            tau = tc-t
            for lower, upper in FIT_WINDOWS:
                use = ((tau >= lower) & (tau <= upper) & (center >= 4)
                       & np.isfinite(radius) & (radius > 0))
                beta, intercept, r2 = fit_line(np.log10(tau[use]), np.log10(radius[use]))
                fits.append({
                    "case": case["case"], "N": case["n"], "amplitude": case["amplitude"],
                    "offset_percent": case["offset_percent"], "tc_method": method,
                    "tc": tc, "alpha_central": alpha, "tc_fit_r2": tc_r2,
                    "tc_fit_points": tc_n, "tau_min": lower, "tau_max": upper,
                    "n": int(use.sum()), "beta_gaussian_r50": beta, "fit_r2": r2,
                })
        ids = np.flatnonzero(center >= 4)[-40:]
        use = ids[np.isfinite(radius[ids]) & (radius[ids] > 0)]
        q, _, q_r2 = (fit_line(np.log10(center[use]), np.log10(radius[use]))
                      if len(use) >= 5 else (np.nan, np.nan, np.nan))
        quality = np.asarray([row["log_shape_r2"] for row in rows])
        core_ids = np.flatnonzero((center >= 4) & np.isfinite(quality))
        scaling.append({
            "case": case["case"], "N": case["n"], "amplitude": case["amplitude"],
            "offset_percent": case["offset_percent"], "n_width_fit": len(use),
            "central_min_fit": float(center[use[0]]) if len(use) else np.nan,
            "central_max_fit": float(center[use[-1]]) if len(use) else np.nan,
            "q_beta_over_alpha": float(-q), "q_fit_r2": q_r2,
            "median_gaussian_shape_r2_center_ge_4": (
                float(np.nanmedian(quality[core_ids])) if len(core_ids) else np.nan),
            "minimum_gaussian_shape_r2_center_ge_4": (
                float(np.nanmin(quality[core_ids])) if len(core_ids) else np.nan),
        })
    return fits, scaling


def plot_radius(cases, path):
    fig, ax = plt.subplots(figsize=(8, 5))
    for case in cases:
        rows = case["rows"]
        t = np.asarray([row["time"] for row in rows])
        center = np.asarray([row["abs_phi_center"] for row in rows])
        radius = np.asarray([row["r50_gaussian"] for row in rows])
        use = (center >= 4) & np.isfinite(radius) & (radius > 0)
        ax.plot(t[use], radius[use], label=case["case"])
    ax.set(xlabel="time t", ylabel=r"Gaussian-implied half-radius $R_{50,G}$",
           title="Gaussian core-radius estimate during collapse")
    ax.grid(alpha=0.25); ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(path, dpi=170); plt.close(fig)


def plot_loglog(cases, fits, path, method="joint_tc_alpha"):
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = plt.get_cmap("tab10").colors
    for color, case in zip(colors, cases):
        rows = case["rows"]
        estimate = estimate_tc_joint(rows)
        if estimate is None:
            continue
        tc = estimate[0]
        t = np.asarray([row["time"] for row in rows])
        center = np.asarray([row["abs_phi_center"] for row in rows])
        radius = np.asarray([row["r50_gaussian"] for row in rows])
        tau = tc-t
        use = ((tau >= 0.1) & (tau <= 1.6) & (center >= 4)
               & np.isfinite(radius) & (radius > 0))
        ax.plot(np.log10(tau[use]), np.log10(radius[use]), ".", ms=3,
                color=color, alpha=0.6, label=case["case"])
        fit = next((row for row in fits if row["case"] == case["case"]
                    and row["tc_method"] == method
                    and row["tau_min"] == 0.2 and row["tau_max"] == 0.6), None)
        if fit and np.isfinite(fit["beta_gaussian_r50"]):
            selected = use & (tau >= 0.2) & (tau <= 0.6)
            beta, intercept, _ = fit_line(np.log10(tau[selected]), np.log10(radius[selected]))
            ax.plot(np.log10(tau[selected]), np.log10(radius[selected]), "o",
                    ms=3, mfc="none", mec=color)
            x = np.linspace(np.log10(0.2), np.log10(0.6), 50)
            ax.plot(x, beta*x+intercept, "--", color=color,
                    label=rf"{case['case']}: $\beta={beta:.3f}$, $R^2={fit['fit_r2']:.3f}$")
    tc_label = r"joint high-field $t_c,\alpha$"
    ax.set(xlabel=r"$\log_{10}(t_c-t)$", ylabel=r"$\log_{10}R_{50,G}$",
           title=rf"Gaussian-width log-log fits ({tc_label})")
    ax.grid(alpha=0.25); ax.legend(fontsize=7.5); fig.tight_layout()
    fig.savefig(path, dpi=170); plt.close(fig)


def plot_quality(cases, path):
    fig, ax = plt.subplots(figsize=(8, 5))
    for case in cases:
        rows = case["rows"]
        time = np.asarray([row["time"] for row in rows])
        center = np.asarray([row["abs_phi_center"] for row in rows])
        quality = np.asarray([row["log_shape_r2"] for row in rows])
        use = (center >= 1) & np.isfinite(quality)
        ax.plot(time[use], quality[use], ".-", ms=2.5, lw=0.8,
                label=case["case"])
    ax.set(xlabel="time t", ylabel="Gaussian log-profile fit $R^2$",
           title="Gaussian fit quality as the core collapses")
    ax.grid(alpha=0.25); ax.legend(fontsize=7.5); fig.tight_layout()
    fig.savefig(path, dpi=170); plt.close(fig)


def plot_profiles(path):
    rows = read_csv(path)
    rows = [row for row in rows if row["case"] == "supercritical_plus_0.01pct_N500"]
    levels = sorted({row["target_abs_phi_center"] for row in rows})
    if not levels:
        return
    cols = min(3, len(levels)); nrows = int(np.ceil(len(levels)/cols))
    fig, axes = plt.subplots(nrows, cols, figsize=(4.2*cols, 3.0*nrows), squeeze=False)
    for ax, level in zip(axes.flat, levels):
        data = [row for row in rows if row["target_abs_phi_center"] == level
                and row["radius"] <= 8.0]
        r = np.asarray([row["radius"] for row in data])
        actual = np.asarray([abs(row["normalized_phi"]) for row in data])
        model = np.asarray([row["gaussian_fit"] for row in data])
        sigma = data[0]["sigma_gaussian"]
        ax.plot(r, actual, label="simulated profile", lw=1.4)
        ax.plot(r, model, "--", label=rf"Gaussian, $\sigma={sigma:.3g}$", lw=1.4)
        ax.set(title=rf"$|\phi(0)|={level:g}$", xlabel="r", ylabel=r"$|\phi(r)|/|\phi(0)|$",
               ylim=(-0.02, 1.05))
        ax.grid(alpha=0.2); ax.legend(fontsize=7)
    for ax in axes.flat[len(levels):]:
        ax.set_visible(False)
    fig.suptitle("Normalized radial profiles and inner Gaussian fits", y=1.01)
    fig.tight_layout(); fig.savefig(HERE/"results"/"gaussian_profile_snapshots.png",
                                    dpi=170, bbox_inches="tight"); plt.close(fig)


def load_cases():
    index = read_csv(OUT/"case_index.csv")
    return [{
        "case": row["case"], "n": int(row["N"]), "amplitude": row["amplitude"],
        "offset_percent": row["offset_percent"], "level": row["stop_level"],
        "rows": read_csv(OUT/row["csv_file"]),
    } for row in index]


def make_analysis(cases):
    fits, scaling = analyze(cases)
    write_csv(OUT/"beta_fit_sensitivity.csv", fits)
    write_csv(OUT/"width_vs_central_field.csv", scaling)
    plot_radius(cases, OUT/"gaussian_r50_vs_time.png")
    plot_quality(cases, OUT/"gaussian_fit_quality.png")
    if (OUT/"gaussian_profile_snapshots.csv").exists():
        plot_profiles(OUT/"gaussian_profile_snapshots.csv")


def search_critical(args):
    history = []

    def classify(amplitude):
        _, _, crossed, failed, event, steps = evolve(
            amplitude, args.n, args.critical_tmax, CRITICAL_LEVEL,
            args.sample_interval, "threshold_trial", record=False,
        )
        if failed:
            raise RuntimeError(f"Numerical failure at A={amplitude:.10f}")
        history.append({"amplitude": amplitude, "collapsed": int(crossed),
                        "event_time": event if crossed else "", "steps": steps})
        print(f"A={amplitude:.10f}: {'cutoff reached' if crossed else 'no cutoff'}", flush=True)
        return crossed

    lo, hi = args.a_min, args.a_max
    while classify(lo):
        lo *= 0.5
    while not classify(hi):
        hi *= 1.5
        if hi > 10:
            raise RuntimeError("Could not bracket a collapsing amplitude")
    while hi-lo > args.critical_tolerance:
        mid = 0.5*(lo+hi)
        if classify(mid):
            hi = mid
        else:
            lo = mid
    write_csv(OUT/"critical_threshold_search.csv", history)
    return 0.5*(lo+hi), lo, hi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=500)
    parser.add_argument("--t-max", type=float, default=100.0)
    parser.add_argument("--sample-interval", type=float, default=SAMPLE_INTERVAL)
    parser.add_argument("--fit-radius", type=float, default=8.0,
                        help="maximum radius used for the inner-profile fit")
    parser.add_argument("--search-critical", action="store_true",
                        help="recompute the finite-horizon Acrit bracket locally")
    parser.add_argument("--critical-tmax", type=float, default=100.0)
    parser.add_argument("--critical-tolerance", type=float, default=1e-8)
    parser.add_argument("--a-min", type=float, default=0.1)
    parser.add_argument("--a-max", type=float, default=1.2)
    parser.add_argument("--analyze-existing", action="store_true")
    args = parser.parse_args()
    if args.n < 20 or min(args.t_max, args.sample_interval, args.fit_radius) <= 0:
        parser.error("Grid size must be >=20 and time, sample spacing, fit radius must be positive")
    OUT.mkdir(parents=True, exist_ok=True)
    if args.analyze_existing:
        make_analysis(load_cases())
        print(f"Updated plots and fits in {OUT}")
        return

    if args.search_critical:
        acrit, lo, hi = search_critical(args)
        source = "recomputed N-grid finite-horizon bisection"
    else:
        acrit, lo, hi = A_CRIT_REFERENCE, np.nan, np.nan
        source = "reference N=500, t<=100 bracket from the standalone half-radius study"
    write_csv(OUT/"run_settings.csv", [{
        "Acrit_used": acrit, "critical_lower_edge": lo, "critical_upper_edge": hi,
        "critical_level": CRITICAL_LEVEL, "critical_horizon": args.critical_tmax,
        "Acrit_source": source, "N_base": args.n, "Rmax": R_MAX,
        "R0": R0, "stretch": STRETCH, "sample_interval": args.sample_interval,
        "fit_radius": args.fit_radius,
        "gaussian_fit_fraction_min": GAUSSIAN_FIT_FRACTIONS[0],
        "gaussian_fit_fraction_max": GAUSSIAN_FIT_FRACTIONS[1],
        "supercritical_cutoff": RUN_LEVEL,
    }])
    cases, all_profiles, index = [], [], []
    specs = [(OFFSETS[0], args.n, 80.0), (OFFSETS[1], args.n, RUN_LEVEL),
             (OFFSETS[2], args.n, RUN_LEVEL)]
    if args.n != 1000:
        specs.append((OFFSETS[0], 1000, 80.0))
    for offset, n, level in specs:
        percent = 100.0*offset
        suffix = f"N{n}" if n != args.n else f"N{args.n}"
        name = f"supercritical_plus_{percent:g}pct_{suffix}"
        amplitude = acrit*(1.0+offset)
        print(f"Running {name}: A={amplitude:.12f}, stop |phi0|={level:g}", flush=True)
        rows, profiles, crossed, failed, event, steps = evolve(
            amplitude, n, args.t_max, level, args.sample_interval, name,
            capture_profiles=(offset == OFFSETS[0]), fit_radius=args.fit_radius,
        )
        if failed:
            raise RuntimeError(f"Numerical failure in {name}")
        file_name = f"{name}.csv"
        write_csv(OUT/file_name, rows)
        if n == args.n and offset == OFFSETS[0]:
            all_profiles.extend(profiles)
        case = {"case": name, "n": n, "amplitude": amplitude,
                "offset_percent": percent, "level": level, "rows": rows}
        cases.append(case)
        index.append({"case": name, "N": n, "amplitude": amplitude,
                      "offset_percent": percent, "stop_level": level,
                      "crossed": int(crossed), "event_time": event,
                      "steps": steps, "csv_file": file_name})
    write_csv(OUT/"case_index.csv", index)
    write_csv(OUT/"gaussian_profile_snapshots.csv", all_profiles)
    make_analysis(cases)
    print(f"Results saved under {OUT}")


if __name__ == "__main__":
    main()
