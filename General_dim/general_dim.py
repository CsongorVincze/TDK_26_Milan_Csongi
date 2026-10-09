#!/usr/bin/env python3
"""Dimension-by-dimension critical and near-collapse scaling study."""

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
DEFAULT_DIMS = tuple(float(d) for d in range(1, 6))
CRITICAL_LEVEL = 12.0
STOP_LEVEL = 80.0
TAU_WINDOW = (0.1, 0.3)
OFFSET_PERCENT = 0.1
SAMPLE_INTERVAL = 0.01
A_CRIT_SEED_HALF_WIDTH = 0.025


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
    t_max: float = 100.0
    sample_interval: float = SAMPLE_INTERVAL


@dataclass
class Grid:
    d: float
    r: np.ndarray
    rx: np.ndarray
    rxx: np.ndarray
    damping: np.ndarray
    dx: float
    dt: float


def make_grid(c, d):
    x = np.linspace(0.0, 1.0, c.n)
    dx = x[1] - x[0]
    r = c.r_max * np.sinh(c.stretch * x) / np.sinh(c.stretch)
    rx = c.r_max * c.stretch * np.cosh(c.stretch * x) / np.sinh(c.stretch)
    rxx = c.stretch**2 * r
    sponge_start = c.sponge_fraction * c.r_max
    damping = c.sponge_strength * np.maximum(
        (r - sponge_start) / (c.r_max - sponge_start), 0.0
    ) ** c.sponge_power
    return Grid(d, r, rx, rxx, damping, dx, c.cfl * rx[0] * dx)


def laplacian(phi, grid):
    """Radial Laplacian in real dimension d on the stretched grid."""
    lap = np.zeros_like(phi)
    dr0 = grid.rx[0] * grid.dx
    lap[0] = 2.0 * grid.d * (phi[1] - phi[0]) / dr0**2
    px = (phi[2:] - phi[:-2]) / (2.0 * grid.dx)
    pxx = (phi[2:] - 2.0 * phi[1:-1] + phi[:-2]) / grid.dx**2
    lap[1:-1] = pxx / grid.rx[1:-1]**2 + (
        (grid.d - 1.0) / (grid.r[1:-1] * grid.rx[1:-1])
        - grid.rxx[1:-1] / grid.rx[1:-1]**3
    ) * px
    return lap


def rhs(phi, pi, grid):
    dphi = pi.copy()
    # Keep both terms from V'(phi)=phi-phi^3, as in Module 4.
    dpi = laplacian(phi, grid) - phi + phi**3 - grid.damping * pi
    dphi[-1] = dpi[-1] = 0.0
    return dphi, dpi


def rk4(phi, pi, dt, grid):
    k1 = rhs(phi, pi, grid)
    k2 = rhs(phi + dt*k1[0]/2.0, pi + dt*k1[1]/2.0, grid)
    k3 = rhs(phi + dt*k2[0]/2.0, pi + dt*k2[1]/2.0, grid)
    k4 = rhs(phi + dt*k3[0], pi + dt*k3[1], grid)
    new_phi = phi + dt*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
    new_pi = pi + dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6.0
    new_phi[-1] = new_pi[-1] = 0.0
    return new_phi, new_pi


def half_radius(phi, radius):
    """First radius where |phi(r)| falls to half of |phi(0)|."""
    values = np.abs(phi)
    if values[0] <= 0.0:
        return float("nan")
    crossings = np.flatnonzero(values[1:] <= 0.5 * values[0])
    if not len(crossings):
        return float("nan")
    i = int(crossings[0] + 1)
    y0, y1 = values[i-1], values[i]
    if y0 == y1:
        return float(radius[i])
    fraction = np.clip((0.5*values[0] - y0) / (y1-y0), 0.0, 1.0)
    return float(radius[i-1] + fraction*(radius[i]-radius[i-1]))


def sample(t, phi, grid):
    return {
        "d": grid.d, "time": float(t), "phi_center": float(phi[0]),
        "abs_phi_center": float(abs(phi[0])),
        "r50": half_radius(phi, grid.r),
    }


def evolve(d, amplitude, c, stop_level, record=True):
    """Evolve one run; optionally save uniform-time phi(0) and R50 samples."""
    grid = make_grid(c, d)
    phi = amplitude * np.exp(-grid.r**2 / (2.0*c.r0**2))
    pi = np.zeros_like(phi)
    phi[-1] = pi[-1] = 0.0
    t, steps = 0.0, 0
    next_sample = c.sample_interval
    rows = [sample(t, phi, grid)] if record else []
    status = "t_max_without_stop_level"

    while t < c.t_max and abs(phi[0]) < stop_level:
        peak = max(float(np.max(np.abs(phi))), 1.0)
        dt = min(grid.dt, 0.15/peak, c.t_max-t)
        if record:
            dt = min(dt, max(next_sample-t, 1e-12))
        new_phi, new_pi = rk4(phi, pi, dt, grid)
        if not np.all(np.isfinite(new_phi)) or not np.all(np.isfinite(new_pi)):
            status = "nonfinite_state"
            break

        crossed = abs(phi[0]) < stop_level <= abs(new_phi[0])
        if crossed:
            change = abs(new_phi[0]) - abs(phi[0])
            fraction = np.clip((stop_level-abs(phi[0]))/change, 0.0, 1.0)
            t += dt*fraction
            phi = phi + fraction*(new_phi-phi)
            pi = pi + fraction*(new_pi-pi)
            phi[0] = np.copysign(stop_level, new_phi[0])
            status = "reached_stop_level"
        else:
            t += dt
            phi, pi = new_phi, new_pi
        steps += 1

        if record and (t >= next_sample-1e-11 or crossed):
            rows.append(sample(t, phi, grid))
            while next_sample <= t+1e-11:
                next_sample += c.sample_interval

    if record and (not rows or rows[-1]["time"] < t-1e-11):
        rows.append(sample(t, phi, grid))
    return rows, {
        "d": float(d), "amplitude": float(amplitude), "final_time": float(t),
        "steps": steps, "max_abs_phi_center": float(abs(phi[0])),
        "status": status,
    }


def search_critical(d, c, a_min, a_max, tolerance):
    """Bisect the finite-horizon boundary for reaching |phi(0)|=12."""
    history = []

    def classify(amplitude):
        _, meta = evolve(d, amplitude, c, CRITICAL_LEVEL, record=False)
        collapsed = meta["status"] == "reached_stop_level"
        history.append({
            "d": d, "amplitude": amplitude, "reached_phi12": int(collapsed),
            "event_time": meta["final_time"] if collapsed else "",
            "horizon": c.t_max, "steps": meta["steps"],
            "status": meta["status"],
        })
        print(f"  d={d:g}, A={amplitude:.9f}: "
              f"{'reached 12' if collapsed else meta['status']}", flush=True)
        if meta["status"] == "nonfinite_state":
            raise RuntimeError(f"Numerical failure in critical search at d={d:g}, A={amplitude:g}.")
        return collapsed

    lo, hi = a_min, a_max
    while classify(lo):
        lo *= 0.5
        if lo < 1e-8:
            raise RuntimeError(f"Could not bracket a non-collapsing amplitude at d={d:g}.")
    while not classify(hi):
        hi *= 1.5
        if hi > 10.0:
            raise RuntimeError(f"Could not bracket a collapsing amplitude at d={d:g}.")
    while hi-lo > tolerance:
        mid = 0.5*(lo+hi)
        if classify(mid):
            hi = mid
        else:
            lo = mid
    return lo, hi, history


def line_fit(x, y):
    """OLS line, fit R2, and one-standard-error slope uncertainty."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    n = len(x)
    if n < 5 or np.ptp(x) == 0:
        return None
    xbar, ybar = float(x.mean()), float(y.mean())
    sxx = float(np.sum((x-xbar)**2))
    slope = float(np.sum((x-xbar)*(y-ybar))/sxx)
    intercept = ybar - slope*xbar
    residual = y - (slope*x+intercept)
    rss = float(np.sum(residual**2))
    total = float(np.sum((y-ybar)**2))
    dof = n-2
    variance = rss/dof
    slope_variance = variance/sxx
    return {
        "slope": slope, "intercept": intercept,
        "r2": 1.0-rss/total if total > 0 else float("nan"),
        "slope_se": math.sqrt(max(0.0, slope_variance)),
        "slope_variance": slope_variance, "residuals": residual,
        "sxx": sxx, "dof": dof, "n": n,
    }


def estimate_tc(rows):
    """Jointly fit tc and the central-field power law to the last 40 tail samples."""
    tail = [row for row in rows if row["abs_phi_center"] >= 4.0][-40:]
    if len(tail) < 8:
        return None
    time = np.asarray([row["time"] for row in tail])
    log_center = np.log10([row["abs_phi_center"] for row in tail])
    left, right = time[-1]+1e-5, time[-1]+1.0

    for _ in range(4):
        trial_tc = np.linspace(left, right, 1201)
        log_tau = np.log10(trial_tc[:, None] - time[None, :])
        centered_x = log_tau - log_tau.mean(axis=1, keepdims=True)
        centered_y = log_center - log_center.mean()
        slopes = np.sum(centered_x*centered_y, axis=1) / np.sum(centered_x**2, axis=1)
        intercepts = log_center.mean() - slopes*log_tau.mean(axis=1)
        residuals = log_center[None, :] - (slopes[:, None]*log_tau+intercepts[:, None])
        rss = np.sum(residuals**2, axis=1)
        rss[slopes >= 0.0] = np.inf
        index = int(np.argmin(rss))
        step = (right-left)/(len(trial_tc)-1)
        tc, alpha, best_rss = float(trial_tc[index]), float(-slopes[index]), float(rss[index])
        left, right = max(time[-1]+1e-6, tc-5*step), tc+5*step

    tau = tc-time
    tail_fit = line_fit(np.log10(tau), log_center)
    if tail_fit is None or tc <= rows[-1]["time"]:
        return None
    jacobian = np.column_stack((
        np.ones_like(tau), -np.log10(tau), -alpha/(np.log(10.0)*tau),
    ))
    dof = len(time)-3
    try:
        covariance = (best_rss/dof) * np.linalg.inv(jacobian.T@jacobian)
        tc_se = math.sqrt(max(0.0, float(covariance[2, 2])))
    except (np.linalg.LinAlgError, ZeroDivisionError):
        tc_se = float("nan")
    return {
        "tc": tc, "tc_se": tc_se, "alpha_tail": alpha,
        "tc_fit_r2": tail_fit["r2"], "tc_fit_points": len(tail),
    }


def fit_scaling(rows, tc):
    """Fit alpha and R50 beta in the requested shared log10(tau) window."""
    lower, upper = TAU_WINDOW
    times = np.asarray([row["time"] for row in rows])
    fields = np.asarray([row["abs_phi_center"] for row in rows])
    radii = np.asarray([row["r50"] for row in rows])
    tau = tc-times
    use = ((tau >= lower) & (tau <= upper) & (fields >= 4.0)
           & np.isfinite(radii) & (radii > 0.0))
    x = np.log10(tau[use])
    alpha_fit = line_fit(x, np.log10(fields[use]))
    beta_fit = line_fit(x, np.log10(radii[use]))
    if alpha_fit is None or beta_fit is None:
        return None, []

    alpha = -alpha_fit["slope"]
    beta = beta_fit["slope"]
    slope_cov = float(np.dot(alpha_fit["residuals"], beta_fit["residuals"]))
    slope_cov /= alpha_fit["dof"] * alpha_fit["sxx"]
    sum_variance = (alpha_fit["slope_variance"] + beta_fit["slope_variance"]
                    - 2.0*slope_cov)
    points = [{
        "time": float(t), "tau": float(rem), "log10_tau": float(np.log10(rem)),
        "phi_center": float(field), "log10_abs_phi_center": float(np.log10(field)),
        "r50": float(radius), "log10_r50": float(np.log10(radius)),
    } for t, rem, field, radius in zip(times[use], tau[use], fields[use], radii[use])]
    return {
        "alpha": alpha, "alpha_se": alpha_fit["slope_se"],
        "alpha_r2": alpha_fit["r2"], "beta": beta,
        "beta_se": beta_fit["slope_se"], "beta_r2": beta_fit["r2"],
        "alpha_plus_beta": alpha+beta,
        "alpha_plus_beta_se": math.sqrt(max(0.0, sum_variance)),
        "fit_points": len(points),
    }, points


def threshold_events(rows):
    """Interpolate standard central-field thresholds for the saved trajectory."""
    levels = (8.0, 10.0, 12.0, 16.0, 20.0, 30.0, 40.0, 60.0, 80.0)
    events = []
    for level in levels:
        for before, after in zip(rows[:-1], rows[1:]):
            y0, y1 = before["abs_phi_center"], after["abs_phi_center"]
            if y0 < level <= y1:
                fraction = (level-y0)/(y1-y0)
                events.append({
                    "d": before["d"], "level": level,
                    "time": before["time"] + fraction*(after["time"]-before["time"]),
                    "phi_center": before["phi_center"] + fraction*(after["phi_center"]-before["phi_center"]),
                    "r50": before["r50"] + fraction*(after["r50"]-before["r50"]),
                })
                break
    return events


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def dimension_ticks(values, max_ticks=11):
    unique = np.asarray(sorted(set(values)), dtype=float)
    indices = np.linspace(0, len(unique)-1, min(max_ticks, len(unique)))
    return unique[np.unique(np.rint(indices).astype(int))]


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key, value in row.items():
            try:
                row[key] = float(value)
            except (TypeError, ValueError):
                pass
    return rows


def load_reusable_run(path, args):
    """Load matching dimensions from a prior run, rejecting changed physics/settings."""
    settings = read_csv(path/"run_settings.csv")[0]
    expected = {
        "N": args.n, "Rmax": 25.0, "stretch": 5.0, "R0": 2.0,
        "sponge_fraction": 0.8, "sponge_strength": 7.0,
        "sponge_power": 3.0, "CFL": 0.4,
        "Acrit_horizon": args.critical_tmax, "Acrit_level": CRITICAL_LEVEL,
        "supercritical_tmax": args.t_max,
        "supercritical_offset_percent": args.offset_percent,
        "supercritical_stop_level": STOP_LEVEL,
        "sample_interval": args.sample_interval,
        "tau_min": TAU_WINDOW[0], "tau_max": TAU_WINDOW[1],
    }
    for name, value in expected.items():
        if name not in settings or not math.isclose(
                float(settings[name]), value, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"Cannot reuse {path}: setting {name} does not match.")
    files = (
        "dimension_summary.csv", "critical_summary.csv", "critical_search.csv",
        "trajectory_samples.csv", "threshold_events.csv", "scaling_fit_points.csv",
    )
    return {name[:-4]: read_csv(path/name) for name in files}


def plot_scaling(summary, path):
    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    styles = (
        ("alpha", "alpha_se", r"$\alpha$", "#7c3aed", "o"),
        ("beta", "beta_se", r"$\beta$ ($R_{50}$)", "#2563eb", "s"),
        ("alpha_plus_beta", "alpha_plus_beta_se", r"$\alpha+\beta$", "#b45309", "D"),
    )
    for value_key, error_key, label, color, marker in styles:
        rows = [row for row in summary if np.isfinite(row.get(value_key, np.nan))]
        if not rows:
            continue
        d = np.asarray([row["d"] for row in rows])
        values = np.asarray([row[value_key] for row in rows])
        errors = np.asarray([row.get(error_key, np.nan) for row in rows])
        ax.errorbar(d, values, yerr=errors, fmt=marker+"-", color=color,
                    ecolor=color, capsize=4, elinewidth=1.7, lw=1.4,
                    markersize=5, label=label)
    ax.set(xlabel="spatial dimension $d$", ylabel="fitted scaling exponent",
           title=r"Near-collapse exponents: $\tau\in[0.1,0.3]$")
    dimensions = [row["d"] for row in summary]
    ax.set_xlim(min(dimensions)-0.1, max(dimensions)+0.1)
    ax.set_xticks(dimension_ticks(dimensions))
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_collapse_times(summary, path):
    rows = [row for row in summary if np.isfinite(row.get("tc", np.nan))]
    if not rows:
        return
    d = np.asarray([row["d"] for row in rows])
    tc = np.asarray([row["tc"] for row in rows])
    errors = np.asarray([row.get("tc_se", np.nan) for row in rows])
    errors = np.where(np.isfinite(errors), errors, 0.0)
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    ax.errorbar(d, tc, yerr=errors, fmt="o-", color="#2166ac",
                ecolor="#333333", capsize=3, lw=1.5, markersize=5)
    ax.set(xlabel="spatial dimension $d$", ylabel="fitted collapse time $t_c$",
           title="Independently estimated collapse time at each dimension")
    ax.set_xlim(min(d)-0.1, max(d)+0.1)
    ax.set_xticks(dimension_ticks(d))
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_loglog_fits(summary, fit_points, path):
    """Show the two fitted log-log relations together for every dimension."""
    dimensions = [row["d"] for row in summary]
    columns = 3
    rows_count = math.ceil(len(dimensions)/columns)
    fig, axes = plt.subplots(rows_count, columns,
                             figsize=(13, 2.8*rows_count), sharex=True, squeeze=False)
    axes = axes.ravel()
    for ax, d in zip(axes, dimensions):
        points = [row for row in fit_points if row["d"] == d]
        ax_r = ax.twinx()
        if points:
            x = np.asarray([row["log10_tau"] for row in points])
            y_phi = np.asarray([row["log10_abs_phi_center"] for row in points])
            y_r50 = np.asarray([row["log10_r50"] for row in points])
            alpha = next(row["alpha"] for row in summary if row["d"] == d)
            beta = next(row["beta"] for row in summary if row["d"] == d)
            alpha_r2 = next(row["alpha_r2"] for row in summary if row["d"] == d)
            beta_r2 = next(row["beta_r2"] for row in summary if row["d"] == d)
            slope_a, intercept_a = np.polyfit(x, y_phi, 1)
            slope_b, intercept_b = np.polyfit(x, y_r50, 1)
            order = np.argsort(x)
            ax.scatter(x, y_phi, s=13, color="#7c3aed", label="central field")
            ax.plot(x[order], (slope_a*x+intercept_a)[order],
                    color="#7c3aed", lw=1.3,
                    label=rf"$\alpha={alpha:.3f}$, $R^2={alpha_r2:.3f}$")
            ax_r.scatter(x, y_r50, s=13, marker="s", color="#2563eb",
                         label=r"$R_{50}$")
            ax_r.plot(x[order], (slope_b*x+intercept_b)[order],
                      color="#2563eb", lw=1.3, ls="--",
                      label=rf"$\beta={beta:.3f}$, $R^2={beta_r2:.3f}$")
        alpha = next(row["alpha"] for row in summary if row["d"] == d)
        beta = next(row["beta"] for row in summary if row["d"] == d)
        ax.set_title(rf"$d={d:g}$: $\alpha={alpha:.3f}$, $\beta={beta:.3f}$",
                     fontsize=10)
        ax.set_xlim(np.log10(TAU_WINDOW[0]), np.log10(TAU_WINDOW[1]))
        ax.set_xlabel(r"$\log_{10}(t_c-t)$")
        ax.set_ylabel(r"$\log_{10}|\phi(0,t)|$", color="#7c3aed", fontsize=8)
        ax_r.set_ylabel(r"$\log_{10}R_{50}$", color="#2563eb", fontsize=8)
        ax.tick_params(axis="y", labelcolor="#7c3aed", labelsize=8)
        ax_r.tick_params(axis="y", labelcolor="#2563eb", labelsize=8)
        ax.tick_params(axis="x", labelsize=8)
        ax.grid(True, alpha=0.2)
    for ax in axes[len(dimensions):]:
        ax.set_visible(False)
    fig.suptitle(r"Central-field $\alpha$ and half-radius $\beta$ fits", y=0.998)
    fig.tight_layout(rect=(0, 0, 1, 0.995))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dimensions", type=float, nargs="+", default=DEFAULT_DIMS)
    parser.add_argument("--n", type=int, default=500)
    parser.add_argument("--a-min", type=float, default=0.1)
    parser.add_argument("--a-max", type=float, default=1.2)
    parser.add_argument("--critical-tmax", type=float, default=100.0)
    parser.add_argument("--critical-tolerance", type=float, default=1e-7)
    parser.add_argument("--t-max", type=float, default=100.0)
    parser.add_argument("--sample-interval", type=float, default=SAMPLE_INTERVAL)
    parser.add_argument("--offset-percent", type=float, default=OFFSET_PERCENT)
    parser.add_argument("--output-dir", type=Path, default=RESULTS)
    parser.add_argument("--reuse-from", type=Path,
                        help="reuse matching dimensions from a compatible results folder")
    args = parser.parse_args()
    if args.n < 20 or args.a_min <= 0 or args.a_max <= args.a_min:
        parser.error("Require N >= 20 and 0 < a_min < a_max.")
    if args.critical_tmax <= 0 or args.t_max <= 0 or args.sample_interval <= 0:
        parser.error("Time horizons and sample interval must be positive.")
    if args.critical_tolerance <= 0 or args.offset_percent <= 0:
        parser.error("Critical tolerance and supercritical offset must be positive.")
    if not args.dimensions or any(d <= 0 for d in args.dimensions):
        parser.error("Dimensions must be positive real values.")
    try:
        cache = load_reusable_run(args.reuse_from, args) if args.reuse_from else None
    except (OSError, ValueError, IndexError) as error:
        parser.error(str(error))

    config = Config(n=args.n, t_max=args.critical_tmax,
                    sample_interval=args.sample_interval)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    critical_history, critical_summary, trajectories = [], [], []
    all_events = []
    fit_points, summary = [], []
    cached_by_d = ({round(float(row["d"]), 8): row
                    for row in cache["dimension_summary"]} if cache else {})

    for d in args.dimensions:
        key = round(float(d), 8)
        if key in cached_by_d:
            summary.append(cached_by_d[key])
            for name, target in (
                    ("critical_summary", critical_summary),
                    ("critical_search", critical_history),
                    ("trajectory_samples", trajectories),
                    ("threshold_events", all_events),
                    ("scaling_fit_points", fit_points)):
                target.extend(row for row in cache[name]
                              if round(float(row["d"]), 8) == key)
            print(f"d={d:g}: reusing matching saved run", flush=True)
            continue

        search_min, search_max = args.a_min, args.a_max
        if cache:
            anchors = sorted(cache["dimension_summary"], key=lambda row: row["d"])
            anchor_d = np.asarray([row["d"] for row in anchors], dtype=float)
            if anchor_d[0] < d < anchor_d[-1]:
                anchor_a = np.asarray([row["Acrit"] for row in anchors], dtype=float)
                seed = float(np.interp(d, anchor_d, anchor_a))
                search_min = max(1e-8, seed-A_CRIT_SEED_HALF_WIDTH)
                search_max = seed+A_CRIT_SEED_HALF_WIDTH
        print(f"Searching A_crit for d={d:g} (finite horizon {config.t_max:g}, "
              f"initial bracket [{search_min:.6f}, {search_max:.6f}])", flush=True)
        lo, hi, history = search_critical(
            float(d), config, search_min, search_max, args.critical_tolerance)
        critical_history.extend(history)
        acrit = 0.5*(lo+hi)
        amplitude = acrit*(1.0+args.offset_percent/100.0)
        run_config = Config(n=args.n, t_max=args.t_max,
                            sample_interval=args.sample_interval)
        rows, meta = evolve(float(d), amplitude, run_config, STOP_LEVEL, record=True)
        all_events.extend(threshold_events(rows))
        for row in rows:
            trajectories.append({
                **row, "Acrit": acrit, "amplitude": amplitude,
                "offset_percent": args.offset_percent,
            })

        tc_fit = estimate_tc(rows) if meta["status"] == "reached_stop_level" else None
        scaling, used_points = (fit_scaling(rows, tc_fit["tc"])
                                if tc_fit else (None, []))
        fit_points.extend({"d": float(d), **point} for point in used_points)
        record = {
            "d": float(d), "N": args.n, "r0": run_config.r0,
            "critical_horizon": args.critical_tmax,
            "critical_level": CRITICAL_LEVEL,
            "Acrit_lower": lo, "Acrit_upper": hi, "Acrit": acrit,
            "offset_percent": args.offset_percent, "amplitude": amplitude,
            "stop_level": STOP_LEVEL, "run_status": meta["status"],
            "collapse_event_time": meta["final_time"] if tc_fit else "",
            "tc": tc_fit["tc"] if tc_fit else float("nan"),
            "tc_se": tc_fit["tc_se"] if tc_fit else float("nan"),
            "tc_fit_r2": tc_fit["tc_fit_r2"] if tc_fit else float("nan"),
            "tc_fit_points": tc_fit["tc_fit_points"] if tc_fit else 0,
            **(scaling or {
                "alpha": float("nan"), "alpha_se": float("nan"),
                "alpha_r2": float("nan"), "beta": float("nan"),
                "beta_se": float("nan"), "beta_r2": float("nan"),
                "alpha_plus_beta": float("nan"),
                "alpha_plus_beta_se": float("nan"), "fit_points": 0,
            }),
        }
        summary.append(record)
        critical_summary.append({
            "d": float(d), "N": args.n, "horizon": args.critical_tmax,
            "collapse_level": CRITICAL_LEVEL, "Acrit_lower": lo,
            "Acrit_upper": hi, "Acrit": acrit,
            "bracket_width": hi-lo,
            "supercritical_offset_percent": args.offset_percent,
            "supercritical_amplitude": amplitude,
        })
        print(f"d={d:g}: Acrit={acrit:.9f}, A={amplitude:.9f}, "
              f"{meta['status']}, t_end={meta['final_time']:.6f}, "
              f"tc={record['tc']:.6f}, alpha={record['alpha']:.5f}, "
              f"beta={record['beta']:.5f}", flush=True)

    write_csv(args.output_dir/"critical_search.csv", critical_history)
    write_csv(args.output_dir/"critical_summary.csv", critical_summary)
    write_csv(args.output_dir/"trajectory_samples.csv", trajectories)
    write_csv(args.output_dir/"dimension_summary.csv", summary)
    write_csv(args.output_dir/"threshold_events.csv", all_events)
    write_csv(args.output_dir/"scaling_fit_points.csv", fit_points)
    write_csv(args.output_dir/"run_settings.csv", [{
        "dimensions": ";".join(f"{d:g}" for d in args.dimensions),
        "N": args.n, "Rmax": 25.0, "stretch": 5.0, "R0": 2.0,
        "sponge_fraction": 0.8, "sponge_strength": 7.0,
        "sponge_power": 3.0, "CFL": 0.4,
        "initial_velocity": 0.0, "Acrit_horizon": args.critical_tmax,
        "supercritical_tmax": args.t_max,
        "Acrit_level": CRITICAL_LEVEL, "supercritical_offset_percent": args.offset_percent,
        "supercritical_stop_level": STOP_LEVEL,
        "sample_interval": args.sample_interval,
        "critical_tolerance": args.critical_tolerance,
        "reused_results_from": str(args.reuse_from) if args.reuse_from else "",
        "intermediate_search_seed": (
            f"linear interpolation of reused Acrit values; +/-{A_CRIT_SEED_HALF_WIDTH:g} bracket"
            if cache else "global amplitude bracket from --a-min and --a-max"),
        "seed_bracket_half_width": A_CRIT_SEED_HALF_WIDTH if cache else "",
        "tau_min": TAU_WINDOW[0], "tau_max": TAU_WINDOW[1],
        "radius_definition": "first r where |phi(r,t)| = 0.5*|phi(0,t)|",
        "tc_estimator": "joint free-tc central-field power-law fit to last 40 samples with |phi0|>=4",
    }])
    plot_scaling(summary, args.output_dir/"alpha_beta_sum_vs_dimension.png")
    plot_loglog_fits(summary, fit_points, args.output_dir/"loglog_fits_by_dimension.png")
    plot_collapse_times(summary, args.output_dir/"collapse_time_vs_dimension.png")


if __name__ == "__main__":
    main()
