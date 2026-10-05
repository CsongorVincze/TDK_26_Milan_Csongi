#!/usr/bin/env python3
"""Solve Module 4's radial single-well PDE and analyze collapse."""

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator, PercentFormatter
import numpy as np


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
    core_radius: float = 8.0
    sample_interval: float = 0.02
    collapse_level: float = 12.0


@dataclass
class Grid:
    x: np.ndarray
    r: np.ndarray
    rx: np.ndarray
    rxx: np.ndarray
    gamma: np.ndarray
    dx: float
    dr_min: float
    dt: float


def make_grid(c):
    x = np.linspace(0.0, 1.0, c.n)
    dx = x[1] - x[0]
    r = c.r_max * np.sinh(c.stretch * x) / np.sinh(c.stretch)
    rx = c.r_max * c.stretch * np.cosh(c.stretch * x) / np.sinh(c.stretch)
    rxx = c.stretch**2 * r
    start = c.sponge_fraction * c.r_max
    gamma = c.sponge_strength * np.maximum(
        (r - start) / (c.r_max - start), 0.0
    ) ** c.sponge_power
    dr_min = rx[0] * dx
    return Grid(x, r, rx, rxx, gamma, dx, dr_min, c.cfl * dr_min)


def laplacian(phi, g):
    """Radial 3D Laplacian; the origin uses Delta(phi)(0)=3 phi_rr(0)."""
    lap = np.zeros_like(phi)
    lap[0] = 6.0 * (phi[1] - phi[0]) / (g.rx[0] * g.dx) ** 2
    px = (phi[2:] - phi[:-2]) / (2.0 * g.dx)
    pxx = (phi[2:] - 2.0 * phi[1:-1] + phi[:-2]) / g.dx**2
    lap[1:-1] = pxx / g.rx[1:-1]**2 + (
        2.0 / (g.r[1:-1] * g.rx[1:-1])
        - g.rxx[1:-1] / g.rx[1:-1]**3
    ) * px
    return lap


def rhs(phi, pi, g):
    dphi = pi.copy()
    dpi = laplacian(phi, g) - (phi - phi**3) - g.gamma * pi
    dphi[-1] = dpi[-1] = 0.0
    return dphi, dpi


def rk4(phi, pi, dt, g):
    k1 = rhs(phi, pi, g)
    k2 = rhs(phi + dt*k1[0]/2, pi + dt*k1[1]/2, g)
    k3 = rhs(phi + dt*k2[0]/2, pi + dt*k2[1]/2, g)
    k4 = rhs(phi + dt*k3[0], pi + dt*k3[1], g)
    out_phi = phi + dt*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6
    out_pi = pi + dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6
    out_phi[-1] = out_pi[-1] = 0.0
    return out_phi, out_pi


def moment(phi, pi, c, g):
    """Return total/core energy and the notes' signed energy-weighted <r^2>."""
    dphi_dr = np.gradient(phi, g.x, edge_order=2) / g.rx
    dphi_dr[0] = 0.0
    rho = 0.5*pi**2 + 0.5*dphi_dr**2 + 0.5*phi**2 - 0.25*phi**4
    total_energy = 4*np.pi*np.trapezoid(rho*g.r**2, x=g.r)
    use = g.r <= c.core_radius
    r, rho = g.r[use], rho[use]
    den = np.trapezoid(rho*r**2, x=r)
    num = np.trapezoid(rho*r**4, x=r)
    scale = np.trapezoid(np.abs(rho)*r**2, x=r)
    if abs(den) <= 1e-12 * max(scale, 1.0):
        return float(total_energy), float("nan"), float("nan")
    r2 = num / den
    if not np.isfinite(r2) or r2 <= 0:
        r2 = float("nan")
    return float(total_energy), 4*np.pi*den, float(r2)


def sample(t, phi, pi, c, g):
    total_energy, core_energy, r2 = moment(phi, pi, c, g)
    return (t, phi[0], pi[0], total_energy, core_energy, r2,
            np.sqrt(r2) if np.isfinite(r2) else float("nan"))


def simulate(a, c, g, t_max, save=False):
    phi = a * np.exp(-g.r**2 / (2*c.r0**2))
    pi = np.zeros_like(phi)
    phi[-1] = 0.0
    t, steps, next_sample = 0.0, 0, c.sample_interval
    rows = [sample(t, phi, pi, c, g)] if save else []
    collapsed = abs(phi[0]) >= c.collapse_level

    while t < t_max and not collapsed:
        peak = max(float(np.max(np.abs(phi))), 1.0)
        dt = min(g.dt, 0.15/peak, t_max - t)
        if save:
            dt = min(dt, max(next_sample - t, 1e-12))
        new_phi, new_pi = rk4(phi, pi, dt, g)
        if not np.all(np.isfinite(new_phi)) or not np.all(np.isfinite(new_pi)):
            return True, t + dt, steps + 1, np.asarray(rows)

        crossed = abs(new_phi[0]) >= c.collapse_level
        if crossed:
            rise = abs(new_phi[0]) - abs(phi[0])
            f = np.clip((c.collapse_level - abs(phi[0])) / rise, 0.0, 1.0)
            t += dt*f
            phi += f*(new_phi - phi)
            pi += f*(new_pi - pi)
            phi[0] = np.copysign(c.collapse_level, new_phi[0])
        else:
            t += dt
            phi, pi = new_phi, new_pi
        steps += 1

        collapsed = crossed
        if save and (t >= next_sample - 1e-12 or crossed):
            rows.append(sample(t, phi, pi, c, g))
            while next_sample <= t + 1e-12:
                next_sample += c.sample_interval

    if save and rows[-1][0] < t - 1e-12:
        rows.append(sample(t, phi, pi, c, g))
    return collapsed, t if collapsed else float("nan"), steps, np.asarray(rows)


def find_threshold(c, g, args):
    history = []

    def classify(a):
        collapsed, tc, steps, _ = simulate(a, c, g, args.threshold_tmax)
        history.append((a, collapsed, tc, steps))
        print(f"A={a:.7f}: {'collapse' if collapsed else 'no collapse by t_max'}")
        return collapsed

    lo = args.a_min
    while classify(lo):
        lo *= 0.5
        if lo < 1e-6:
            raise RuntimeError("Could not find a non-collapsing lower bracket.")
    hi = args.a_max
    while not classify(hi):
        hi *= 1.5
        if hi > 10:
            raise RuntimeError("Could not find a collapsing upper bracket.")
    while hi - lo > args.threshold_tolerance:
        mid = (lo + hi)/2
        if classify(mid):
            hi = mid
        else:
            lo = mid
    return lo, hi, history


def line_fit(x, y):
    slope, intercept = np.polyfit(x, y, 1)
    residual = np.sum((y - (slope*x + intercept))**2)
    total = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - residual/total if total > 0 else float("nan")
    return float(slope), float(intercept), float(r2)


def estimate_tc(data, level):
    t, center = data[:, 0], np.abs(data[:, 1])
    ids = np.flatnonzero(center >= max(2.0, 0.2*level))[-40:]
    if len(ids) < 5:
        raise RuntimeError("Too few high-amplitude samples to estimate t_c.")
    slope, intercept, r2 = line_fit(t[ids], 1.0/center[ids])
    tc = -intercept/slope
    if slope >= 0 or tc <= t[ids[-1]]:
        raise RuntimeError("Inverse-amplitude fit did not extrapolate to collapse.")
    return float(tc), r2, ids


def fit_beta(data, tc, c, g, count):
    tau = tc - data[:, 0]
    radius = data[:, 6]
    use = (np.isfinite(radius) & (radius >= 3*g.dr_min)
           & (tau >= 2*g.dt) & (tau <= 1.0)
           & (np.abs(data[:, 1]) >= 1.2))
    ids = np.flatnonzero(use)[-count:]
    if len(ids) < 5:
        raise RuntimeError("Too few resolved core-radius samples for the beta fit.")
    # Fit the stated relation directly in explicit base-10 log coordinates.
    beta, log_scale, r2 = line_fit(np.log10(tau[ids]), np.log10(radius[ids]))
    return beta, log_scale, r2, ids


def save_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


TIME_SERIES_HEADER = (
    "time", "phi_center", "pi_center", "total_energy", "core_energy",
    "mean_r2", "core_rms_radius",
)


def project_path(path):
    """Use a project-relative path in summaries when possible."""
    path = path.resolve()
    try:
        return str(path.relative_to(Path(__file__).resolve().parent.parent))
    except ValueError:
        return str(path)


def save_threshold_plot(out, history, filename="module4_threshold_search.png"):
    if not history:
        return
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.scatter([h[0] for h in history],
               [h[2] if h[1] else np.nan for h in history], s=22)
    ax.set(xlabel="initial amplitude $A_0$",
           ylabel="time to collapse proxy", title="Critical-amplitude search")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / filename, dpi=150)
    plt.close(fig)


def save_critical_plot(out, midpoint, upper_edge, acrit, amplitude,
                       collapse_level, midpoint_event, event_time):
    """Compare the threshold estimate with the nearby collapsing trajectory."""
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].plot(midpoint[:, 0], midpoint[:, 1], color="navy",
                 label=rf"$A_0=A_{{crit}}={acrit:.10f}$")
    axes[0].plot(upper_edge[:, 0], upper_edge[:, 1], color="darkorange",
                 label=rf"supercritical bracket edge $A_0={amplitude:.10f}$")
    axes[0].axhline(collapse_level, color="firebrick", linestyle=":")
    axes[0].axhline(-collapse_level, color="firebrick", linestyle=":")
    if np.isfinite(midpoint_event):
        axes[0].axvline(midpoint_event, color="navy", linestyle="--",
                        label=f"Acrit collapse: t={midpoint_event:.4f}")
    if np.isfinite(event_time):
        axes[0].axvline(event_time, color="darkorange", linestyle="--",
                        label=f"upper-edge collapse: t={event_time:.4f}")
    axes[0].legend(fontsize=8)
    axes[0].set_ylabel(r"central field $\phi(0,t)$")
    axes[0].set_title("Critical-amplitude estimate and nearby collapse")
    axes[1].semilogy(
        midpoint[:, 0], np.maximum(np.abs(midpoint[:, 1]), 1e-8),
        color="navy", label=r"$A_0=A_{crit}$",
    )
    axes[1].semilogy(
        upper_edge[:, 0], np.maximum(np.abs(upper_edge[:, 1]), 1e-8),
        color="darkorange", label="supercritical bracket edge",
    )
    axes[1].set(xlabel="time $t$", ylabel=r"$|\phi(0,t)|$ (log scale)")
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "module4_critical_central_field.png", dpi=150)
    plt.close(fig)


def save_plots(out, history, data, tc, beta, log_scale, beta_r2, ids):
    save_threshold_plot(out, history)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(data[:, 0], data[:, 1], color="navy")
    ax.axvline(tc, color="firebrick", linestyle=":", label=rf"$t_c={tc:.5f}$")
    ax.set(xlabel="time $t$", ylabel=r"central field $\phi(0,t)$",
           title="Supercritical central field")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "module4_central_field.png", dpi=150)
    plt.close(fig)

    tau = tc - data[ids, 0]
    radius = data[ids, 6]
    log_tau, log_radius = np.log10(tau), np.log10(radius)
    order = np.argsort(log_tau)
    log_tau, log_radius = log_tau[order], log_radius[order]
    reliable = beta > 0 and beta_r2 >= 0.9
    label = (rf"power law: $\beta={beta:.4f}$, $R^2={beta_r2:.3f}$"
             if reliable else rf"OLS diagnostic: $\beta={beta:.4f}$, $R^2={beta_r2:.3f}$")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    fit_x = np.linspace(log_tau.min(), log_tau.max(), 150)
    ax.plot(log_tau, log_radius, "o", ms=4, label="fit samples")
    ax.plot(fit_x, beta*fit_x + log_scale, "--", label=label)
    ax.set(xlabel=r"$\log_{10}(t_c-t)$",
           ylabel=r"$\log_{10}(\sqrt{\langle r^2\rangle})$",
           title="Core-radius scaling (explicit log coordinates)")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "module4_core_contraction.png", dpi=150)
    plt.close(fig)


def run_critical_mode(args, c, g):
    """Search a long-horizon threshold, then plot its collapsing bracket edge."""
    args.threshold_tmax = args.critical_tmax
    args.threshold_tolerance = args.critical_tolerance
    lo, hi, history = find_threshold(c, g, args)
    acrit = (lo + hi) / 2
    amplitude = hi  # upper bracket edge is the closest tested collapsing value
    print(f"Acrit≈{acrit:.10f}; bracket=[{lo:.10f}, {hi:.10f}]")
    print(f"Running the supercritical bracket edge A0={amplitude:.10f}.")

    save_csv(args.output_dir / "module4_critical_threshold_search.csv",
             ("amplitude", "collapsed", "event_time", "steps"), history)
    save_threshold_plot(args.output_dir, history,
                        "module4_critical_threshold_search.png")
    midpoint_run = simulate(acrit, c, g, args.critical_tmax, save=True)
    upper_run = simulate(amplitude, c, g, args.critical_tmax, save=True)
    midpoint_collapsed, midpoint_time, midpoint_steps, midpoint_data = midpoint_run
    collapsed, event_time, steps, data = upper_run
    save_csv(args.output_dir / "module4_critical_midpoint_data.csv",
             TIME_SERIES_HEADER, midpoint_data)
    csv_path = args.output_dir / "module4_critical_data.csv"
    save_csv(csv_path, TIME_SERIES_HEADER, data)
    save_critical_plot(args.output_dir, midpoint_data, data, acrit, amplitude,
                       c.collapse_level, midpoint_time, event_time)

    summary = {
        "Acrit_estimate": acrit,
        "Acrit_bracket": [lo, hi],
        "run_amplitude": amplitude,
        "run_amplitude_offset": amplitude - acrit,
        "threshold_horizon": args.critical_tmax,
        "threshold_tolerance": args.critical_tolerance,
        "collapse_level": c.collapse_level,
        "Acrit_midpoint_collapse_observed": bool(midpoint_collapsed),
        "Acrit_midpoint_collapse_time": (
            float(midpoint_time) if midpoint_collapsed else None
        ),
        "Acrit_midpoint_integration_steps": int(midpoint_steps),
        "collapse_observed": bool(collapsed),
        "collapse_proxy_time": float(event_time) if collapsed else None,
        "integration_steps": int(steps),
        "csv_file": project_path(csv_path),
        "Acrit_midpoint_csv_file": project_path(
            args.output_dir / "module4_critical_midpoint_data.csv"
        ),
    }
    (args.output_dir / "module4_critical_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    midpoint_status = (f"collapse at t={midpoint_time:.5f}" if midpoint_collapsed
                       else "no collapse by the time limit")
    edge_status = (f"collapse at t={event_time:.5f}" if collapsed
                   else "no collapse by the time limit")
    print(f"Acrit midpoint: {midpoint_status}; upper bracket edge: {edge_status}.")
    print(f"Critical-run data and plot saved in {args.output_dir.resolve()}")


def save_supercritical_sweep_plot(out, runs, acrit, collapse_level, t_max):
    """Plot the central-field runs and collapse time across the amplitude sweep."""
    offsets = [run["relative_offset"] for run in runs]
    cmap = plt.get_cmap("viridis")
    norm = LogNorm(vmin=min(offsets), vmax=max(offsets))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    event_limit = max(
        run["event_time"] if run["collapsed"] else t_max for run in runs
    )

    for run in runs:
        color = cmap(norm(run["relative_offset"]))
        data = run["data"]
        axes[0].plot(data[:, 0], data[:, 1], color=color, lw=1.0)
        axes[1].scatter(
            100 * run["relative_offset"],
            run["event_time"] if run["collapsed"] else t_max,
            color=color, marker="o" if run["collapsed"] else "x", s=34,
            zorder=3,
        )

    axes[0].axhline(collapse_level, color="firebrick", ls=":", lw=1)
    axes[0].axhline(-collapse_level, color="firebrick", ls=":", lw=1)
    axes[0].set(
        xlabel="time $t$", ylabel=r"central field $\phi(0,t)$",
        title=rf"Supercritical runs ($A_{{crit}}\approx{acrit:.8f}$)",
    )
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar = fig.colorbar(mappable, ax=axes[0], pad=0.02)
    colorbar.set_label(r"relative excess $(A-A_{crit})/A_{crit}$ (%)")
    colorbar.set_ticks(offsets)
    colorbar.set_ticklabels([f"{100*offset:g}" for offset in offsets])

    axes[1].set_xscale("log")
    axes[1].set(
        xlabel=r"amplitude excess $(A-A_{crit})/A_{crit}$ (%)",
        ylabel="time to collapse cutoff",
        title="Collapse time versus supercritical offset",
        ylim=(0, event_limit * 1.12),
    )
    axes[1].xaxis.set_major_formatter(
        PercentFormatter(xmax=100, decimals=2)
    )
    if any(not run["collapsed"] for run in runs):
        axes[1].text(
            0.03, 0.04, r"$\times$: cutoff not reached by $t_{max}$",
            transform=axes[1].transAxes, fontsize=8,
        )
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "module4_supercritical_sweep.png", dpi=160)
    plt.close(fig)


def save_beta_axes(ax, run, g, show_legend=False):
    """Draw one log-log core-radius fit on an existing axis."""
    data, tc = run["data"], run["tc"]
    tau = tc - data[:, 0]
    radius = data[:, 6]
    context = (
        np.isfinite(radius) & (radius >= 3*g.dr_min)
        & (tau >= 2*g.dt) & (tau <= 1.0)
        & (np.abs(data[:, 1]) >= 1.2)
    )
    fit_ids = run["fit_ids"]
    fit_tau, fit_radius = tau[fit_ids], radius[fit_ids]
    order = np.argsort(fit_tau)
    fit_tau, fit_radius = fit_tau[order], fit_radius[order]
    log_tau = np.log10(tau[context])
    log_radius = np.log10(radius[context])
    fit_x, fit_y = np.log10(fit_tau), np.log10(fit_radius)
    line_x = np.linspace(fit_x.min(), fit_x.max(), 100)
    color = "#2563eb"
    ax.plot(log_tau, log_radius, ".", ms=3, color="#94a3b8",
            alpha=0.6, label="eligible near-collapse samples")
    ax.plot(fit_x, fit_y, "o", ms=3.8, color=color,
            label="points used in fit")
    ax.plot(
        line_x, run["beta_ols"]*line_x + run["log_scale"],
        "--", color="#dc2626", lw=1.6,
        label=rf"fit, $\beta={run['beta_ols']:.4f}$",
    )
    status = "accepted" if run["beta_accepted"] else "undetermined"
    ax.set_title(
        rf"$A_0={run['amplitude']:.6f}$ (+{100*run['relative_offset']:g}%)",
        fontsize=9,
    )
    ax.text(
        0.04, 0.96,
        rf"$t_c={tc:.5f}$" + "\n"
        + rf"$\beta={run['beta_ols']:.4f},\ R^2={run['beta_r2']:.3f}$" + "\n"
        + f"{status}; n={len(fit_ids)}",
        transform=ax.transAxes, va="top", fontsize=7.5,
        bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "#cbd5e1", "alpha": 0.88},
    )
    ax.grid(True, which="both", alpha=0.22)
    if show_legend:
        ax.legend(fontsize=7.5, loc="lower right")


def save_supercritical_beta_plots(out, runs, acrit, g):
    """Save one beta-fit panel per run, plus a compact comparison grid."""
    for run in runs:
        if not run["beta_available"]:
            continue
        fig, ax = plt.subplots(figsize=(6.2, 4.6))
        save_beta_axes(ax, run, g, show_legend=True)
        ax.set(
            xlabel=r"$\log_{10}(t_c-t)$",
            ylabel=r"$\log_{10}(\sqrt{\langle r^2\rangle})$",
            title=ax.get_title() + "\n"
            + rf"$A_{{crit}}={acrit:.8f}$; inverse-amplitude $R^2={run['tc_r2']:.3f}$"
            + "\nexplicit base-10 log coordinates",
        )
        fig.tight_layout()
        fig.savefig(out / run["beta_plot_file"].split("/")[-1], dpi=160)
        plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(14, 7), sharex=True, sharey=True)
    for ax, run in zip(axes.flat, runs):
        if run["beta_available"]:
            save_beta_axes(ax, run, g)
        else:
            ax.set_title(f"A0={run['amplitude']:.6f}")
            ax.text(0.5, 0.5, run["beta_status"], ha="center", va="center",
                    transform=ax.transAxes, fontsize=8)
            ax.grid(True, alpha=0.22)
        ax.set_xticks((-0.6, -0.4, -0.2, 0.0))
        ax.xaxis.set_minor_locator(NullLocator())
    fig.legend(
        handles=(
            Line2D([], [], marker=".", ls="", color="#94a3b8"),
            Line2D([], [], marker="o", ls="", color="#2563eb"),
            Line2D([], [], ls="--", color="#dc2626"),
        ),
        labels=("eligible samples", "fit points", "OLS power law"),
        loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=3,
        frameon=False, fontsize=8,
    )
    fig.supxlabel(r"$\log_{10}(t_c-t)$", y=0.07)
    fig.supylabel(r"$\log_{10}(\sqrt{\langle r^2\rangle})$")
    fig.suptitle("Supercritical core-radius fits in explicit base-10 log coordinates", y=0.99)
    fig.tight_layout(rect=(0.04, 0.13, 1, 0.95))
    fig.savefig(out / "module4_supercritical_beta_fits.png", dpi=170)
    plt.close(fig)


def run_supercritical_sweep(args, c, g):
    """Search Acrit, simulate a short logarithmic sweep above it, and save CSVs."""
    if args.critical_amplitude is None:
        args.threshold_tmax = args.critical_tmax
        args.threshold_tolerance = args.critical_tolerance
        lo, hi, history = find_threshold(c, g, args)
        acrit = (lo + hi) / 2
        bracket = [lo, hi]
        acrit_source = "fresh threshold search"
        save_csv(
            args.output_dir / "module4_supercritical_threshold_search.csv",
            ("amplitude", "collapsed", "event_time", "steps"), history,
        )
    else:
        acrit, bracket, history = args.critical_amplitude, None, []
        acrit_source = "supplied estimate"
        saved_critical = args.output_dir / "module4_critical_summary.json"
        if saved_critical.exists():
            prior = json.loads(saved_critical.read_text(encoding="utf-8"))
            if np.isclose(prior.get("Acrit_estimate", np.nan), acrit,
                          rtol=0, atol=1e-12):
                bracket = prior.get("Acrit_bracket")
                acrit_source = project_path(saved_critical)
        print(f"Using supplied Acrit={acrit:.10f}; skipping threshold search.")
    offsets = sorted(set(args.sweep_relative_offsets))

    runs = []
    for index, offset in enumerate(offsets, start=1):
        amplitude = acrit * (1 + offset)
        collapsed, event_time, steps, data = simulate(
            amplitude, c, g, args.sweep_tmax, save=True,
        )
        csv_path = args.output_dir / (
            f"module4_supercritical_{index:02d}_A{amplitude:.9f}.csv"
        )
        save_csv(csv_path, TIME_SERIES_HEADER, data)
        run = {
            "relative_offset": offset, "amplitude": amplitude,
            "collapsed": collapsed, "event_time": event_time,
            "steps": steps, "csv_path": csv_path, "data": data,
            "beta_available": False, "beta_accepted": False,
        }
        if collapsed:
            try:
                tc, tc_r2, _ = estimate_tc(data, c.collapse_level)
                beta, log_scale, beta_r2, fit_ids = fit_beta(
                    data, tc, c, g, args.fit_points,
                )
                run.update({
                    "beta_available": True, "tc": tc, "tc_r2": tc_r2,
                    "beta_ols": beta, "log_scale": log_scale,
                    "beta_r2": beta_r2, "fit_ids": fit_ids,
                    "beta_accepted": beta > 0 and beta_r2 >= 0.9,
                    "fit_start": float(data[fit_ids[0], 0]),
                    "fit_end": float(data[fit_ids[-1], 0]),
                    "beta_status": (
                        "accepted" if beta > 0 and beta_r2 >= 0.9
                        else "undetermined; no reliable positive-slope fit"
                    ),
                })
            except RuntimeError as error:
                run["beta_status"] = f"fit unavailable: {error}"
        else:
            run["beta_status"] = "not fit: collapse cutoff was not reached"
        run["beta_plot_file"] = (
            f"module4_supercritical_beta_{index:02d}_A{amplitude:.9f}.png"
            if run["beta_available"] else None
        )
        runs.append(run)
        status = f"collapse at t={event_time:.5f}" if collapsed else "no cutoff reached"
        print(
            f"Sweep {index}/{len(offsets)}: A={amplitude:.9f} "
            f"(+{100*offset:.4g}%): {status} ({steps} steps); "
            f"beta {run['beta_status']}"
        )

    summary_path = args.output_dir / "module4_supercritical_sweep_summary.csv"
    save_csv(
        summary_path,
        ("relative_offset", "relative_offset_percent", "amplitude",
         "collapsed", "event_time", "steps", "csv_file"),
        ((run["relative_offset"], 100*run["relative_offset"],
          run["amplitude"], run["collapsed"],
          run["event_time"] if run["collapsed"] else "",
          run["steps"], project_path(run["csv_path"])) for run in runs),
    )
    save_supercritical_sweep_plot(
        args.output_dir, runs, acrit, c.collapse_level, args.sweep_tmax,
    )
    save_supercritical_beta_plots(args.output_dir, runs, acrit, g)
    beta_summary_path = args.output_dir / "module4_supercritical_beta_summary.csv"
    save_csv(
        beta_summary_path,
        ("relative_offset_percent", "amplitude", "collapse_time", "estimated_tc",
         "tc_fit_r2", "beta", "beta_ols_diagnostic", "beta_fit_r2",
         "beta_fit_points", "fit_start", "fit_end", "beta_status",
         "beta_plot_file", "time_series_csv"),
        ((100*run["relative_offset"], run["amplitude"],
          run["event_time"] if run["collapsed"] else "",
          run.get("tc", ""), run.get("tc_r2", ""),
          run.get("beta_ols", "") if run["beta_accepted"] else "",
          run.get("beta_ols", ""), run.get("beta_r2", ""),
          len(run.get("fit_ids", ())), run.get("fit_start", ""),
          run.get("fit_end", ""), run["beta_status"],
          run["beta_plot_file"] or "", project_path(run["csv_path"]))
         for run in runs),
    )
    summary = {
        "Acrit_estimate": acrit, "Acrit_bracket": bracket,
        "Acrit_source": acrit_source,
        "threshold_horizon": args.critical_tmax,
        "threshold_tolerance": args.critical_tolerance,
        "sweep_horizon": args.sweep_tmax,
        "collapse_level": c.collapse_level,
        "relative_offsets": offsets,
        "runs": [
            {
                "relative_offset": float(run["relative_offset"]),
                "amplitude": float(run["amplitude"]),
                "collapsed": bool(run["collapsed"]),
                "event_time": (float(run["event_time"])
                               if run["collapsed"] else None),
                "steps": int(run["steps"]),
                "csv_file": project_path(run["csv_path"]),
                "estimated_tc": float(run["tc"]) if run["beta_available"] else None,
                "tc_fit_r2": float(run["tc_r2"]) if run["beta_available"] else None,
                "beta": float(run["beta_ols"])
                        if run["beta_available"] and run["beta_accepted"] else None,
                "beta_ols_diagnostic": float(run["beta_ols"])
                                       if run["beta_available"] else None,
                "beta_fit_r2": float(run["beta_r2"])
                               if run["beta_available"] else None,
                "beta_fit_points": len(run.get("fit_ids", ())),
                "beta_status": run["beta_status"],
                "beta_plot_file": run["beta_plot_file"],
            } for run in runs
        ],
        "summary_csv": project_path(summary_path),
        "plot_file": project_path(args.output_dir / "module4_supercritical_sweep.png"),
        "beta_summary_csv": project_path(beta_summary_path),
        "beta_comparison_plot": project_path(
            args.output_dir / "module4_supercritical_beta_fits.png"
        ),
        "config": asdict(c),
    }
    (args.output_dir / "module4_supercritical_sweep_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8",
    )
    print(f"Sweep summary and plot saved in {args.output_dir.resolve()}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--mode", choices=("full", "threshold", "collapse", "critical", "sweep"),
        default="full",
    )
    default_results = Path(__file__).resolve().parent / "results"
    p.add_argument("--output-dir", type=Path, default=default_results)
    p.add_argument("--csv-name", default="module4_collapse_data.csv")
    p.add_argument("--n", type=int, default=500)
    p.add_argument("--r-max", type=float, default=25.0)
    p.add_argument("--r0", type=float, default=2.0)
    p.add_argument("--core-radius", type=float, default=8.0)
    p.add_argument("--threshold-tmax", type=float, default=40.0)
    p.add_argument("--collapse-tmax", type=float, default=40.0)
    p.add_argument("--critical-tmax", type=float, default=100.0)
    p.add_argument("--critical-tolerance", type=float, default=1e-8)
    p.add_argument("--sweep-tmax", type=float, default=100.0)
    p.add_argument("--critical-amplitude", type=float,
                   help="use an existing Acrit estimate instead of searching again")
    p.add_argument(
        "--sweep-relative-offsets", type=float, nargs="+",
        default=(1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.10, 0.25),
        metavar="FRACTION",
    )
    p.add_argument("--collapse-threshold", type=float, default=12.0)
    p.add_argument("--sample-interval", type=float, default=0.02)
    p.add_argument("--a-min", type=float, default=0.1)
    p.add_argument("--a-max", type=float, default=1.2)
    p.add_argument("--threshold-tolerance", type=float, default=5e-5)
    p.add_argument("--supercritical-amplitude", type=float)
    p.add_argument("--fit-points", type=int, default=30)
    args = p.parse_args()

    if args.mode == "critical" and (
        args.critical_tmax <= 0 or args.critical_tolerance <= 0
    ):
        p.error("--critical-tmax and --critical-tolerance must be positive.")
    if args.mode == "sweep" and (
        args.critical_tmax <= 0 or args.critical_tolerance <= 0
        or args.sweep_tmax <= 0 or len(set(args.sweep_relative_offsets)) < 2
        or any(offset <= 0 for offset in args.sweep_relative_offsets)
        or (args.critical_amplitude is not None and args.critical_amplitude <= 0)
    ):
        p.error(
            "Sweep requires positive horizons/tolerance and at least two "
            "distinct positive relative offsets."
        )

    c = Config(n=args.n, r_max=args.r_max, r0=args.r0,
               core_radius=args.core_radius,
               sample_interval=args.sample_interval,
               collapse_level=args.collapse_threshold)
    if c.n < 20 or c.r0 <= 0 or not 0 < c.core_radius < c.r_max:
        p.error("Require n >= 20 and 0 < R0, R_core < r_max.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    g = make_grid(c)
    history = []

    if args.mode == "critical":
        run_critical_mode(args, c, g)
        return
    if args.mode == "sweep":
        run_supercritical_sweep(args, c, g)
        return

    if args.mode != "collapse":
        lo, hi, history = find_threshold(c, g, args)
        acrit = (lo + hi)/2
        print(f"Acrit={acrit:.4f}; bracket=[{lo:.7f}, {hi:.7f}]")
        save_csv(args.output_dir / "module4_threshold_search.csv",
                 ("amplitude", "collapsed", "event_time", "steps"), history)
        if args.mode == "threshold":
            result = {"Acrit": acrit, "bracket": [lo, hi],
                      "horizon": args.threshold_tmax,
                      "collapse_level": c.collapse_level,
                      "config": asdict(c)}
            (args.output_dir / "module4_threshold_summary.json").write_text(
                json.dumps(result, indent=2) + "\n", encoding="utf-8")
            return
        a = args.supercritical_amplitude
        if a is None:
            a = hi + max(0.02, 0.05*acrit)
    else:
        if args.supercritical_amplitude is None:
            p.error("--mode collapse requires --supercritical-amplitude.")
        a = args.supercritical_amplitude
        lo = hi = acrit = float("nan")

    collapsed, event_time, steps, data = simulate(
        a, c, g, args.collapse_tmax, save=True)
    if not collapsed:
        raise RuntimeError(f"A0={a:g} did not reach |phi(0,t)|="
                           f"{c.collapse_level:g} by t={args.collapse_tmax:g}.")
    csv_path = args.output_dir / args.csv_name
    save_csv(csv_path, TIME_SERIES_HEADER, data)
    tc, tc_r2, tc_ids = estimate_tc(data, c.collapse_level)
    beta, log_scale, beta_r2, ids = fit_beta(data, tc, c, g, args.fit_points)
    save_plots(args.output_dir, history, data, tc, beta, log_scale, beta_r2, ids)

    accepted = beta > 0 and beta_r2 >= 0.9
    summary = {
        "Acrit": acrit, "Acrit_bracket": [lo, hi],
        "threshold_horizon": args.threshold_tmax,
        "collapse_level": c.collapse_level,
        "supercritical_amplitude": a,
        "collapse_proxy_time": event_time,
        "estimated_tc": tc, "tc_fit_r2": tc_r2,
        "beta": beta if accepted else None,
        "beta_ols_diagnostic": beta, "beta_fit_r2": beta_r2,
        "beta_status": "accepted" if accepted else "undetermined; no reliable positive-slope power law",
        "beta_fit_samples": len(ids),
        "beta_fit_time_window": [float(data[ids[0], 0]), float(data[ids[-1], 0])],
        "core_radius_cutoff": c.core_radius,
        "csv_file": project_path(csv_path),
        "config": asdict(c), "integration_steps": steps,
    }
    (args.output_dir / "module4_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")
    print(f"tc={tc:.6f} (fit R2={tc_r2:.4f})")
    if accepted:
        print(f"beta={beta:.4f} (fit R2={beta_r2:.4f})")
    else:
        print(f"beta undetermined; slope={beta:.4f}, fit R2={beta_r2:.4f}")
    print(f"Results: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
