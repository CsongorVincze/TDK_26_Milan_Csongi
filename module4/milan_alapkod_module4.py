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
    beta, log_scale, r2 = line_fit(np.log(tau[ids]), np.log(radius[ids]))
    return beta, log_scale, r2, ids


def save_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def save_plots(out, history, data, tc, beta, log_scale, beta_r2, ids):
    if history:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.scatter([h[0] for h in history],
                   [h[2] if h[1] else np.nan for h in history], s=22)
        ax.set(xlabel="initial amplitude $A_0$",
               ylabel="time to collapse proxy",
               title="Critical-amplitude search")
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(out / "module4_threshold_search.png", dpi=150)
        plt.close(fig)

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
    order = np.argsort(tau)
    tau, radius = tau[order], radius[order]
    fit_tau = np.geomspace(tau.min(), tau.max(), 150)
    reliable = beta > 0 and beta_r2 >= 0.9
    label = (rf"power law: $\beta={beta:.4f}$, $R^2={beta_r2:.3f}$"
             if reliable else rf"OLS diagnostic: $\beta={beta:.4f}$, $R^2={beta_r2:.3f}$")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(tau, radius, "o", ms=4, label="energy-weighted core radius")
    ax.loglog(fit_tau, np.exp(log_scale)*fit_tau**beta, "--", label=label)
    ax.set(xlabel=r"remaining time $t_c-t$",
           ylabel=r"$\sqrt{\langle r^2\rangle}$",
           title="Core-radius scaling")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "module4_core_contraction.png", dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=("full", "threshold", "collapse"), default="full")
    default_results = Path(__file__).resolve().parent / "results"
    p.add_argument("--output-dir", type=Path, default=default_results)
    p.add_argument("--csv-name", default="module4_collapse_data.csv")
    p.add_argument("--n", type=int, default=500)
    p.add_argument("--r-max", type=float, default=25.0)
    p.add_argument("--r0", type=float, default=2.0)
    p.add_argument("--core-radius", type=float, default=8.0)
    p.add_argument("--threshold-tmax", type=float, default=40.0)
    p.add_argument("--collapse-tmax", type=float, default=40.0)
    p.add_argument("--collapse-threshold", type=float, default=12.0)
    p.add_argument("--sample-interval", type=float, default=0.02)
    p.add_argument("--a-min", type=float, default=0.1)
    p.add_argument("--a-max", type=float, default=1.2)
    p.add_argument("--threshold-tolerance", type=float, default=5e-5)
    p.add_argument("--supercritical-amplitude", type=float)
    p.add_argument("--fit-points", type=int, default=30)
    args = p.parse_args()

    c = Config(n=args.n, r_max=args.r_max, r0=args.r0,
               core_radius=args.core_radius,
               sample_interval=args.sample_interval,
               collapse_level=args.collapse_threshold)
    if c.n < 20 or c.r0 <= 0 or not 0 < c.core_radius < c.r_max:
        p.error("Require n >= 20 and 0 < R0, R_core < r_max.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    g = make_grid(c)
    history = []

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
    save_csv(csv_path,
             ("time", "phi_center", "pi_center", "total_energy", "core_energy",
              "mean_r2", "core_rms_radius"), data)
    csv_reference = csv_path.resolve()
    try:
        project_root = Path(__file__).resolve().parent.parent
        csv_reference = csv_reference.relative_to(project_root)
    except ValueError:
        pass
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
        "csv_file": str(csv_reference),
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
