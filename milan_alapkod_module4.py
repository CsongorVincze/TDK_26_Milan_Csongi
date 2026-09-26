#!/usr/bin/env python3
"""Module 4 simulation: threshold search and collapse-core contraction.

This is a focused copy of ``milan_alapkod.py`` for the single-well model
V(phi) = phi**2/2 - phi**4/4.  It keeps the stretched spherical grid and
sponge layer, and adds a reproducible threshold search plus the energy-
weighted core-radius analysis requested in Oscillon_dynamics_notes.

Only NumPy and Matplotlib are required.  Run with ``--help`` for controls.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass(frozen=True)
class Config:
    """Numerical and physical parameters for one simulation."""

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
    collapse_threshold: float = 12.0


def make_grid(config: Config):
    """Build the same sinh-stretched radial grid used by the base script."""
    x = np.linspace(0.0, 1.0, config.n)
    dx = x[1] - x[0]
    sinh_stretch = np.sinh(config.stretch)
    r = config.r_max * np.sinh(config.stretch * x) / sinh_stretch
    dr_dx = config.r_max * config.stretch * np.cosh(
        config.stretch * x
    ) / sinh_stretch
    d2r_dx2 = config.r_max * config.stretch**2 * np.sinh(
        config.stretch * x
    ) / sinh_stretch

    sponge_start = config.sponge_fraction * config.r_max
    gamma = np.zeros_like(r)
    sponge_mask = r > sponge_start
    gamma[sponge_mask] = config.sponge_strength * (
        (r[sponge_mask] - sponge_start)
        / (config.r_max - sponge_start)
    ) ** config.sponge_power

    grid = {
        "x": x,
        "dx": dx,
        "r": r,
        "dr_dx": dr_dx,
        "d2r_dx2": d2r_dx2,
        "gamma": gamma,
        "dr_min": dr_dx[0] * dx,
        "base_dt": config.cfl * dr_dx[0] * dx,
    }
    return grid


def potential(phi):
    """Module 4 single-well potential V(phi) = phi^2/2 - phi^4/4."""
    return 0.5 * phi**2 - 0.25 * phi**4


def potential_derivative(phi):
    """Derivative V'(phi) = phi - phi^3."""
    return phi - phi**3


def radial_laplacian(phi, grid):
    """Spherically symmetric 3D Laplacian on the stretched radial grid."""
    lap = np.zeros_like(phi)
    dx = grid["dx"]
    r = grid["r"]
    dr_dx = grid["dr_dx"]
    d2r_dx2 = grid["d2r_dx2"]

    # For a smooth spherical field, Delta(phi)(0) = 3 phi_rr(0).
    lap[0] = 6.0 * (phi[1] - phi[0]) / (dr_dx[0] * dx) ** 2

    phi_x = (phi[2:] - phi[:-2]) / (2.0 * dx)
    phi_xx = (phi[2:] - 2.0 * phi[1:-1] + phi[:-2]) / dx**2
    lap[1:-1] = (
        phi_xx / dr_dx[1:-1] ** 2
        + (
            2.0 / (r[1:-1] * dr_dx[1:-1])
            - d2r_dx2[1:-1] / dr_dx[1:-1] ** 3
        )
        * phi_x
    )

    # Homogeneous Dirichlet boundary: the boundary value is held at zero.
    lap[-1] = 0.0
    return lap


def rhs(phi, pi, grid):
    """Right-hand side of phi_tt - Delta(phi) + V'(phi) = 0 with sponge."""
    dphi = pi.copy()
    dpi = radial_laplacian(phi, grid) - potential_derivative(phi)
    dpi -= grid["gamma"] * pi
    dphi[-1] = 0.0
    dpi[-1] = 0.0
    return dphi, dpi


def rk4_step(phi, pi, dt, grid):
    """One classical fourth-order Runge-Kutta step for (phi, pi)."""
    k1_phi, k1_pi = rhs(phi, pi, grid)
    k2_phi, k2_pi = rhs(
        phi + 0.5 * dt * k1_phi,
        pi + 0.5 * dt * k1_pi,
        grid,
    )
    k3_phi, k3_pi = rhs(
        phi + 0.5 * dt * k2_phi,
        pi + 0.5 * dt * k2_pi,
        grid,
    )
    k4_phi, k4_pi = rhs(phi + dt * k3_phi, pi + dt * k3_pi, grid)

    phi_next = phi + (dt / 6.0) * (
        k1_phi + 2.0 * k2_phi + 2.0 * k3_phi + k4_phi
    )
    pi_next = pi + (dt / 6.0) * (
        k1_pi + 2.0 * k2_pi + 2.0 * k3_pi + k4_pi
    )
    phi_next[-1] = 0.0
    pi_next[-1] = 0.0
    return phi_next, pi_next


def initial_state(amplitude, config, grid):
    """Gaussian initial data phi=A exp(-r^2/(2 R0^2)), pi=0."""
    phi = amplitude * np.exp(-grid["r"] ** 2 / (2.0 * config.r0**2))
    pi = np.zeros_like(phi)
    phi[-1] = 0.0
    return phi, pi


def core_energy_moment(phi, pi, config, grid):
    """Return core energy and literal energy-weighted <r^2> from the notes.

    The notes leave R_core unspecified.  This implementation uses a fixed
    spherical core of radius ``config.core_radius``; 4*pi cancels in the ratio.
    Negative local potential energy is retained as written in the definition.
    """
    x = grid["x"]
    r = grid["r"]
    dphi_dr = np.gradient(phi, x, edge_order=2) / grid["dr_dx"]
    dphi_dr[0] = 0.0
    density = 0.5 * pi**2 + 0.5 * dphi_dr**2 + potential(phi)
    mask = r <= config.core_radius
    rc = r[mask]
    rho = density[mask]
    denominator = np.trapezoid(rho * rc**2, x=rc)
    numerator = np.trapezoid(rho * rc**4, x=rc)
    absolute_scale = np.trapezoid(np.abs(rho) * rc**2, x=rc)
    if abs(denominator) <= 1.0e-12 * max(absolute_scale, 1.0):
        return float("nan"), float("nan")
    mean_r2 = numerator / denominator
    if not np.isfinite(mean_r2) or mean_r2 <= 0.0:
        mean_r2 = float("nan")
    return 4.0 * np.pi * denominator, float(mean_r2)


def sample_state(t, phi, pi, config, grid):
    """Make one compact diagnostic row."""
    core_energy, mean_r2 = core_energy_moment(phi, pi, config, grid)
    return (
        float(t),
        float(phi[0]),
        float(pi[0]),
        core_energy,
        mean_r2,
        float(np.sqrt(mean_r2)) if np.isfinite(mean_r2) else float("nan"),
    )


def simulate(amplitude, config, grid, t_max, store=False, progress=False):
    """Evolve one Gaussian and stop at collapse proxy or ``t_max``.

    The time step obeys the grid CFL limit and is reduced as the field grows
    to resolve the increasingly fast local nonlinear timescale.  Collapse is
    classified when the central field first reaches ``collapse_threshold``.
    """
    phi, pi = initial_state(amplitude, config, grid)
    t = 0.0
    collapsed = bool(abs(phi[0]) >= config.collapse_threshold)
    next_sample = config.sample_interval
    samples = [sample_state(t, phi, pi, config, grid)] if store else []
    last_report = time.perf_counter()
    steps = 0

    while t < t_max and not collapsed:
        max_phi = float(np.max(np.abs(phi)))
        dt = min(grid["base_dt"], 0.15 / max(max_phi, 1.0))
        dt = min(dt, t_max - t)
        if store:
            dt = min(dt, max(next_sample - t, 1.0e-12))
        old_center = float(abs(phi[0]))
        phi_next, pi_next = rk4_step(phi, pi, dt, grid)
        t_next = t + dt

        if not np.all(np.isfinite(phi_next)) or not np.all(np.isfinite(pi_next)):
            collapsed = True
            t = t_next
            phi, pi = phi_next, pi_next
            break

        crossed = abs(phi_next[0]) >= config.collapse_threshold
        if crossed:
            rise = abs(phi_next[0]) - old_center
            fraction = (
                (config.collapse_threshold - old_center) / rise
                if rise > 0.0
                else 1.0
            )
            fraction = float(np.clip(fraction, 0.0, 1.0))
            t = t + dt * fraction
            # Store a linearly interpolated state at the same event time.
            # This keeps the last CSV row's timestamp and field synchronized.
            phi = phi + fraction * (phi_next - phi)
            pi = pi + fraction * (pi_next - pi)
            phi[0] = np.copysign(config.collapse_threshold, phi_next[0])
        else:
            t = t_next
            phi, pi = phi_next, pi_next
        steps += 1

        if store and (t >= next_sample - 1.0e-12 or crossed):
            samples.append(sample_state(t, phi, pi, config, grid))
            while next_sample <= t + 1.0e-12:
                next_sample += config.sample_interval

        collapsed = crossed or float(np.max(np.abs(phi))) >= config.collapse_threshold
        if progress and time.perf_counter() - last_report >= 1.0:
            print(
                f"  A={amplitude:.6f}: t={t:.2f}/{t_max:g}, "
                f"phi(0,t)={phi[0]:.4g}, steps={steps}",
                flush=True,
            )
            last_report = time.perf_counter()

    if store and (not samples or samples[-1][0] < t - 1.0e-12):
        samples.append(sample_state(t, phi, pi, config, grid))
    return {
        "amplitude": float(amplitude),
        "collapsed": bool(collapsed),
        "event_time": float(t) if collapsed else float("nan"),
        "t_end": float(t),
        "steps": steps,
        "samples": np.asarray(samples, dtype=float),
    }


def find_threshold(config, grid, t_max, a_min, a_max, tolerance):
    """Bisection on the collapse/no-collapse classification."""
    history = []

    def classify(amplitude):
        result = simulate(amplitude, config, grid, t_max=t_max)
        row = {
            "amplitude": float(amplitude),
            "collapsed": bool(result["collapsed"]),
            "event_time": result["event_time"],
            "steps": result["steps"],
        }
        history.append(row)
        outcome = "collapse" if result["collapsed"] else "no collapse by t_max"
        print(f"  A={amplitude:.7f}: {outcome}", flush=True)
        return result["collapsed"]

    lower = float(a_min)
    while classify(lower):
        lower *= 0.5
        if lower < 1.0e-6:
            raise RuntimeError("Could not find a non-collapsing lower bracket.")

    upper = float(a_max)
    while not classify(upper):
        upper *= 1.5
        if upper > 10.0:
            raise RuntimeError("Could not find a collapsing upper bracket.")

    while upper - lower > tolerance:
        middle = 0.5 * (lower + upper)
        if classify(middle):
            upper = middle
        else:
            lower = middle
    return lower, upper, history


def estimate_collapse_time(samples, collapse_threshold):
    """Estimate t_c by linear extrapolation of 1/phi(0,t) to zero.

    This uses the cubic-wave near-collapse scaling phi(0,t) ~ (t_c-t)^-1.
    The fit is restricted to the increasing high-amplitude tail.
    """
    t = samples[:, 0]
    center = np.abs(samples[:, 1])
    indices = np.flatnonzero(center >= max(2.0, 0.2 * collapse_threshold))
    if indices.size < 5:
        raise RuntimeError("Not enough high-amplitude samples to estimate t_c.")
    indices = indices[-min(40, indices.size):]
    fit_t = t[indices]
    inverse_phi = 1.0 / center[indices]
    slope, intercept = np.polyfit(fit_t, inverse_phi, 1)
    if not np.isfinite(slope) or slope >= 0.0:
        raise RuntimeError("The inverse-amplitude fit did not approach blowup.")
    t_c = -intercept / slope
    if t_c <= fit_t[-1]:
        raise RuntimeError("The extrapolated collapse time is not after the data.")
    fitted = slope * fit_t + intercept
    residual = np.sum((inverse_phi - fitted) ** 2)
    total = np.sum((inverse_phi - np.mean(inverse_phi)) ** 2)
    fit_r2 = 1.0 - residual / total if total > 0.0 else float("nan")
    return float(t_c), float(fit_r2), indices


def fit_contraction(samples, t_c, config, grid, fit_points):
    """Fit R_energy ~ (t_c-t)^beta on the resolved, contracting tail."""
    t = samples[:, 0]
    center = np.abs(samples[:, 1])
    radius = samples[:, 5]
    tau = t_c - t
    mask = (
        np.isfinite(radius)
        & (radius >= 3.0 * grid["dr_min"])
        & (tau >= 2.0 * grid["base_dt"])
        & (tau <= 1.0)
        & (center >= 1.2)
    )
    indices = np.flatnonzero(mask)
    if indices.size < 5:
        raise RuntimeError(
            "Fewer than five resolved energy-radius samples are available "
            "for the requested log-log fit. Try a larger supercritical "
            "amplitude, smaller --sample-interval, or larger --core-radius."
        )
    indices = indices[-fit_points:]
    log_tau = np.log(tau[indices])
    log_radius = np.log(radius[indices])
    beta, log_scale = np.polyfit(log_tau, log_radius, 1)
    fitted = beta * log_tau + log_scale
    residual = np.sum((log_radius - fitted) ** 2)
    total = np.sum((log_radius - np.mean(log_radius)) ** 2)
    fit_r2 = 1.0 - residual / total if total > 0.0 else float("nan")
    return float(beta), float(log_scale), float(fit_r2), indices


def write_threshold_csv(path, history):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("amplitude", "collapsed", "event_time", "steps"),
        )
        writer.writeheader()
        writer.writerows(history)


def write_collapse_csv(path, samples):
    columns = (
        "time",
        "phi_center",
        "pi_center",
        "core_energy",
        "mean_r2",
        "core_rms_radius",
    )
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(samples)


def plot_results(output_dir, history, samples, t_c, beta, beta_r2,
                 log_scale, fit_indices):
    """Save threshold, central-field, and log-log contraction plots."""
    figure, axis = plt.subplots(figsize=(7.5, 4.5))
    amplitudes = [row["amplitude"] for row in history]
    collapse_times = [
        row["event_time"] if row["collapsed"] else np.nan for row in history
    ]
    axis.scatter(amplitudes, collapse_times, color="navy", s=24)
    axis.set_xlabel("initial amplitude $A_0$")
    axis.set_ylabel("time to collapse proxy")
    axis.set_title("Module 4 critical-amplitude bisection")
    axis.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(output_dir / "module4_threshold_search.png", dpi=160)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(7.5, 4.5))
    axis.plot(samples[:, 0], samples[:, 1], color="navy", linewidth=1.0)
    axis.axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
    axis.axvline(t_c, color="firebrick", linestyle=":", linewidth=1.0,
                 label=rf"$t_c={t_c:.5f}$")
    axis.set_xlabel("time $t$")
    axis.set_ylabel(r"central field $\phi(0,t)$")
    axis.set_title("Supercritical central-field collapse")
    axis.legend()
    axis.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(output_dir / "module4_central_field.png", dpi=160)
    plt.close(figure)

    selected = samples[fit_indices]
    tau = t_c - selected[:, 0]
    radius = selected[:, 5]
    order = np.argsort(tau)
    tau = tau[order]
    radius = radius[order]
    fit_tau = np.geomspace(float(np.min(tau)), float(np.max(tau)), 200)

    figure, axis = plt.subplots(figsize=(7.5, 4.8))
    axis.loglog(tau, radius, "o", ms=4, label="energy-weighted core radius")
    fit_is_reliable = beta > 0.0 and beta_r2 >= 0.9
    fit_label = (
        rf"power-law fit: $\beta={beta:.4f}$, $R^2={beta_r2:.3f}$"
        if fit_is_reliable
        else rf"OLS diagnostic: $\beta={beta:.4f}$, $R^2={beta_r2:.3f}$"
    )
    axis.loglog(
        fit_tau,
        np.exp(log_scale) * fit_tau**beta,
        "--",
        label=fit_label,
    )
    axis.set_xlabel(r"remaining time $t_c-t$")
    axis.set_ylabel(r"$\sqrt{\langle r^2\rangle}$")
    axis.set_title("Module 4 core-contraction exponent")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_dir / "module4_core_contraction.png", dpi=160)
    plt.close(figure)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("full", "threshold", "collapse"),
                        default="full")
    parser.add_argument("--output-dir", type=Path, default=Path("module4_results"))
    parser.add_argument("--n", type=int, default=500,
                        help="radial grid points; default matches the base script")
    parser.add_argument("--r-max", type=float, default=25.0)
    parser.add_argument("--r0", type=float, default=2.0,
                        help="Module 4 initial Gaussian radius")
    parser.add_argument("--core-radius", type=float, default=8.0,
                        help="fixed R_core in the energy-weighted moment")
    parser.add_argument("--threshold-tmax", type=float, default=40.0,
                        help="classification horizon for threshold bisection")
    parser.add_argument("--collapse-tmax", type=float, default=40.0,
                        help="maximum duration of the supercritical run")
    parser.add_argument("--collapse-threshold", type=float, default=12.0,
                        help="central-field value used to identify collapse")
    parser.add_argument("--sample-interval", type=float, default=0.02)
    parser.add_argument("--a-min", type=float, default=0.1)
    parser.add_argument("--a-max", type=float, default=1.2)
    parser.add_argument("--threshold-tolerance", type=float, default=5.0e-5)
    parser.add_argument("--supercritical-amplitude", type=float, default=None,
                        help="override the automatic amplitude above Acrit")
    parser.add_argument("--fit-points", type=int, default=30,
                        help="maximum resolved samples in the beta fit")
    return parser.parse_args()


def main():
    args = parse_args()
    config = Config(
        n=args.n,
        r_max=args.r_max,
        r0=args.r0,
        core_radius=args.core_radius,
        sample_interval=args.sample_interval,
        collapse_threshold=args.collapse_threshold,
    )
    if config.n < 20 or config.r0 <= 0.0 or config.core_radius <= 0.0:
        raise SystemExit("Require n >= 20 and positive R0 and core radius.")
    if config.core_radius >= config.r_max:
        raise SystemExit("core radius must be smaller than r_max.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    grid = make_grid(config)

    print("Module 4: V(phi)=phi^2/2-phi^4/4, R0=%.4g" % config.r0)
    print(
        f"grid: N={config.n}, r_max={config.r_max:g}, "
        f"dr_min={grid['dr_min']:.6g}, base_dt={grid['base_dt']:.6g}; "
        f"R_core={config.core_radius:g}"
    )

    threshold_history = []
    if args.mode in ("full", "threshold"):
        print("Searching for the subcritical/supercritical boundary:")
        a_sub, a_super, threshold_history = find_threshold(
            config,
            grid,
            args.threshold_tmax,
            args.a_min,
            args.a_max,
            args.threshold_tolerance,
        )
        a_crit = 0.5 * (a_sub + a_super)
        write_threshold_csv(
            args.output_dir / "module4_threshold_search.csv",
            threshold_history,
        )
        print(
            f"Acrit = {a_crit:.4f} (bisection bracket "
            f"[{a_sub:.7f}, {a_super:.7f}])"
        )
        if args.mode == "threshold":
            summary = {
                "Acrit_midpoint": a_crit,
                "subcritical_bracket": a_sub,
                "supercritical_bracket": a_super,
                "threshold_definition": (
                    "lower amplitude does not reach |phi(0,t)|="
                    "collapse_threshold by threshold_tmax; upper amplitude does"
                ),
                "threshold_tmax": args.threshold_tmax,
                "collapse_threshold": config.collapse_threshold,
                "config": asdict(config),
            }
            (args.output_dir / "module4_threshold_summary.json").write_text(
                json.dumps(summary, indent=2) + "\n", encoding="utf-8"
            )
            return
        amplitude = args.supercritical_amplitude
        if amplitude is None:
            amplitude = a_super + max(0.02, 0.05 * a_crit)
    else:
        if args.supercritical_amplitude is None:
            raise SystemExit("--mode collapse requires --supercritical-amplitude.")
        amplitude = args.supercritical_amplitude
        a_crit = float("nan")
        a_sub = float("nan")
        a_super = float("nan")

    print(f"Running supercritical collapse case at A0={amplitude:.7f}:")
    result = simulate(
        amplitude,
        config,
        grid,
        t_max=args.collapse_tmax,
        store=True,
        progress=True,
    )
    if not result["collapsed"]:
        raise RuntimeError(
            f"A0={amplitude:.7f} did not reach |phi(0,t)|="
            f"{config.collapse_threshold:g} by t={args.collapse_tmax:g}. "
            "Increase --collapse-tmax or --supercritical-amplitude."
        )

    samples = result["samples"]
    write_collapse_csv(args.output_dir / "module4_collapse_data.csv", samples)
    t_c, tc_r2, tc_indices = estimate_collapse_time(
        samples, config.collapse_threshold
    )
    beta, log_scale, beta_r2, fit_indices = fit_contraction(
        samples, t_c, config, grid, args.fit_points
    )
    plot_results(
        args.output_dir,
        threshold_history,
        samples,
        t_c,
        beta,
        beta_r2,
        log_scale,
        fit_indices,
    )

    fit_is_reliable = beta > 0.0 and beta_r2 >= 0.9
    summary = {
        "Acrit_midpoint": a_crit,
        "Acrit_bracket": [a_sub, a_super],
        "Acrit_operational_definition": (
            "lower amplitude does not reach |phi(0,t)|=collapse_threshold "
            "by threshold_tmax; upper amplitude does"
        ),
        "threshold_classification_horizon": args.threshold_tmax,
        "supercritical_amplitude": amplitude,
        "collapse_threshold_phi_center": config.collapse_threshold,
        "collapse_proxy_time": result["event_time"],
        "estimated_tc_from_inverse_center_field": t_c,
        "tc_estimation_method": (
            "linear extrapolation of 1/abs(phi(0,t)) over the last "
            "high-amplitude samples"
        ),
        "tc_fit_r_squared": tc_r2,
        "beta_energy_weighted_core_radius": beta if fit_is_reliable else None,
        "beta_raw_loglog_slope_diagnostic": beta,
        "beta_fit_r_squared": beta_r2,
        "beta_fit_status": (
            "accepted power-law fit"
            if fit_is_reliable
            else "undetermined: signed energy-weighted radius does not show "
                 "a reliable positive-slope power law"
        ),
        "beta_fit_sample_count": int(len(fit_indices)),
        "beta_fit_time_window": [
            float(samples[fit_indices[0], 0]),
            float(samples[fit_indices[-1], 0]),
        ],
        "beta_fit_tau_window": [
            float(np.min(t_c - samples[fit_indices, 0])),
            float(np.max(t_c - samples[fit_indices, 0])),
        ],
        "beta_fit_r_squared_acceptance_threshold": 0.9,
        "tc_fit_sample_count": int(len(tc_indices)),
        "core_radius_cutoff": config.core_radius,
        "config": asdict(config),
    }
    (args.output_dir / "module4_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    print(f"collapse proxy time: {result['event_time']:.6f}")
    print(f"estimated tc: {t_c:.6f} (inverse-amplitude fit R^2={tc_r2:.5f})")
    if fit_is_reliable:
        print(f"beta: {beta:.4f} (log-log fit R^2={beta_r2:.5f})")
    else:
        print(
            f"beta: undetermined; formal log-log slope={beta:.4f}, "
            f"R^2={beta_r2:.5f} (no reliable positive-slope power law)"
        )
    print(f"Wrote results to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
