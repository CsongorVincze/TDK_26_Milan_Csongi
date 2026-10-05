#!/usr/bin/env python3
"""Animate the near-critical subcritical radial profile and core diagnostics."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

from milan_alapkod_module4 import Config, make_grid, rk4, sample


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE / "visualizations" / "subcritical_near_critical"
HEADER = (
    "time", "phi_center", "pi_center", "total_energy", "core_energy",
    "mean_r2", "core_rms_radius",
)


def evolve(amplitude, config, grid, t_max, sample_dt, profile_dt):
    """Evolve the same RK4 system, retaining sparse radial-profile snapshots."""
    phi = amplitude * np.exp(-grid.r**2 / (2*config.r0**2))
    pi = np.zeros_like(phi)
    phi[-1] = pi[-1] = 0.0
    t, steps = 0.0, 0
    rows = [sample(t, phi, pi, config, grid)]
    profile_times, profiles = [t], [phi.copy()]
    next_sample, next_profile = sample_dt, profile_dt
    cutoff_reached = abs(phi[0]) >= config.collapse_level

    while t < t_max and not cutoff_reached:
        peak = max(float(np.max(np.abs(phi))), 1.0)
        dt = min(grid.dt, 0.15/peak, t_max-t)
        dt = min(dt, max(next_sample-t, 1e-12),
                 max(next_profile-t, 1e-12))
        new_phi, new_pi = rk4(phi, pi, dt, grid)
        if not np.all(np.isfinite(new_phi)) or not np.all(np.isfinite(new_pi)):
            break

        cutoff_reached = abs(new_phi[0]) >= config.collapse_level
        if cutoff_reached:
            rise = abs(new_phi[0]) - abs(phi[0])
            fraction = np.clip(
                (config.collapse_level-abs(phi[0]))/rise, 0.0, 1.0,
            )
            t += dt*fraction
            phi += fraction*(new_phi-phi)
            pi += fraction*(new_pi-pi)
            phi[0] = np.copysign(config.collapse_level, new_phi[0])
        else:
            t += dt
            phi, pi = new_phi, new_pi
        steps += 1

        if t >= next_sample-1e-11 or cutoff_reached:
            rows.append(sample(t, phi, pi, config, grid))
            while next_sample <= t+1e-11:
                next_sample += sample_dt
        if t >= next_profile-1e-11 or cutoff_reached:
            profile_times.append(t)
            profiles.append(phi.copy())
            while next_profile <= t+1e-11:
                next_profile += profile_dt

    if rows[-1][0] < t-1e-11:
        rows.append(sample(t, phi, pi, config, grid))
    if profile_times[-1] < t-1e-11:
        profile_times.append(t)
        profiles.append(phi.copy())
    return (np.asarray(rows), np.asarray(profile_times), np.asarray(profiles),
            steps, cutoff_reached, t)


def make_animation(out, amplitude, acrit, config, grid, data, profile_times,
                   profiles, t_end, before, after):
    center_abs = np.abs(data[:, 1])
    peak_index = int(np.argmax(center_abs))
    peak_time = float(data[peak_index, 0])
    start = max(0.0, peak_time-before)
    end = min(t_end, peak_time+after)
    frames = np.flatnonzero((profile_times >= start) & (profile_times <= end))
    if not len(frames):
        frames = np.array([int(np.argmin(np.abs(profile_times-peak_time)))])
    start, end = float(profile_times[frames[0]]), float(profile_times[frames[-1]])

    core = grid.r <= config.core_radius
    radii = grid.r[core]
    shown_profiles = profiles[frames, :][:, core]
    y_limit = max(1.25, 1.12*float(np.max(np.abs(shown_profiles))))
    near = (data[:, 0] >= start) & (data[:, 0] <= end)
    time, center = data[near, 0], data[near, 1]
    rms = data[near, 6]

    fig = plt.figure(figsize=(8.5, 6), dpi=65)
    layout = fig.add_gridspec(2, 2, height_ratios=(2.1, 1))
    ax_profile = fig.add_subplot(layout[0, :])
    ax_center = fig.add_subplot(layout[1, 0])
    ax_rms = fig.add_subplot(layout[1, 1])

    profile, = ax_profile.plot(radii, shown_profiles[0], color="#1d4ed8", lw=1.8)
    ax_profile.axhline(0, color="black", lw=0.7, alpha=0.5)
    ax_profile.axhline(config.collapse_level, color="#dc2626", ls=":", lw=0.9)
    ax_profile.axhline(-config.collapse_level, color="#dc2626", ls=":", lw=0.9)
    ax_profile.set(
        xlim=(0, config.core_radius), ylim=(-y_limit, y_limit),
        xlabel="radius $r$", ylabel=r"field $\phi(r,t)$",
        title="Radial field profile inside the measured core",
    )
    stamp = ax_profile.text(
        0.02, 0.95, "", transform=ax_profile.transAxes, va="top",
        bbox={"boxstyle": "round", "fc": "white", "ec": "#cbd5e1", "alpha": 0.9},
    )

    ax_center.plot(time, center, color="#0f766e", lw=1.2)
    center_cursor = ax_center.axvline(start, color="#dc2626", ls="--")
    ax_center.axhline(0, color="black", lw=0.6, alpha=0.4)
    ax_center.set(xlim=(start, end), xlabel="time $t$",
                  ylabel=r"$\phi(0,t)$", title="Central field")

    ax_rms.plot(time, rms, color="#7c3aed", lw=1.2)
    rms_cursor = ax_rms.axvline(start, color="#dc2626", ls="--")
    ax_rms.set(
        xlim=(start, end), xlabel="time $t$",
        ylabel=r"$\sqrt{\langle r^2\rangle}$",
        title=rf"Signed-energy RMS radius ($R_{{core}}={config.core_radius:g}$)",
    )
    for ax in (ax_center, ax_rms):
        ax.grid(alpha=0.25)
    fig.suptitle(
        rf"Near-critical subcritical run: $A_0={amplitude:.12f}$ "
        + rf"($A_{{crit}}\approx{acrit:.12f}$)",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    def update(frame_number):
        profile_index = frames[frame_number]
        frame_time = profile_times[profile_index]
        data_index = int(np.argmin(np.abs(data[:, 0]-frame_time)))
        profile.set_ydata(profiles[profile_index, core])
        center_cursor.set_xdata((frame_time, frame_time))
        rms_cursor.set_xdata((frame_time, frame_time))
        stamp.set_text(
            f"t={frame_time:.3f}; phi(0,t)={data[data_index, 1]:.3f}; "
            f"RMS radius={data[data_index, 6]:.3f}"
        )
        return profile, center_cursor, rms_cursor, stamp

    animation = FuncAnimation(fig, update, frames=len(frames), blit=False)
    gif_path = out / "subcritical_core_evolution.gif"
    animation.save(gif_path, writer=PillowWriter(fps=15))
    plt.close(fig)
    return gif_path, peak_time, float(center_abs[peak_index]), [start, end], len(frames)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amplitude", type=float,
                        help="subcritical A0; defaults to the saved lower bracket edge")
    parser.add_argument("--t-max", type=float, default=100.0)
    parser.add_argument("--sample-interval", type=float, default=0.02)
    parser.add_argument("--profile-interval", type=float, default=0.2)
    parser.add_argument("--window-before", type=float, default=6.0)
    parser.add_argument("--window-after", type=float, default=12.0)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    critical_path = HERE / "results" / "module4_critical_summary.json"
    if not critical_path.exists():
        parser.error("Run the Module 4 critical search first; its summary is missing.")
    critical = json.loads(critical_path.read_text(encoding="utf-8"))
    acrit = float(critical["Acrit_estimate"])
    amplitude = (float(critical["Acrit_bracket"][0])
                 if args.amplitude is None else args.amplitude)
    if args.t_max <= 0 or args.sample_interval <= 0 or args.profile_interval <= 0:
        parser.error("Time horizon and sampling intervals must be positive.")
    if amplitude >= acrit:
        parser.error("Choose A0 below the saved Acrit estimate for a subcritical run.")

    config = Config(sample_interval=args.sample_interval)
    grid = make_grid(config)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Integrating subcritical A0={amplitude:.12f} to t={args.t_max:g}...",
          flush=True)
    data, profile_times, profiles, steps, cutoff, t_end = evolve(
        amplitude, config, grid, args.t_max,
        args.sample_interval, args.profile_interval,
    )
    if not len(data):
        raise RuntimeError("The run produced no diagnostic samples.")

    csv_path = args.output_dir / "subcritical_near_critical_timeseries.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        writer.writerows(data)
    print("Creating focused radial-profile animation...", flush=True)
    gif_path, peak_time, peak_center, window, frame_count = make_animation(
        args.output_dir, amplitude, acrit, config, grid, data,
        profile_times, profiles, t_end, args.window_before, args.window_after,
    )

    summary = {
        "amplitude": amplitude,
        "Acrit_estimate": acrit,
        "amplitude_offset_percent": 100*(amplitude/acrit-1),
        "classification": ("collapse cutoff reached" if cutoff
                           else f"subcritical through t={t_end:g}"),
        "collapse_cutoff": config.collapse_level,
        "maximum_abs_phi_center": peak_center,
        "time_of_maximum_abs_phi_center": peak_time,
        "animation_window": window,
        "animation_frames": frame_count,
        "integration_steps": steps,
        "timeseries_csv": csv_path.name,
        "animation_gif": gif_path.name,
        "config": {"n": config.n, "r_max": config.r_max,
                   "r0": config.r0, "core_radius": config.core_radius,
                   "sample_interval": config.sample_interval,
                   "profile_interval": args.profile_interval},
    }
    (args.output_dir / "subcritical_near_critical_summary.json").write_text(
        json.dumps(summary, indent=2)+"\n", encoding="utf-8",
    )
    print(f"A0={amplitude:.12f}; max |phi(0,t)|={peak_center:.6f} at t={peak_time:.3f}")
    print(f"{summary['classification']}; saved {frame_count} frames to {gif_path}")


if __name__ == "__main__":
    main()
