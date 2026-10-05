#!/usr/bin/env python3
"""Animate four subcritical, one threshold, and four supercritical runs."""

import csv
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig_module4_fun")

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

MODULE4 = HERE.parent
sys.path.insert(0, str(MODULE4))
from milan_alapkod_module4 import Config, make_grid, rk4  # noqa: E402


# N=500 operational threshold from the Module 4 threshold search.
A_CRIT = 1.202609700895846
COLLAPSE_LEVEL = 12.0
DISPLAY_T_MAX = 100.0
FRAME_COUNT = 201
CORE_RADIUS = 8.0

COLORS = {
    "subcritical": "#2878b5",
    "critical": "#d69e00",
    "supercritical": "#d1495b",
}

CASE_SPECS = (
    ("subcritical", -5.0, 0, 0),
    ("subcritical", -0.5, 0, 1),
    ("subcritical", -0.05, 0, 2),
    ("subcritical", -0.01, 1, 0),
    ("critical", 0.0, 1, 1),
    ("supercritical", 0.01, 1, 2),
    ("supercritical", 0.05, 2, 0),
    ("supercritical", 0.5, 2, 1),
    ("supercritical", 5.0, 2, 2),
)


def evolve_profiles(amplitude, config, grid, frame_times, plot_mask):
    """Return sampled radial profiles and the first |phi(0)|=12 crossing."""
    phi = amplitude * np.exp(-grid.r**2 / (2 * config.r0**2))
    pi = np.zeros_like(phi)
    phi[-1] = 0.0
    profiles = np.empty((len(frame_times), int(np.count_nonzero(plot_mask))),
                        dtype=np.float32)
    profiles[0] = phi[plot_mask]
    frame_index = 1
    t = 0.0
    collapse_time = None
    failed = False

    while frame_index < len(frame_times):
        target_time = frame_times[frame_index]
        peak = max(float(np.max(np.abs(phi))), 1.0)
        dt = min(grid.dt, 0.15 / peak, target_time - t)
        new_phi, new_pi = rk4(phi, pi, dt, grid)

        if not np.all(np.isfinite(new_phi)) or not np.all(np.isfinite(new_pi)):
            failed = True
            profiles[frame_index:] = phi[plot_mask]
            break

        if abs(new_phi[0]) >= config.collapse_level:
            rise = abs(new_phi[0]) - abs(phi[0])
            fraction = np.clip(
                (config.collapse_level - abs(phi[0])) / rise
                if rise > 0 else 1.0,
                0.0, 1.0,
            )
            t += dt * fraction
            phi = phi + fraction * (new_phi - phi)
            pi = pi + fraction * (new_pi - pi)
            phi[0] = np.copysign(config.collapse_level, new_phi[0])
            collapse_time = t
            profiles[frame_index:] = phi[plot_mask]
            break

        t += dt
        phi, pi = new_phi, new_pi
        if t >= target_time - 1e-12:
            profiles[frame_index] = phi[plot_mask]
            frame_index += 1

    return profiles, collapse_time, failed


def make_cases(config, grid, frame_times):
    plot_mask = grid.r <= CORE_RADIUS
    cases = []
    for category, offset_percent, row, col in CASE_SPECS:
        amplitude = A_CRIT * (1.0 + offset_percent / 100.0)
        profiles, collapse_time, failed = evolve_profiles(
            amplitude, config, grid, frame_times, plot_mask,
        )
        if category == "critical":
            label = "ESTIMATED CRITICAL · 0%"
        else:
            label = f"{category.upper()} · {offset_percent:+g}%"
        cases.append({
            "category": category,
            "offset_percent": offset_percent,
            "amplitude": amplitude,
            "row": row,
            "col": col,
            "label": label,
            "profiles": profiles,
            "collapse_time": collapse_time,
            "failed": failed,
        })
        end_state = (f"cutoff at t={collapse_time:.3f}"
                     if collapse_time is not None else "no cutoff")
        if failed:
            end_state = "numerical failure"
        print(f"{label}: A={amplitude:.10f}, {end_state}", flush=True)
    return cases, plot_mask


def save_summary(cases, path):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("category", "relative_offset_percent", "amplitude",
                         "collapse_cutoff_time", "numerical_failure"))
        for case in cases:
            writer.writerow((
                case["category"], case["offset_percent"], case["amplitude"],
                "" if case["collapse_time"] is None else case["collapse_time"],
                case["failed"],
            ))


def render_animation(cases, radii, frame_times, output_path):
    fig, axes = plt.subplots(3, 3, figsize=(12.5, 10.0), sharex=True, sharey=True)
    lines = {}
    center_texts = {}
    status_texts = {}

    for case in cases:
        ax = axes[case["row"], case["col"]]
        color = COLORS[case["category"]]
        initial_profile = case["profiles"][0]
        ax.set_facecolor({
            "subcritical": "#f1f7fb",
            "critical": "#fff8df",
            "supercritical": "#fcf1f2",
        }[case["category"]])
        (line,) = ax.plot(radii, initial_profile, color=color, lw=1.6)
        (center_marker,) = ax.plot([0], [initial_profile[0]], "o",
                                   color=color, ms=4)
        lines[(case["row"], case["col"])] = (line, center_marker)
        ax.set_title(case["label"], loc="left", fontsize=9, color=color,
                     fontweight="bold")
        ax.set_xlim(0, CORE_RADIUS)
        ax.set_ylim(-12.8, 12.8)
        ax.set_xticks((0, 2, 4, 6, 8))
        ax.set_yticks((-10, -5, 0, 5, 10))
        ax.grid(alpha=0.2, lw=0.6)
        center_texts[(case["row"], case["col"])] = ax.text(
            0.97, 0.92, f"φ(0,t) = {initial_profile[0]:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
        )
        status_texts[(case["row"], case["col"])] = ax.text(
            0.03, 0.06, "", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8, color=color,
        )

    for row in range(3):
        axes[row, 0].set_ylabel("field  φ(r,t)")
    for col in range(3):
        axes[2, col].set_xlabel("radius  r")

    fig.suptitle("The collapse threshold: 4 below · 1 estimated critical · 4 above",
                 fontsize=15, fontweight="bold", y=0.985)
    clock = fig.text(0.5, 0.945, "t = 0.00", ha="center", fontsize=11)
    fig.text(
        0.5, 0.035,
        rf"Shared axes · $A_{{crit}}={A_CRIT:.10f}$ · stop marker: $|\phi(0,t)|={COLLAPSE_LEVEL:g}$",
        ha="center", fontsize=9,
    )
    legend = [
        Line2D([], [], color=COLORS["subcritical"], lw=2, label="subcritical"),
        Line2D([], [], color=COLORS["critical"], lw=2, label="estimated threshold"),
        Line2D([], [], color=COLORS["supercritical"], lw=2, label="supercritical"),
    ]
    fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.065),
               ncol=3, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0.035, 0.10, 0.98, 0.925))

    def update(frame_index):
        now = frame_times[frame_index]
        clock.set_text(f"t = {now:.2f}")
        for case in cases:
            key = (case["row"], case["col"])
            profile = case["profiles"][frame_index]
            line, marker = lines[key]
            line.set_ydata(profile)
            marker.set_ydata([profile[0]])
            center_texts[key].set_text(f"φ(0,t) = {profile[0]:.2f}")

            collapse_time = case["collapse_time"]
            if case["failed"]:
                status, status_color = "numerical issue", "#8b1e2d"
            elif collapse_time is None:
                status, status_color = f"no cutoff by t={DISPLAY_T_MAX:g}", COLORS["subcritical"]
            elif now >= collapse_time:
                status = f"cutoff reached · t={collapse_time:.2f}"
                status_color = COLORS["supercritical"]
            else:
                status = f"cutoff in {collapse_time - now:.1f} time units"
                status_color = COLORS[case["category"]]
            status_texts[key].set_text(status)
            status_texts[key].set_color(status_color)
        return []

    animation = FuncAnimation(fig, update, frames=len(frame_times), interval=75,
                              repeat=True, blit=False)
    animation.save(output_path, writer=PillowWriter(fps=13), dpi=85)
    plt.close(fig)


def main():
    config = Config(n=500, collapse_level=COLLAPSE_LEVEL)
    grid = make_grid(config)
    frame_times = np.linspace(0.0, DISPLAY_T_MAX, FRAME_COUNT)
    cases, plot_mask = make_cases(config, grid, frame_times)
    save_summary(cases, RESULTS / "case_summary.csv")
    render_animation(cases, grid.r[plot_mask], frame_times,
                     RESULTS / "critical_transition_3x3.gif")
    print(f"Saved animation and case summary under {RESULTS}")


if __name__ == "__main__":
    main()
