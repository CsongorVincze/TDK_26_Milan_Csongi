#!/usr/bin/env python3
"""Plot the existing alapkod trace around its late-time central-field decay."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PROJECT_DIR = Path(__file__).resolve().parent
DEFAULT_SOURCE = PROJECT_DIR / "results" / "milan_alapkod_data_20260930_105152.csv"
RESULTS_DIR = PROJECT_DIR / "results"
BASELINE_WINDOW = (4000.0, 4700.0)
SEARCH_START = 4700.0
ONSET_FRACTION, END_FRACTION = 0.75, 0.10
ENVELOPE_WINDOW = 100.0


def moving_rms(signal, sample_width):
    """Smooth the oscillating central field with a centered RMS window."""
    width = max(3, int(sample_width))
    if width % 2 == 0:
        width += 1
    kernel = np.ones(width) / width
    return np.sqrt(np.convolve(signal**2, kernel, mode="same")), width


def crossing_time(time, envelope, threshold, start_at):
    hits = np.flatnonzero((time >= start_at) & (envelope <= threshold))
    if not len(hits):
        raise RuntimeError(f"Central-field RMS never fell below {threshold:.6g}.")
    return float(time[hits[0]])


def save_plot(path, time, center, mean_r2, envelope, baseline_rms,
              decay_start, decay_end):
    plot_start, plot_end = decay_start - 400.0, decay_end + 250.0
    visible = (time >= plot_start) & (time <= plot_end)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.axvspan(decay_start, decay_end, color="#f59e0b", alpha=0.11,
               label="75% to 10% decay interval")
    ax.plot(time[visible], center[visible], color="#334155", lw=0.8,
            alpha=0.8, label=r"central field $\phi(0,t)$")
    ax.plot(time[visible], envelope[visible], color="#0f766e", lw=2,
            label=f"{ENVELOPE_WINDOW:g}-time-unit central-field RMS")
    ax.axhline(ONSET_FRACTION * baseline_rms, color="#d97706", ls=":", lw=1.2,
               label="75% of pre-decay RMS")
    ax.axhline(END_FRACTION * baseline_rms, color="#7c3aed", ls=":", lw=1.2,
               label="10% of pre-decay RMS")
    ax.axvline(decay_start, color="#d97706", ls="--", lw=1.3,
               label=f"decay starts: t={decay_start:.1f}")
    ax.axvline(decay_end, color="#7c3aed", ls="--", lw=1.3,
               label=f"decay endpoint: t={decay_end:.1f}")
    ax.set(
        xlabel="time t", ylabel=r"$\phi(0,t)$ and central-field RMS",
        title="Baseline long-time decay and energy-weighted radial moment",
        xlim=(plot_start, plot_end),
    )
    ax.grid(True, alpha=0.22)

    moment_ax = ax.twinx()
    moment_ax.plot(time[visible], mean_r2[visible], color="#7c3aed", lw=1.8,
                   label=r"$\langle r^2\rangle(t)$")
    moment_ax.set_ylabel(r"$\langle r^2\rangle$ (length$^2$)", color="#7c3aed")
    moment_ax.tick_params(axis="y", colors="#7c3aed")
    handles, labels = ax.get_legend_handles_labels()
    moment_handles, moment_labels = moment_ax.get_legend_handles_labels()
    fig.legend(handles + moment_handles, labels + moment_labels,
               loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=4,
               fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=RESULTS_DIR)
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f"Input CSV not found: {args.input}")

    data = np.genfromtxt(args.input, delimiter=",", names=True)
    time, center = data["time"], data["phi_center"]
    mean_r2 = data["mean_r2"]
    dt = float(np.median(np.diff(time)))
    envelope, width = moving_rms(center, ENVELOPE_WINDOW / dt)
    baseline = ((time >= BASELINE_WINDOW[0]) & (time < BASELINE_WINDOW[1]))
    baseline_rms = float(np.median(envelope[baseline]))
    decay_start = crossing_time(
        time, envelope, ONSET_FRACTION * baseline_rms, SEARCH_START,
    )
    decay_end = crossing_time(
        time, envelope, END_FRACTION * baseline_rms, decay_start,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = args.output_dir / "long_time_decay_moments.png"
    save_plot(plot_path, time, center, mean_r2, envelope, baseline_rms,
              decay_start, decay_end)
    plot_start, plot_end = decay_start - 400.0, decay_end + 250.0
    visible = (time >= plot_start) & (time <= plot_end)
    window_path = args.output_dir / "long_time_decay_moments_window.csv"
    with window_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("time", "phi_center", "central_rms_100", "mean_r2"))
        writer.writerows(zip(time[visible], center[visible],
                             envelope[visible], mean_r2[visible]))

    nearest_5000 = int(np.argmin(np.abs(time - 5000.0)))
    moment_interval = (time >= decay_start) & (time <= plot_end)
    peak_index = np.flatnonzero(moment_interval)[
        int(np.argmax(mean_r2[moment_interval]))
    ]
    summary = {
        "source_csv": str(args.input.resolve()),
        "central_field": "phi_center = phi(t, r=0)",
        "moment_column": "mean_r2 = <r^2>(t)",
        "core_radius_cutoff": 8.0,
        "decay_marker_method": {
            "central_field_envelope": "centered moving RMS",
            "envelope_window_time_units": ENVELOPE_WINDOW,
            "baseline_time_window": list(BASELINE_WINDOW),
            "baseline_rms": baseline_rms,
            "start_fraction": ONSET_FRACTION,
            "start_time": decay_start,
            "endpoint_fraction": END_FRACTION,
            "endpoint_time": decay_end,
        },
        "near_time_5000": {
            "time": float(time[nearest_5000]),
            "phi_center": float(center[nearest_5000]),
            "mean_r2": float(mean_r2[nearest_5000]),
        },
        "largest_mean_r2_in_decay_view": {
            "time": float(time[peak_index]),
            "mean_r2": float(mean_r2[peak_index]),
        },
        "plot_file": str(plot_path.resolve()),
        "window_csv": str(window_path.resolve()),
        "beta_fit_performed": False,
    }
    summary_path = args.output_dir / "long_time_decay_moments_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(
        f"Decay start t={decay_start:.3f}, 10% endpoint t={decay_end:.3f}; "
        f"<r^2>(nearest t=5000)={mean_r2[nearest_5000]:.5f}."
    )
    print(f"Plot: {plot_path.resolve()}")


if __name__ == "__main__":
    main()
