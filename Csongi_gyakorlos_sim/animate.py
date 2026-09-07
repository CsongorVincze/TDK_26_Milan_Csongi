#!/usr/bin/env python3
"""Create a basic animation from output/snapshots.dat."""

from pathlib import Path
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
import numpy as np


ROOT = Path(__file__).resolve().parent


def read_snapshots(path):
    snapshots = []
    current_t = None
    radii = []
    values = []

    def save_current():
        if radii:
            snapshots.append(
                (current_t, np.asarray(radii), np.asarray(values))
            )

    with path.open() as stream:
        for line in stream:
            fields = line.split()
            if not fields:
                save_current()
                current_t = None
                radii = []
                values = []
                continue

            t, radius, phi = map(float, fields)
            current_t = t
            radii.append(radius)
            values.append(phi)

    save_current()
    return snapshots


def read_metadata(path):
    if not path.exists():
        return [], "unknown"

    lines = [line.strip() for line in path.read_text().splitlines()
             if line.strip()]
    date = next((line.split("=", 1)[1] for line in lines
                 if line.startswith("date=")), "unknown")
    return lines, date


def resolve_run_directory(value):
    directory = Path(value)
    if not directory.is_absolute():
        directory = ROOT / directory
    return directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run",
        default="output",
        help="run directory, e.g. output/sweep/run_001_A0p40_R02p00",
    )
    args = parser.parse_args()

    run_directory = resolve_run_directory(args.run)
    snapshot_file = run_directory / "snapshots.dat"
    metadata_file = run_directory / "run_info.txt"
    if not snapshot_file.exists():
        raise SystemExit(
            f"{snapshot_file} does not exist; run the C++ simulation first"
        )

    snapshots = read_snapshots(snapshot_file)
    if not snapshots:
        raise SystemExit("the snapshot file is empty")
    metadata_lines, run_date = read_metadata(metadata_file)

    radius = snapshots[0][1]
    amplitude = max(
        max(np.max(np.abs(values)) for _, _, values in snapshots), 1.0e-8
    )

    fig, ax = plt.subplots(figsize=(8, 4.5))
    fig.subplots_adjust(bottom=0.34)
    (line,) = ax.plot(radius, snapshots[0][2], color="navy", linewidth=2)
    time_label = ax.text(0.02, 0.93, "", transform=ax.transAxes)
    if metadata_lines:
        fig.text(0.01, 0.01, "\n".join(metadata_lines),
                 fontsize=5.5, family="monospace", va="bottom")

    ax.set_xlim(radius[0], radius[-1])
    ax.set_ylim(-1.15 * amplitude, 1.15 * amplitude)
    ax.set_xlabel("r")
    ax.set_ylabel(r"$\phi(r,t)$")
    ax.set_title("Radial scalar-field evolution")
    ax.axhline(0.0, color="gray", linewidth=0.7)

    def update(frame):
        t, r, values = snapshots[frame]
        line.set_data(r, values)
        time_label.set_text(f"t = {t:.3f}")
        return line, time_label

    animation = FuncAnimation(
        fig, update, frames=len(snapshots), interval=40, blit=True
    )

    video_dir = run_directory / "video"
    video_dir.mkdir(parents=True, exist_ok=True)
    mp4_path = video_dir / f"oscillon_{run_date}.mp4"
    try:
        animation.save(mp4_path, writer=FFMpegWriter(fps=25), dpi=150)
        print(f"wrote {mp4_path}")
    except Exception as error:
        gif_path = video_dir / f"oscillon_{run_date}.gif"
        print(f"MP4 export failed ({error}); writing {gif_path}")
        animation.save(gif_path, writer=PillowWriter(fps=25))
        print(f"wrote {gif_path}")


if __name__ == "__main__":
    main()
