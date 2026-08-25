#!/usr/bin/env python3
"""Create a basic animation from output/snapshots.dat."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
import numpy as np


ROOT = Path(__file__).resolve().parent
SNAPSHOT_FILE = ROOT / "output" / "snapshots.dat"
VIDEO_DIR = ROOT / "output" / "video"


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


def main():
    if not SNAPSHOT_FILE.exists():
        raise SystemExit(
            f"{SNAPSHOT_FILE} does not exist; run the C++ simulation first"
        )

    snapshots = read_snapshots(SNAPSHOT_FILE)
    if not snapshots:
        raise SystemExit("the snapshot file is empty")

    radius = snapshots[0][1]
    amplitude = max(
        max(np.max(np.abs(values)) for _, _, values in snapshots), 1.0e-8
    )

    fig, ax = plt.subplots(figsize=(8, 4.5))
    (line,) = ax.plot(radius, snapshots[0][2], color="navy", linewidth=2)
    time_label = ax.text(0.02, 0.93, "", transform=ax.transAxes)

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

    VIDEO_DIR.mkdir(parents=True, exist_ok=True)
    mp4_path = VIDEO_DIR / "oscillon.mp4"
    try:
        animation.save(mp4_path, writer=FFMpegWriter(fps=25), dpi=150)
        print(f"wrote {mp4_path}")
    except Exception as error:
        gif_path = VIDEO_DIR / "oscillon.gif"
        print(f"MP4 export failed ({error}); writing {gif_path}")
        animation.save(gif_path, writer=PillowWriter(fps=25))
        print(f"wrote {gif_path}")


if __name__ == "__main__":
    main()
