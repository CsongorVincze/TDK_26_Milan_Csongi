#!/usr/bin/env python3
"""
Visualizes the output of the C++ radial oscillon simulation.

Reads from output/data:
  snapshots.dat   lines of "t r phi", blank line between time snapshots
  diagnostics.dat lines of "t phi(0,t) E(t)"
  run_info.txt    exact parameters and run date

Usage:
  python3 python/plot_oscillon.py                  # interactive plots
  python3 python/plot_oscillon.py --no-show        # save dated PNGs
  python3 python/plot_oscillon.py --no-show --save-mp4  # save dated MP4 too
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "output" / "data"
IMAGE_DIR = PROJECT_ROOT / "output" / "images"
VIDEO_DIR = PROJECT_ROOT / "output" / "video"

SNAPSHOT_FILE = DATA_DIR / "snapshots.dat"
DIAGNOSTIC_FILE = DATA_DIR / "diagnostics.dat"
METADATA_FILE = DATA_DIR / "run_info.txt"


def read_metadata(path):
    """Returns the exact metadata lines and the run date."""
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]
    date = next((line.split("=", 1)[1] for line in lines
                 if line.startswith("date=")), "unknown")
    return lines, date


def read_snapshots(path):
    """Parses the snapshot file into a list of (t, r_array, phi_array)."""
    snaps = []
    t_cur, rs, phis = None, [], []

    def flush():
        if rs and len(rs) == len(phis):
            snaps.append((t_cur, np.array(rs), np.array(phis)))

    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                flush()
                t_cur, rs, phis = None, [], []
                continue
            t, r, phi = map(float, parts)
            if t_cur is not None and not np.isclose(t, t_cur):
                flush()
                rs, phis = [], []
            t_cur = t
            rs.append(r)
            phis.append(phi)
    flush()
    return snaps


def animate(snaps, metadata_lines, run_date, save_mp4=False):
    """Animated evolution of phi(r) in time."""
    if not snaps:
        raise ValueError("no snapshots found")

    r_max = max(s[1].max() for s in snaps)
    amplitude = max(max(np.max(np.abs(s[2])) for s in snaps), 1e-8)

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    fig.subplots_adjust(bottom=0.42)
    (line,) = ax.plot([], [], lw=2, color="crimson")
    label = ax.text(0.02, 0.95, "", transform=ax.transAxes)
    fig.text(0.01, 0.01, "\n".join(metadata_lines),
             fontsize=5.5, family="monospace", va="bottom")

    ax.set_xlim(0, r_max)
    ax.set_ylim(-1.15 * amplitude, 1.15 * amplitude)
    ax.set_xlabel("r")
    ax.set_ylabel(r"$\phi(r)$")
    ax.set_title("Radial scalar field profile")
    ax.axhline(0.0, color="gray", lw=0.5)

    def update(k):
        t, r, phi = snaps[k]
        line.set_data(r, phi)
        label.set_text(f"t = {t:.2f}")
        return line, label

    ani = FuncAnimation(fig, update, frames=len(snaps), interval=40, blit=True)

    if save_mp4:
        video_path = VIDEO_DIR / f"oscillon_{run_date}.mp4"
        try:
            ani.save(video_path, fps=25, dpi=150)
            print(f"wrote {video_path}")
        except Exception as e:
            gif_path = VIDEO_DIR / f"oscillon_{run_date}.gif"
            print(f"mp4 export failed ({e}); falling back to {gif_path}")
            ani.save(gif_path, writer=PillowWriter(fps=25))
            print(f"wrote {gif_path}")
    return fig, ani


def plot_diagnostics(metadata_lines, run_date, show=True):
    """Center amplitude phi(0,t) and total energy E(t)."""
    d = np.atleast_2d(np.loadtxt(DIAGNOSTIC_FILE))
    t, phi0, energy = d[:, 0], d[:, 1], d[:, 2]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 5.5), sharex=True)
    ax1.plot(t, phi0, color="navy", lw=1.2)
    ax1.set_ylabel(r"$\phi(0,t)$")
    ax1.set_title("Central field amplitude")

    ax2.plot(t, energy, color="darkgreen", lw=1.2)
    ax2.set_xlabel("t")
    ax2.set_ylabel("E(t)")
    ax2.set_title("Total energy (decreases in the sponge layer)")

    fig.tight_layout(rect=(0, 0.28, 1, 1))
    fig.text(0.01, 0.01, "\n".join(metadata_lines),
             fontsize=5.5, family="monospace", va="bottom")
    image_path = IMAGE_DIR / f"diagnostics_{run_date}.png"
    fig.savefig(image_path, dpi=150)
    print(f"wrote {image_path}")
    if show:
        plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-show", action="store_true", help="headless mode, only save pngs")
    parser.add_argument("--save-mp4", action="store_true", help="render the animation to a file")
    args = parser.parse_args()

    IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    VIDEO_DIR.mkdir(parents=True, exist_ok=True)
    metadata_lines, run_date = read_metadata(METADATA_FILE)
    snapshots = read_snapshots(SNAPSHOT_FILE)
    print(f"loaded {len(snapshots)} snapshots")

    fig_anim, _ = animate(snapshots, metadata_lines, run_date, save_mp4=args.save_mp4)
    if args.no_show:
        image_path = IMAGE_DIR / f"oscillon_evolution_{run_date}.png"
        fig_anim.savefig(image_path, dpi=150)
        print(f"wrote {image_path}")

    plot_diagnostics(metadata_lines, run_date, show=not args.no_show)

    if not args.no_show:
        plt.show()   # opens both windows
