#!/usr/bin/env python3
"""Run the C++ oscillon solver and plot the field at the origin."""

from pathlib import Path
import argparse
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, FuncAnimation
import numpy as np


SIM_DIR = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--A", type=float, default=0.8,
                        help="initial amplitude")
    parser.add_argument("--R0", type=float, default=2.86,
                        help="Gaussian radius")
    parser.add_argument("--T", type=float, default=1000.0,
                        help="final simulation time")
    parser.add_argument("--snapshot-every", type=int, default=1000,
                        help="full spatial snapshot interval for the video")
    parser.add_argument("--origin-every", type=int, default=25,
                        help="origin-output interval")
    parser.add_argument("--N", type=int, default=250,
                        help="number of spatial points (use 500 for convergence)")
    parser.add_argument("--output", default="output/long_lived",
                        help="directory for solver output and the plot")
    parser.add_argument("--solver", default="oscillon",
                        help="compiled C++ solver path")
    args = parser.parse_args()

    solver = Path(args.solver)
    if not solver.is_absolute():
        solver = SIM_DIR / solver
    output_directory = Path(args.output)
    if not output_directory.is_absolute():
        output_directory = SIM_DIR / output_directory
    output_directory.mkdir(parents=True, exist_ok=True)

    command = [
        str(solver), "single", str(args.A), str(args.R0), str(args.T),
        str(args.snapshot_every), str(args.origin_every),
        str(output_directory), str(args.N),
    ]
    subprocess.run(command, cwd=SIM_DIR, check=True)

    origin_file = output_directory / "origin.dat"
    data = np.loadtxt(origin_file)
    data = np.atleast_2d(data)
    time, displacement = data[:, 0], data[:, 1]

    figure, axis = plt.subplots(figsize=(10, 4.5))
    axis.plot(time, displacement, color="navy", linewidth=0.7)
    axis.set_xlabel("time")
    axis.set_ylabel(r"origin displacement $\phi(0,t)$")
    axis.set_title("Oscillon displacement at the origin")
    axis.grid(alpha=0.25)
    figure.tight_layout()

    plot_file = output_directory / "origin_displacement.png"
    figure.savefig(plot_file, dpi=160)
    print(f"wrote {plot_file}")

    energy_file = output_directory / "energy.dat"
    energy_data = np.loadtxt(energy_file)
    energy_data = np.atleast_2d(energy_data)
    energy_time, total_energy = energy_data[:, 0], energy_data[:, 1]

    energy_figure, energy_axis = plt.subplots(figsize=(10, 4.5))
    energy_axis.plot(energy_time, total_energy, color="darkgreen", linewidth=0.9)
    energy_axis.set_xlabel("time")
    energy_axis.set_ylabel("total energy E(t)")
    energy_axis.set_title("Total field energy")
    energy_axis.grid(alpha=0.25)
    energy_figure.tight_layout()

    energy_plot_file = output_directory / "energy.png"
    energy_figure.savefig(energy_plot_file, dpi=160)
    print(f"wrote {energy_plot_file}")

    snapshots = read_snapshots(output_directory / "snapshots.dat")
    radius = snapshots[0][1]
    field_scale = max(
        max(np.max(np.abs(values)) for _, _, values in snapshots), 1.0e-8
    )
    video_figure, video_axis = plt.subplots(figsize=(8, 4.5))
    (profile_line,) = video_axis.plot(
        radius, snapshots[0][2], color="navy", linewidth=2
    )
    time_label = video_axis.text(0.02, 0.93, "", transform=video_axis.transAxes)
    video_axis.set_xlim(radius[0], radius[-1])
    video_axis.set_ylim(-1.15 * field_scale, 1.15 * field_scale)
    video_axis.set_xlabel("r")
    video_axis.set_ylabel(r"$\phi(r,t)$")
    video_axis.set_title("Oscillon evolution at accelerated playback")
    video_axis.axhline(0.0, color="gray", linewidth=0.7)

    def update(frame):
        current_time, current_radius, values = snapshots[frame]
        profile_line.set_data(current_radius, values)
        time_label.set_text(f"t = {current_time:.3f}")
        return profile_line, time_label

    animation = FuncAnimation(
        video_figure, update, frames=len(snapshots), interval=40, blit=True
    )
    video_file = output_directory / "oscillon_speedup.mp4"
    animation.save(video_file, writer=FFMpegWriter(fps=25), dpi=140)
    print(f"wrote {video_file}")


def read_snapshots(path):
    snapshots = []
    current_time = None
    radii = []
    values = []

    def save_current():
        if radii:
            snapshots.append(
                (current_time, np.asarray(radii), np.asarray(values))
            )

    with path.open() as stream:
        for line in stream:
            fields = line.split()
            if not fields:
                save_current()
                current_time = None
                radii = []
                values = []
                continue
            current_time, radius, value = map(float, fields)
            radii.append(radius)
            values.append(value)

    save_current()
    return snapshots


if __name__ == "__main__":
    main()
