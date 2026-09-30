"""Evolve the radial double-well field and save its diagnostics."""

import csv
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as itg
from matplotlib.lines import Line2D


# User-adjustable simulation parameters
N, R_max, alpha = 500, 25, 5
d = 3
A, R_0 = 2.0, 1.915
C_CFL = 0.4
R_sponge_fraction, gamma_0, p = 0.8, 7, 3
t0, T, time_samples = 0, 7000, 10000
solver_method, solver_rtol, solver_atol = "RK45", 1e-6, 1e-9

# Output and plotting options
show_progress = True
progress_interval_seconds = 5.0
extra_plots, save_figures, show_plots = True, True, True
results_directory = Path(__file__).resolve().parent / "results"
extra_plot_directory = results_directory / "extra_plots"
potential_plot_min, potential_plot_max, potential_plot_samples = -2, 2, 100
potential_plot_file = "potential_functions.png"
grid_plot_file = "grid_and_sponge.png"
initial_profile_file = "initial_profile.png"
plot_profile_count, plot_amplitude_factor = 5, 5
output_file = results_directory / "test.jpg"
energy_plot_file = results_directory / "energy_time.png"
core_radius = 8.0
csv_output_file = results_directory / (
    f"milan_alapkod_data_{datetime.now():%Y%m%d_%H%M%S}.csv"
)
color_map, plot_dpi = "winter", 150
end = -1

R_sponge = R_sponge_fraction * R_max
date_tag = datetime.now().strftime("%m-%d")


@dataclass
class Grid:
    """Stretched radial grid and coefficients reused by the PDE solver."""

    x: np.ndarray
    r: np.ndarray
    dr_dx: np.ndarray
    d2r_dx2: np.ndarray
    gamma: np.ndarray
    dx: float
    dr_min: float
    inv_dx2: float
    inv_2dx: float
    interior_inv_dr_dx2: np.ndarray
    interior_radial_coefficient: np.ndarray
    origin_laplacian_coefficient: float
    outer_inv_dr_dx2: float
    outer_radial_coefficient: float


def make_grid():
    """Build r(x), its derivatives, the sponge, and fixed FD coefficients."""
    dx = 1 / (N - 1)
    x = np.linspace(0, 1, N)
    sinh_alpha = np.sinh(alpha)
    r = R_max * np.sinh(alpha * x) / sinh_alpha
    dr_dx = R_max * alpha * np.cosh(alpha * x) / sinh_alpha
    d2r_dx2 = R_max * alpha**2 * np.sinh(alpha * x) / sinh_alpha

    gamma = np.zeros_like(x)
    in_sponge = r > R_sponge
    gamma[in_sponge] = gamma_0 * (
        (r[in_sponge] - R_sponge) / (R_max - R_sponge)
    ) ** p

    return Grid(
        x=x, r=r, dr_dx=dr_dx, d2r_dx2=d2r_dx2, gamma=gamma, dx=dx,
        dr_min=alpha * R_max / sinh_alpha * dx,
        inv_dx2=1 / dx**2,
        inv_2dx=1 / (2 * dx),
        interior_inv_dr_dx2=1 / dr_dx[1:-1]**2,
        interior_radial_coefficient=(
            (d - 1) / (r[1:-1] * dr_dx[1:-1])
            - d2r_dx2[1:-1] / dr_dx[1:-1]**3
        ),
        origin_laplacian_coefficient=2 * d / dr_dx[0]**2,
        outer_inv_dr_dx2=1 / dr_dx[-1]**2,
        outer_radial_coefficient=(
            (d - 1) / (r[-1] * dr_dx[-1])
            - d2r_dx2[-1] / dr_dx[-1]**3
        ),
    )


def potential(phi):
    """Shifted double-well potential, with the outer vacuum at phi=0."""
    return 0.25 * phi**2 * (phi - 2)**2


def potential_derivative(phi):
    return phi * (phi - 1) * (phi - 2)


class Progress:
    """Throttled solver progress display with an approximate time remaining."""

    def __init__(self):
        self.started_at = time.perf_counter()
        self.last_update = 0.0

    def update(self, current_t, force=False):
        if not show_progress:
            return
        now = time.perf_counter()
        if not force and now - self.last_update < progress_interval_seconds:
            return

        fraction = np.clip((current_t - t0) / max(T - t0, 1e-15), 0, 1)
        elapsed = now - self.started_at
        bar = "=" * int(30 * fraction) + " " * (30 - int(30 * fraction))
        eta = f"ETA {elapsed * (1 - fraction) / fraction / 60:.1f} min" if fraction else "ETA --"
        print(
            f"\r[{bar}] {100 * fraction:6.2f}% t={current_t:.1f}/{T:g} "
            f"| elapsed {elapsed / 60:.1f} min | {eta}",
            end="", flush=True,
        )
        self.last_update = now


def system(t, state, grid, progress):
    """First-order PDE system: phi_t=Pi and damped radial wave evolution."""
    progress.update(t)
    phi, Pi = state[:N], state[N:]

    # Enforce the outer Dirichlet values at each RHS evaluation.
    phi[-1] = Pi[-1] = 0
    derivative = np.empty_like(state)
    dphi, dPi = derivative[:N], derivative[N:]
    dphi[:] = Pi
    dPi.fill(0.0)

    # Regularity at r=0, transformed-coordinate radial Laplacian in the bulk,
    # and the one-sided outer-boundary stencil, respectively.
    dPi[0] = (
        grid.origin_laplacian_coefficient * (phi[1] - phi[0]) * grid.inv_dx2
        - potential_derivative(phi[0])
    )
    dPi[1:-1] = (
        grid.interior_inv_dr_dx2 * (phi[2:] - 2 * phi[1:-1] + phi[:-2])
        * grid.inv_dx2
        + grid.interior_radial_coefficient * (phi[2:] - phi[:-2])
        * grid.inv_2dx
        - potential_derivative(phi[1:-1])
        - grid.gamma[1:-1] * Pi[1:-1]
    )
    dPi[-1] = (
        grid.outer_inv_dr_dx2 * (0 - 2 * phi[-1] + phi[-2]) * grid.inv_dx2
        + grid.outer_radial_coefficient * (0 - phi[-1]) * grid.inv_2dx
        - potential_derivative(phi[-1]) - grid.gamma[-1] * Pi[-1]
    )
    return derivative


def initial_profile(grid):
    """Gaussian initial field and zero initial velocity."""
    phi = A * np.exp(-grid.r**2 / (2 * R_0**2))
    return phi, np.zeros_like(phi)


def evolve(grid, progress):
    """Integrate from t0 to T with the configured CFL-limited step size."""
    phi, Pi = initial_profile(grid)
    initial_state = np.concatenate((phi, Pi))
    dt = C_CFL * grid.dr_min
    times = np.linspace(t0, T, time_samples)
    solution = itg.solve_ivp(
        lambda t, state: system(t, state, grid, progress),
        (t0, T), initial_state, method=solver_method, max_step=dt,
        t_eval=times, rtol=solver_rtol, atol=solver_atol,
    )
    progress.update(solution.t[-1], force=True)
    if show_progress:
        elapsed = time.perf_counter() - progress.started_at
        print(f"\nSolver finished in {elapsed / 60:.1f} min.")
    return solution, phi, Pi


def energy_per_point(phi, Pi, grid):
    """Local energy density, including the radial derivative of phi."""
    dphi_dr = np.empty_like(phi)
    dphi_dr[0] = 0.0
    dphi_dr[1:-1] = (
        (phi[2:] - phi[:-2]) * grid.inv_2dx / grid.dr_dx[1:-1]
    )
    dphi_dr[-1] = (
        (3 * phi[-1] - 4 * phi[-2] + phi[-3])
        * grid.inv_2dx / grid.dr_dx[-1]
    )
    return 0.5 * Pi**2 + 0.5 * dphi_dr**2 + potential(phi)


def calculate_diagnostics(solution, grid):
    """Compute total/core energy and the signed-energy-weighted core radius."""
    count = solution.t.size
    energy = np.empty(count)
    core_energy = np.empty(count)
    mean_r2 = np.full(count, np.nan)
    core_mask = grid.r <= core_radius
    core_r = grid.r[core_mask]

    for j in range(count):
        phi, Pi = solution.y[:N, j], solution.y[N:, j]
        density = energy_per_point(phi, Pi, grid)
        energy[j] = 4 * np.pi * itg.simpson(
            density * grid.r**(d - 1) * grid.dr_dx, x=grid.x
        )
        weighted = density[core_mask] * core_r**(d - 1)
        denominator = itg.simpson(weighted, x=core_r)
        numerator = itg.simpson(weighted * core_r**2, x=core_r)
        scale = itg.simpson(np.abs(weighted), x=core_r)
        core_energy[j] = 4 * np.pi * denominator
        if abs(denominator) > 1e-12 * max(scale, 1.0):
            candidate = numerator / denominator
            if np.isfinite(candidate) and candidate > 0:
                mean_r2[j] = candidate
    return energy, core_energy, mean_r2


def write_csv(solution, diagnostics):
    """Write one row per stored solver time for later analysis/fitting."""
    energy, core_energy, mean_r2 = diagnostics
    path = Path(csv_output_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow((
            "time", "phi_center", "pi_center", "total_energy", "core_energy",
            "mean_r2", "core_rms_radius",
        ))
        for j, current_t in enumerate(solution.t):
            r2 = mean_r2[j]
            writer.writerow((
                current_t, solution.y[0, j], solution.y[N, j], energy[j],
                core_energy[j], r2, np.sqrt(r2) if np.isfinite(r2) else np.nan,
            ))
    print(f"Wrote time-series data to {path.resolve()}")


def dated_filename(filename):
    filename = Path(filename)
    return filename.with_name(f"{filename.stem}_{date_tag}{filename.suffix}")


def save_figure(figure, filename):
    if save_figures:
        filename = dated_filename(filename)
        filename.parent.mkdir(parents=True, exist_ok=True)
        figure.tight_layout()
        figure.savefig(filename, dpi=plot_dpi, bbox_inches="tight")


def plot_extra_diagnostics(grid, phi0, Pi0):
    """Plot the potential, stretched grid/sponge, and initial field."""
    values = np.linspace(potential_plot_min, potential_plot_max,
                         potential_plot_samples)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(values, potential(values), label=r"$V(\phi)$")
    ax.plot(values, potential_derivative(values), label=r"$V'(\phi)$")
    ax.set(title=r"Potential functions $V(\phi)$ and $V'(\phi)$",
           xlabel=r"$\phi$", ylabel="Values")
    ax.legend()
    save_figure(fig, Path(extra_plot_directory) / potential_plot_file)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].plot(grid.x, grid.r, label=r"$r(x)$")
    axes[0].plot(grid.x, grid.dr_dx, label=r"$r'(x)$")
    axes[0].plot(grid.x, grid.d2r_dx2, label=r"$r''(x)$")
    axes[0].set(title="Grid mapping", xlabel="x", ylabel="Values")
    axes[0].legend()
    axes[1].plot(grid.x, grid.gamma, label=r"$\gamma(x)$")
    axes[1].set(title="Sponge-layer damping", xlabel="x", ylabel=r"$\gamma(x)$")
    axes[1].legend()
    save_figure(fig, Path(extra_plot_directory) / grid_plot_file)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(grid.r, phi0, label=r"$\phi(0,r)$")
    ax.plot(grid.r, Pi0, label=r"$\dot{\phi}(0,r)$")
    ax.set(title="Initial profile", xlabel="r", ylabel="Values")
    ax.legend()
    save_figure(fig, Path(extra_plot_directory) / initial_profile_file)


def plot_results(solution, grid, energy):
    """Plot the central field, selected profiles, and total energy history."""
    cmap = plt.get_cmap(color_map)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].plot(solution.t[:end], solution.y[0, :end], label=r"$\phi(t,0)$")
    axes[0].set(title=r"Time evolution of $\phi_0$", xlabel="t",
                ylabel=r"$\phi(t,0)$", xlim=(t0, T),
                ylim=(-plot_amplitude_factor * A, plot_amplitude_factor * A))
    metadata = Line2D(
        [], [], color="none", linestyle="none",
        label=(rf"$A={A}, R_0={R_0}$" + "\nshifted double-well potential\n"
               + rf"$N={N}, R_{{max}}={R_max}, \alpha={alpha}$" + "\n"
               + rf"$R_{{sponge}}={R_sponge}, \gamma_0={gamma_0}, p={p}$"),
    )
    axes[0].legend(handles=[metadata], handlelength=0, handletextpad=0)

    times = solution.t[:end]
    step = max(1, int(len(times) / plot_profile_count))
    for i in range(0, len(times), step):
        axes[1].plot(
            grid.r, solution.y[:N, i], color=cmap(1 - i / len(times)),
            linewidth=1, label=f"$t={solution.t[i]:.2f}$",
        )
    axes[1].plot(grid.r, solution.y[:N, -1], color=cmap(0),
                 label=f"$t={solution.t[-1]:.2f}$")
    axes[1].set(title=r"Spatial profile $\phi(t,r)$", xlabel="r",
                ylabel=r"$\phi(t,r)$", xlim=(0, R_max),
                ylim=(-plot_amplitude_factor * A, plot_amplitude_factor * A))
    axes[1].legend(loc="upper right", fontsize=8)
    save_figure(fig, output_file)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(solution.t, energy, color="tab:purple")
    ax.set(title="Total energy of the system", xlabel="t", ylabel="E(t)")
    ax.grid(True, alpha=0.3)
    save_figure(fig, energy_plot_file)


def main():
    """Run setup → evolution → diagnostics/CSV → figures."""
    grid = make_grid()
    progress = Progress()
    solution, phi0, Pi0 = evolve(grid, progress)
    diagnostics = calculate_diagnostics(solution, grid)
    write_csv(solution, diagnostics)
    if extra_plots:
        plot_extra_diagnostics(grid, phi0, Pi0)
    plot_results(solution, grid, diagnostics[0])
    if show_plots:
        plt.show()


if __name__ == "__main__":
    main()
