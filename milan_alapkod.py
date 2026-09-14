from datetime import datetime
from pathlib import Path

import numpy as np
import scipy.integrate as itg
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.lines import Line2D

# ============================================================================
# User-adjustable parameters
# ============================================================================

# Diagnostics and plotting
extra_plots = True
save_figures = True
show_plots = True
extra_plot_directory = "extra_plots"
potential_plot_min = -2
potential_plot_max = 2
potential_plot_samples = 100
potential_plot_file = "potential_functions.png"
grid_plot_file = "grid_and_sponge.png"
initial_profile_file = "initial_profile.png"
plot_profile_count = 5
plot_amplitude_factor = 5
output_file = "test.jpg"
color_map = "winter"
plot_dpi = 150

# Grid parameters
N = 100
R_max = 10
alpha = 5

# Sponge-layer parameters
R_sponge_fraction = 0.75
gamma_0 = 7
p = 3

# Potential parameters
m = 1
mu = 1

# Courant-Friedrichs-Lewy stability condition
C_CFL = 0.4

# Spatial dimension
d = 3

# Initial profile parameters
A = 0.2
R_0 = 2.5

# Time-integration parameters
t0 = 0
T = 5
time_samples = 100
solver_method = "RK45"
solver_rtol = 1e-6
solver_atol = 1e-9

# Display parameters
end = -1

# ============================================================================
# Derived quantities
# ============================================================================

dx = 1/(N-1)
R_sponge = R_sponge_fraction * R_max
date_tag = datetime.now().strftime("%m-%d")


def dated_filename(filename):
    filename = Path(filename)
    return filename.with_name(f"{filename.stem}_{date_tag}{filename.suffix}")


def save_figure(figure, filename):
    """Save one figure when figure saving is enabled."""
    if not save_figures:
        return

    filename = dated_filename(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(filename, dpi=plot_dpi, bbox_inches="tight")

# Self-interaction potential

def V(phi):
    return 1/2*m**2*phi**2 - 1/4*mu*phi**4
def dV_dphi(phi):
    return m**2*phi - mu*phi**3

if extra_plots:
    phi_plot = np.linspace(potential_plot_min, potential_plot_max,
                           potential_plot_samples)
    Vlist = np.zeros_like(phi_plot)
    dVlist = np.zeros_like(phi_plot)
    for i, l in np.ndenumerate(phi_plot):
        Vlist[i] = V(l)
        dVlist[i] = dV_dphi(l)
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    ax.plot(phi_plot, Vlist)
    ax.plot(phi_plot, dVlist)
    ax.set_title('Potential functions $V(\\phi)$ and $V\'(\\phi)$')
    ax.set_xlabel('$\\phi$')
    ax.set_ylabel('Values')
    ax.legend(['$V(\\phi)$', '$V\'(\\phi)$'])
    save_figure(fig, Path(extra_plot_directory) / potential_plot_file)

# Non-uniform coordinate mapping

x = np.linspace(0, 1, N)
r = np.zeros_like(x)
dr_dx = np.zeros_like(x)
d2r_dx2 = np.zeros_like(x)
gamma = np.zeros_like(x)

for i in range(N):
    r[i] = R_max * (np.sinh(alpha * x[i]) / np.sinh(alpha))
    dr_dx[i] = R_max * (alpha * np.cosh(alpha * x[i]) / np.sinh(alpha))
    d2r_dx2[i] = R_max * (alpha**2 * np.sinh(alpha * x[i]) / np.sinh(alpha))
    if r[i] > R_sponge:
        gamma[i] = gamma_0 * ((r[i] - R_sponge) / (R_max - R_sponge))**p

dr_min = alpha*R_max/np.sinh(alpha)*dx

# Coefficients used by the spatial finite-difference operator. These depend
# only on the fixed grid, so calculate them once instead of in every RHS call.
inv_dx2 = 1/dx**2
inv_2dx = 1/(2*dx)
interior_inv_dr_dx2 = 1/dr_dx[1:-1]**2
interior_radial_coefficient = (
    (d-1)/(r[1:-1]*dr_dx[1:-1])
    - d2r_dx2[1:-1]/dr_dx[1:-1]**3
)
origin_laplacian_coefficient = 2*d/dr_dx[0]**2
outer_inv_dr_dx2 = 1/dr_dx[-1]**2
outer_radial_coefficient = (
    (d-1)/(r[-1]*dr_dx[-1])
    - d2r_dx2[-1]/dr_dx[-1]**3
)
interior_gamma = gamma[1:-1]

if extra_plots:
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))

    axs[0].plot(x, r)
    axs[0].plot(x, dr_dx)
    axs[0].plot(x, d2r_dx2)
    axs[0].set_title('Grid mapping functions r(x), r\'(x), r\"(x)')
    axs[0].set_xlabel('x')
    axs[0].set_ylabel('Values')
    axs[0].legend(['r(x)', 'r\'(x)', 'r\"(x)'])

    axs[1].plot(x, gamma)
    axs[1].set_title('Sponge layer damping factor $\\gamma(x)$')
    axs[1].set_xlabel('x')
    axs[1].set_ylabel('$\\gamma(x)$')
    axs[1].legend(['$\\gamma(x)$ = $\\gamma_0 \\left(\\frac{r(x) - R_{sponge}}{R_{max} - R_{sponge}}\\right)^p$'])
    save_figure(fig, Path(extra_plot_directory) / grid_plot_file)

# System of ODEs

def system(t, state):
    phi = state[:N]
    Pi = state[N:]

    #boundary conditions
    phi[-1] = 0
    Pi[-1] = 0

    # Use one derivative buffer instead of separately allocating dphi and dPi
    # and concatenating them at the end of every RHS evaluation.
    dstate = np.empty_like(state)
    dphi = dstate[:N]
    dPi = dstate[N:]
    dphi[:] = Pi
    dPi.fill(0.0)

    #? ezt itt nem teljesen ertem (csongi)
    dPi[0] = (
    origin_laplacian_coefficient * (phi[1]-phi[0]) * inv_dx2
    - dV_dphi(phi[0])
    )
    
    dPi[1:-1] = (
    interior_inv_dr_dx2 * (phi[2:]-2*phi[1:-1]+phi[:-2]) * inv_dx2
    + interior_radial_coefficient * (phi[2:]-phi[:-2]) * inv_2dx
    - dV_dphi(phi[1:-1])
    - interior_gamma*Pi[1:-1]
    )

    #? jo ezt talan ertem (csongi)
    dPi[-1] = (
    outer_inv_dr_dx2 * (0-2*phi[-1]+phi[-2]) * inv_dx2
    + outer_radial_coefficient * (0-phi[-1]) * inv_2dx
    - dV_dphi(phi[-1])
    - gamma[-1]*Pi[-1]
    )

    return dstate

# Initial profile

def phi0(A, R_0, r):
    return A * np.exp(-(r**2)/(2*R_0**2))

phi = np.zeros_like(r)
dphi_dt = np.zeros_like(r)
for i, ri in np.ndenumerate(r):
    phi[i] = phi0(A, R_0, ri)

if extra_plots:
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    ax.plot(r, phi)
    ax.plot(r, dphi_dt)
    ax.set_title('Initial profile of $\\phi(r)$ and $\\dot{\\phi}(r)$')
    ax.set_xlabel('r')
    ax.set_ylabel('Values')
    ax.legend(['$\\phi(0,r)$', '$\\dot{\\phi}(0,r)$'])
    save_figure(fig, Path(extra_plot_directory) / initial_profile_file)

# Solving the system of ODEs

init = np.concatenate([phi, dphi_dt])
dt = C_CFL*dr_min
t_eval = np.linspace(t0, T, time_samples)

solution = itg.solve_ivp(
    system,
    [t0, T],
    init,
    method=solver_method,
    max_step=dt,
    t_eval=t_eval,
    rtol=solver_rtol,
    atol=solver_atol,
)

cmap = plt.get_cmap(color_map)

fig, axs = plt.subplots(1, 2, figsize=(12, 5))
line, = axs[0].plot(solution.t[:end], solution.y[0,:end], label='$\\phi(t,0)$')
axs[0].set_title('Time evolution of $\\phi_0$')
axs[0].set_xlabel('t')
axs[0].set_ylabel('$\\phi(t,0)$')
axs[0].set_xlim([t0, T])
axs[0].set_ylim([-plot_amplitude_factor*A, plot_amplitude_factor*A])
params = Line2D([0], [0], color="none", marker="", linestyle="none", 
                label=rf"$A={A}, R_0={R_0}$"+"\n"
                +rf"$d={d}, m={m}, \lambda={mu}$"+"\n"
                +rf"$N={N}, R_{{max}}={R_max}, \alpha={alpha}$"+"\n"
                +rf"$R_{{sponge}}={R_sponge}, \gamma_0={gamma_0}, p={p}$")
axs[0].legend(handles=[params], handlelength=0, handletextpad=0)

for i in range(0, len(solution.t[:end]),
               max(1, int(len(solution.t[:end])/plot_profile_count))):
    axs[1].plot(r, solution.y[:N,i], color=cmap(1-i/len(solution.t[:end])), linewidth=1, label=f'$t={solution.t[i]:.2f}$')
axs[1].plot(r, solution.y[:N,-1], color=cmap(0), label=f'$t={solution.t[-1]:.2f}$')
axs[1].set_title('Spatial profile of $\\phi(t,r)$ at different timestamps')
axs[1].set_xlabel('r')
axs[1].set_ylabel('$\\phi(t,r)$')
axs[1].set_xlim([0, R_max])
axs[1].set_ylim([-plot_amplitude_factor*A, plot_amplitude_factor*A])
axs[1].legend(loc='upper right', fontsize=8)

save_figure(fig, output_file)

if show_plots:
    plt.show()
