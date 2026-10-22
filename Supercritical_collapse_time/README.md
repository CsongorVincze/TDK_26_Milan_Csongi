# Dense supercritical amplitude sweep

**Latest precision and t=100 audit:** see [the detailed report](results/precision_longtime/REPORT.md).
Plots now use the raw ratio `(A-Ac)/Ac`. The original reference and search
below are retained as provenance; the refined N=500 reference is
`1.2026097002647038`, and spatial convergence remains unresolved.

This study reruns the Module 4 **3D** Gaussian critical-amplitude search and
reproduces the earlier two-panel central-field/collapse-cutoff-time plot.
All new source, logs, figures, and numerical outputs live in this folder.

## Refined critical amplitude

For the original N=500 spatial grid and CFL=0.4:

- Non-crossing lower edge: **1.202609700203175**.
- Crossing upper edge: **1.202609700302128**.
- Midpoint reference: **1.2026097002526517**.
- Bracket width: **9.8953e-11**; midpoint half-width: **4.9477e-11**.
- Both endpoints round to **1.202610** at six decimal places.

The definition is reaching `|phi(0,t)|=12` by `t=100`. The upper edge reaches
that cutoff at `t=21.1989075357`; the lower edge does not. Both endpoint
classifications remain the same at CFL=0.2, and at horizons 100 and 200.
The interval is a refined **local finite-horizon numerical transition**.
Its width measures bisection precision for this discretization, not error
relative to a continuum critical solution. Six-decimal continuum accuracy is
not established. Existing N=1000 diagnostics in
`../module4/Scale_invariance/scaling_diagnostics/report.md` already show
resolution dependence and non-monotone outcomes.

The midpoint itself does not reach the cutoff by t=100; it is retained as the
zero-offset reference, with its full trajectory saved. The two bracket-edge
trajectories are also saved through t=200 (or their cutoff).

## Sweep and figures

The amplitude is `A=Acrit*(1+delta)` for `0 <= delta <= 0.2`.
The initial 342 runs include zero, 181 logarithmically spaced points from
`delta=1e-9` to `0.01`, and 161 linearly spaced points from `0.01` to `0.2`,
with the repeated endpoint removed. Two initially non-crossing positive
offsets receive additional 81-point local scans within +/-10% of each offset.
The final sweep contains **502 amplitudes** (457 cutoff crossings and 45
non-crossings, including the midpoint reference); the count is also recorded
in `results/settings.json`.

The sampled non-crossing positive offsets lie in two bands: approximately
**3.13068e-5 to 3.37400e-5**, and **2.36272e-3 to 2.43432e-3** (raw ratios). These
are sampled ranges, not separately bisected band boundaries. Outcomes are
therefore demonstrably non-monotone even on the N=500 grid.

- [Reproduced two-panel plot](results/supercritical_sweep.png).
- [Cutoff times distinguished by final field sign](results/collapse_time_detail.png).
- [Full-horizon non-crossing trajectories](results/noncrossing_trajectories.png).
- [Every amplitude and its outcome](results/sweep_summary.csv).
- [Critical search trials](results/critical_search_N500_cfl0.4.csv).
- [Critical bracket](results/critical_summary.json).
- [Time-step and horizon checks](results/critical_sensitivity.csv).
- [Independent solver comparison](results/solver_validation.csv).
- [Complete settings](results/settings.json).
- Raw central-field samples: `results/trajectories/run_*.csv`.

The main figure zooms the time axes to show collapse trajectories clearly.
Crosses at the top of the right panel denote **no crossing by t=100**, not
measured event times. Zero offset is drawn as a black trajectory on the left;
it cannot appear on the logarithmic offset axis. The separate full-horizon
figure shows the surviving trajectories through t=100. Narrow non-crossing
windows above the reference mean “supercritical” here labels amplitude above
the chosen reference and does not guarantee a cutoff crossing.

The rightmost amplitude is **1.443131640303182**, exactly 20% above the
reference to floating-point precision. Its cutoff time is **3.74557653367**.
These times are field-cutoff proxies, not extrapolated singularity times.

## Model and numerical method

`phi_tt = phi_rr + (2/r) phi_r - phi + phi^3 - gamma(r) phi_t`,
with `phi(r,0)=A exp(-r^2/(2 R0^2))`, `R0=2`, and zero initial velocity.
The grid, origin stencil, outer boundary, sponge, adaptive step bound and RK4
scheme match `module4/milan_alapkod_module4.py`: N=500, Rmax=25, stretch=5,
sponge starts at 0.8 Rmax with strength 7 and power 3, and CFL=0.4.

`solver.cpp` accelerates the existing equations without fast-math compiler
flags. Two independent comparisons against the original unsampled Python
solver agree in cutoff time within 1.5e-13. Recorded central-field values are
linearly interpolated every 0.02 time units without changing integration
steps; the exact interpolated cutoff point is appended. This avoids the
original solver's recording-dependent step clipping, so search and sweep
use the same integration path. Near-critical trajectories can consequently
differ from old saved trajectories that clipped time steps to output times.
Nonfinite states raise an error rather than being classified as collapse.

## Reproduce

From the project root (requires g++, Python, NumPy, and Matplotlib):

```bash
MPLCONFIGDIR=/tmp/supercritical-mpl python3 Supercritical_collapse_time/study.py
```

The script compiles the C++ solver, validates it against Module 4, reruns the
bisection to 1e-10, performs the initial sweep and local refinements, makes the
figures, and checks bracket sensitivity. Six concurrent solver processes are
used for independent sweep amplitudes. No SciPy or additional packages are
needed. The main implementation and all output settings are kept here.

## English TDK figure style

All three figures use **TDK_plot_skill**: English labels and decimal points,
consistent navy typography, pale blue-gray panels, subtle grids, and legible
colors and marker shapes. Each figure is exported as a 300 dpi PNG only.
The non-crossing trajectories are grouped into three comparable panels.
The project-local `tdk.mplstyle` preserves the style for reproducibility.

Regenerate all figures from existing CSV files without rerunning the solver:

```bash
MPLCONFIGDIR=/tmp/supercritical-mpl python3 Supercritical_collapse_time/study.py --plot-only
```
