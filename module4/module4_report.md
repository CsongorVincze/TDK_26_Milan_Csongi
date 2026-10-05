# Module 4 simulation and results

This report documents the implementation based on the Module 4 questions in
[Oscillon_dynamics_notes_20260917.pdf](../Oscillon_dynamics_notes_20260917.pdf).

## What changed

[milan_alapkod_module4.py](milan_alapkod_module4.py) is a focused version of
the original simulation for the single-well potential
`V(phi) = 0.5 phi^2 - 0.25 phi^4`, with derivative
`V'(phi) = phi - phi^3`.

It evolves the spherically symmetric 3D equation on the original stretched
radial grid, with the outer sponge and a fourth-order Runge–Kutta integrator.
The initial data are `phi(r, 0) = A0 exp(-r^2 / (2 R0^2))` and
`dphi/dt(r, 0) = 0`, using the requested `R0 = 2`.

The script now has one compact workflow for threshold bisection, a supercritical
run, energy-weighted core-radius data, collapse-time estimation, plots, and a
JSON summary. The original [double-well script](../milan_alapkod/milan_alapkod.py)
has also been reorganized into named stages; its potential, stretched grid,
sponge and boundary treatment, RK45 setup, energy diagnostics, plots, and CSV
output are retained.

## Method and assumptions

- Threshold bisection classifies a run as collapsing when
  `abs(phi(0, t))` reaches 12 before `t = 40`. The reported threshold is thus an
  operational finite-time threshold; it is not proof of the infinite-time
  boundary.
- `--mode critical` searches again using `--critical-tmax` (default 100) and
  `--critical-tolerance` (default `1e-8`). It simulates both the bracket
  midpoint `Acrit` and the upper (collapsing) bracket edge, then plots them
  together. The upper edge is within the configured tolerance of `Acrit` and
  reaches the collapse cutoff; the exact midpoint may or may not collapse
  within the finite run time.
- The notes define `R_core` but do not specify its value. The default is a
  fixed `R_core = 8`, which is recorded with each run.
- The moment uses the signed local energy density
  `rho = 0.5 (dphi/dt)^2 + 0.5 (dphi/dr)^2 + V(phi)`. The measured quantity is
  `<r^2>(t) = integral_0^Rcore [r^4 rho(r,t) dr] / integral_0^Rcore [r^2 rho(r,t) dr]`.
- The code estimates `t_c` by fitting `1 / abs(phi(0,t))` against time in the
  high-amplitude tail, then fits `sqrt(<r^2>)` against `t_c - t` on a log-log
  scale. It accepts `beta` only for a positive slope and `R^2 >= 0.9`; otherwise
  it reports the raw slope as a diagnostic and labels `beta` undetermined.

## Results

| Quantity | Result |
|---|---:|
| Critical amplitude, `N=500` | `A_crit = 1.2026` |
| Bisection bracket, `N=500` | `[1.2026031, 1.2026291]` |
| Critical amplitude, `N=250` check | `1.2026` to four decimals |
| Supercritical amplitude used | `A0 = 1.2627599` |
| Collapse-proxy crossing time | `4.472967` |
| Estimated `t_c` | `4.730324`, inverse-amplitude fit `R^2 = 0.98944` |
| Raw log-log slope for the prescribed radius | `0.00547`, fit `R^2 = 0.04681` |

The separate near-critical search (`--critical-tmax 100`, tolerance `1e-8`)
estimated `Acrit = 1.2026097009`, with bracket
`[1.2026096977, 1.2026097041]`. The midpoint run reached the collapse cutoff
at `t = 20.685919`; the upper bracket edge did so at `t = 20.396726`.
These are horizon-100 results and are distinct from the horizon-40 threshold
reported in the table above.

The two grid resolutions agree at the requested four decimal places, but their
bisection brackets differ. At the shared bracket endpoint, the resolutions
classify the amplitude differently, so finer threshold digits are not
converged.

The core radius in the selected fit window first decreases and then increases
(approximately 2.91 to 2.83 to 2.92). The log-log fit is poor, so the data do
not support a reliable contraction exponent `beta` for the signed
energy-weighted moment. The negative quartic term makes the local energy
density signed during collapse, which can make this weighted radius sensitive
to cancellations. A different positive weighting or a moving core boundary
would define a different observable and would need to be justified separately.

An independent review checked the potential, PDE, radial Laplacian, moment
formula, threshold bracket, and the fit-quality conclusion.

## Running and output files

From the project root, run the full Module 4 workflow with:

```bash
python module4/milan_alapkod_module4.py
```

Use `--mode threshold` to run only the threshold search. Module 4 outputs are
saved under [module4/results](results/):

- [Run summary](results/module4_summary.json)
- [Threshold search data](results/module4_threshold_search.csv)
- [Collapse time series](results/module4_collapse_data.csv)
- [Core-radius fit plot](results/module4_core_contraction.png)
- [Central-field plot](results/module4_central_field.png)
- [Threshold-search plot](results/module4_threshold_search.png)

To inspect a delayed, near-critical collapse, run:

```bash
python module4/milan_alapkod_module4.py --mode critical
```

This performs a fresh threshold search with a 100-unit horizon, then evolves
the threshold midpoint and closest tested supercritical amplitude. It writes
`module4_critical_midpoint_data.csv`, `module4_critical_data.csv`,
`module4_critical_summary.json`, `module4_critical_central_field.png`, and a
separate critical threshold-search plot in `module4/results/`. Increase the
horizon with `--critical-tmax 200` if
you want to examine a longer-lived solution; this also makes the threshold
search more expensive. Tighten `--critical-tolerance` further if you want to
probe closer to the threshold (and accept a longer bisection).

Current critical-run outputs: [comparison plot](results/module4_critical_central_field.png),
[midpoint data](results/module4_critical_midpoint_data.csv),
[collapsing-edge data](results/module4_critical_data.csv), and
[critical-run summary](results/module4_critical_summary.json).

## Supercritical amplitude sweep

To see how collapse changes above the horizon-100 critical estimate, the new
`--mode sweep` simulates eight amplitudes from 0.01% through 25% above
`Acrit = 1.2026097009`. The 25% endpoint is well above the earlier default
supercritical run (about 5% above the estimate). These runs reuse the existing
critical summary instead of repeating the threshold search; omitting
`--critical-amplitude` makes the sweep perform a fresh search.

| Excess over Acrit | Initial A | Time to `abs(phi(0,t)) = 12` |
|---:|---:|---:|
| 0.01% | 1.202729962 | 15.14935 |
| 0.05% | 1.203211006 | 11.76284 |
| 0.1% | 1.203812311 | 11.39535 |
| 0.5% | 1.208622749 | 7.99400 |
| 1% | 1.214635798 | 7.64010 |
| 5% | 1.262740186 | 4.47320 |
| 10% | 1.322870671 | 4.06919 |
| 25% | 1.503262126 | 3.69545 |

All eight trajectories reached the Module 4 collapse cutoff within the
100-time-unit horizon. In this sweep, the measured collapse time decreases as
the amplitude is raised. The combined [sweep plot](results/module4_supercritical_sweep.png),
[event-time/data index](results/module4_supercritical_sweep_summary.csv), and
[run settings](results/module4_supercritical_sweep_summary.json) are saved in
`module4/results/`; each amplitude also has its own time-series CSV.

The sweep now also creates a beta-fit plot for every trajectory and an
eight-panel [comparison figure](results/module4_supercritical_beta_fits.png).
These figures plot the transformed coordinates explicitly:
`log10(t_c - t)` horizontally and `log10(sqrt(<r^2>))` vertically (the plotted
axes themselves are linear). This makes the straight-line fit visible without
relying on logarithmic tick scaling.
For each run, `t_c` is estimated from the high-amplitude inverse-central-field
tail, then the same log-log fit is applied to the last 30 resolved samples with
`2*dt <= t_c - t <= 1` and `abs(phi_center) >= 1.2`:

```text
log10(sqrt(<r^2>)) = beta*log10(t_c - t) + constant
```

The slope `beta` is independent of the logarithm base; the code and plots now
use base 10 consistently so the plotted coordinates are exactly those fitted.

The [beta summary CSV](results/module4_supercritical_beta_summary.csv) records
the slope, fit quality, estimated `t_c`, and each individual plot. The current
slopes below are diagnostics only; none passes the existing acceptance rule
(`beta > 0` and fit `R^2 >= 0.9`):

| Excess over Acrit | Diagnostic beta | Fit R^2 | Status |
|---:|---:|---:|---|
| 0.01% | 0.0037 | 0.039 | undetermined |
| 0.05% | 0.0391 | 0.740 | undetermined |
| 0.1% | -0.0126 | 0.271 | undetermined |
| 0.5% | 0.0309 | 0.533 | undetermined |
| 1% | -0.0043 | 0.029 | undetermined |
| 5% | 0.0056 | 0.048 | undetermined |
| 10% | -0.0463 | 0.851 | undetermined |
| 25% | -0.1103 | 0.994 | undetermined |

In particular, the 25%-above-threshold run is nearly linear in log-log space,
but its negative slope means the measured core RMS radius grows toward the
collapse cutoff rather than shrinking. It therefore is not a positive
critical-collapse exponent under this observable.

## Near-critical subcritical animation

The [animation](visualizations/subcritical_near_critical/subcritical_core_evolution.gif)
shows a separate run just below the saved critical bracket. It plots the radial
field profile alongside `phi(0,t)` and the signed-energy RMS radius. This run
does not reach the `|phi(0,t)| = 12` cutoff through `t = 100`; its largest
central-field magnitude is `5.871` at `t = 20.10`.

In the displayed window, the RMS radius drops from about `2.93` at `t = 14.20`
to a minimum of `2.05` at `t = 20.54`, then grows to `4.74` by `t = 32.00`.
So this near-critical subcritical solution does show a contracting interval,
followed by a turn-around and expansion. The [CSV](visualizations/subcritical_near_critical/subcritical_near_critical_timeseries.csv),
[run summary](visualizations/subcritical_near_critical/subcritical_near_critical_summary.json),
[animation script](create_subcritical_animation.py), and a short
[README](visualizations/subcritical_near_critical/README.md) are saved with it.
The radial profile is included because the signed-energy moment can be
misleading on its own.

## Short continuation past the `phi = 12` cutoff

As a controlled check, the nearest supercritical sweep run (`A0 = 1.202729962`)
was continued to `|phi(0,t)| = 20`, leaving the default cutoff unchanged. The
RMS radius was `2.907` at the `phi = 12` event and `3.022` at `phi = 20`
(`mean_r2` grows from `8.453` to `9.135`); thus the signed-energy radius still
grew during this extra interval. To compare cutoff choices fairly, both runs
were sampled at `0.01`. The inverse-amplitude estimate of `t_c` is `15.351`
with cutoff 12 and `15.321` with cutoff 20. Their high-amplitude fit R² values
are both about `0.995`, but their beta diagnostics differ (`-0.0158` versus
`-0.0439`), and neither is accepted. Since the chosen tail depends on the
cutoff, this is a sensitivity warning, not evidence for a stable exponent.
The [cutoff-12 reference CSV and plots](visualizations/supercritical_continue_to_phi20/cutoff_12_reference/)
and [continued-to-20 data and plots](visualizations/supercritical_continue_to_phi20/)
are saved separately. This is only a short, fixed-grid diagnostic; continuing
further should be checked for spatial-resolution convergence.

Recreate the sweep with the saved critical estimate using:

```bash
python module4/milan_alapkod_module4.py --mode sweep \
  --critical-amplitude 1.202609700895846
```

Run the original double-well workflow from the project root with
`python milan_alapkod/milan_alapkod.py`. Its CSV and figures are saved under
`milan_alapkod/results/` (with the additional diagnostic figures in
`milan_alapkod/results/extra_plots/`).

## CSV data for later fitting

Both scripts now export the same time-series columns:
`time`, `phi_center`, `pi_center`, `total_energy`, `core_energy`, `mean_r2`,
and `core_rms_radius`. The radius fields are computed from the signed local
energy density and the configured core cutoff.

- [milan_alapkod.py](../milan_alapkod/milan_alapkod.py) writes one row for every
  stored solver time to a dated
  `milan_alapkod_data_YYYYMMDD_HHMMSS.csv` file under
  `milan_alapkod/results/`.
- [milan_alapkod_module4.py](milan_alapkod_module4.py) writes one row per saved
  sample, plus the collapse
  endpoint, to `module4/results/module4_collapse_data.csv` by default. Its
  default interval is 0.02; use `--sample-interval` to request denser samples,
  or `--csv-name` to choose another filename.
- `--mode critical` writes the near-critical run to
  `module4/results/module4_critical_data.csv` and the exact midpoint run to
  `module4/results/module4_critical_midpoint_data.csv`, with the same columns.

The time series are ready for fits of the center field against time or of
`core_rms_radius` against `t_c - t`. The original script samples 10,000 times
over its configured run, while the Module 4 interval can be adjusted directly.
