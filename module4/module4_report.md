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

The time series are ready for fits of the center field against time or of
`core_rms_radius` against `t_c - t`. The original script samples 10,000 times
over its configured run, while the Module 4 interval can be adjusted directly.
