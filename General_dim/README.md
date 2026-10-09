# Collapse scaling in general dimension

This standalone study extends the Module 4 radial scalar-field model to real
positive spatial dimension `d`. It is now located next to `module4/`, and does
not import code or data from the Module 4 results folders. The radial equation
derivation is in [`derivation.tex`](derivation.tex).

The default dimension sweep is `d = 1, 2, 3, 4, 5`. For every dimension the
script independently searches its own finite-horizon critical amplitude
`A_crit(d)`, then simulates one run at the same relative offset,
`A = 1.001 A_crit(d)` (`+0.1%`). The critical boundary is operationally defined
by reaching `|phi(0,t)| = 12` before `t = 100`; it is a finite-horizon numerical
threshold, not a claim about an exact mathematical critical amplitude. The
bisection bracket tolerance defaults to `1e-7`.

## Current run

All five selected runs reached `|phi(0,t)| = 80` and gave fits in the requested
window. The table reports one-standard-error fit uncertainties. The collapse
time is recalculated independently for each run. Since the amplitudes are
tuned to each dimension's own threshold, this is a near-critical comparison,
not a fixed-amplitude comparison.

| `d` | `A_crit` | run amplitude `A` | fitted `t_c` | `alpha` | `beta` from `R50` | `alpha + beta` |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.037210 | 1.038247 | 8.452499 | 0.9441 ± 0.0016 | 0.4154 ± 0.0021 | 1.3595 ± 0.0038 |
| 2 | 1.135935 | 1.137071 | 6.597640 | 0.8574 ± 0.0036 | 0.3909 ± 0.0021 | 1.2483 ± 0.0057 |
| 3 | 1.202610 | 1.203812 | 11.517886 | 0.8589 ± 0.0030 | 0.4472 ± 0.0001 | 1.3061 ± 0.0031 |
| 4 | 1.290820 | 1.292110 | 7.571426 | 0.6545 ± 0.0038 | 0.4437 ± 0.0019 | 1.0982 ± 0.0020 |
| 5 | 1.374216 | 1.375591 | 6.285726 | 0.7250 ± 0.0033 | 0.4915 ± 0.0021 | 1.2165 ± 0.0013 |

The fitted `t_c` is not constant across dimensions; it peaks at `d=3` in this
set of threshold-offset runs. The exponent sum also varies and is not forced
to equal one. Treat the displayed errors as regression-fit scatter only; they
do not capture shifts from changing the grid, critical horizon, offset, or
remaining-time window.

## Dense scan: step 0.2

A denser sweep is saved separately in
[`results/d_step_0.2/`](results/d_step_0.2/), with a run summary in its
[`README.md`](results/d_step_0.2/README.md), covering all 21 values from `d=1`
to `d=5` in increments of `0.2`. It uses the same `N=500`, 100-unit critical
horizon, `+0.1%` offset, central-field cutoff, and `[0.1, 0.3]` fit window as
the five-point run. The five integer-dimension trajectories and threshold
searches were reused from the matching run above; every intermediate
dimension's `A_crit` was independently re-bracketed and bisected. For those
intermediate searches, a linear interpolation of the known integer-dimension
thresholds only provided an initial bracket; it did not set the resulting
`A_crit`. Their bracket tolerance is `1e-5` (the reused integer brackets are
tighter).

In this finer scan, beta rises overall from about `0.415` at `d=1` to `0.492`
at `d=5`, while alpha declines overall but has a marked dip at `d=4`. The sum
also dips there to about `1.098`; that feature should be checked against grid
resolution before interpreting it physically. The independently fitted `t_c`
reaches its maximum, about `12.532`, near `d=2.6`.

- [Dense alpha/beta/sum plot](results/d_step_0.2/alpha_beta_sum_vs_dimension.png)
- [All 21 log-log fit panels](results/d_step_0.2/loglog_fits_by_dimension.png)
- [Dense collapse-time plot](results/d_step_0.2/collapse_time_vs_dimension.png)
- [Dense summary CSV](results/d_step_0.2/dimension_summary.csv)
- [Dense raw trajectory samples](results/d_step_0.2/trajectory_samples.csv)
- [Exact fit points](results/d_step_0.2/scaling_fit_points.csv)

The supercritical run continues to `|phi(0,t)| = 80`. Its collapse time is
estimated separately at each `d` by jointly fitting a free `t_c` and central
field power law to the final 40 regularly sampled points with `|phi(0,t)| >= 4`
(the same tail-fit convention used in the Module 4 beta study). With this
dimension-specific `t_c`, the exponent fits use only
`0.1 <= tau = t_c - t <= 0.3`:

```text
log10 |phi(0,t)| = -alpha log10(tau) + C_phi
log10 R50(t)     =  beta log10(tau) + C_R
```

`R50` is the first radius where `|phi(r,t)|` has fallen to half its central
value. The plot shows `alpha`, this half-central-field `beta`, and their sum
versus `d`. Small vertical bars are the ordinary-least-squares one-standard-
error slope uncertainties; the sum's uncertainty includes the covariance of
the two fits, which use the same time samples. These bars describe residual
fit scatter only, not systematic errors from resolution, critical-search
horizon, fit-window choice, or time-sample correlation.

Run from the project root:

```bash
python3 General_dim/general_dim.py
```

The script accepts alternate positive real dimensions, grid size, or a
different common offset. For example:

```bash
python3 General_dim/general_dim.py --dimensions 1 1.5 2 2.5 3 3.5 4 4.5 5
```

All dimensions share the same `results/` folder; no per-dimension directories
are created. Main outputs:

- [`alpha_beta_sum_vs_dimension.png`](results/alpha_beta_sum_vs_dimension.png):
  requested three-series scaling plot.
- [`loglog_fits_by_dimension.png`](results/loglog_fits_by_dimension.png): the
  central-field and `R50` log-log data with their separate fitted slopes in
  the selected interval.
- [`collapse_time_vs_dimension.png`](results/collapse_time_vs_dimension.png):
  the separately estimated `t_c(d)` values, with nonlinear-fit uncertainty.
- [`dimension_summary.csv`](results/dimension_summary.csv): per-dimension
  `A_crit`, chosen amplitude, `t_c`, `alpha`, `beta`, their errors, and fit
  quality.
- [`critical_summary.csv`](results/critical_summary.csv) and
  [`critical_search.csv`](results/critical_search.csv): per-dimension brackets
  and every finite-horizon bisection trial.
- [`trajectory_samples.csv`](results/trajectory_samples.csv): sampled central
  fields and `R50` for all dimensions in one file.
- [`scaling_fit_points.csv`](results/scaling_fit_points.csv): the exact points
  used by each log-log exponent fit.
- [`threshold_events.csv`](results/threshold_events.csv) and
  [`run_settings.csv`](results/run_settings.csv): interpolated field-threshold
  crossings and reproducibility settings.

The solver retains the Module 4 potential terms `-phi + phi^3`, RK4 time
stepping, stretched radial grid, outer sponge, and fixed outer boundary. A
non-integer `d` is an analytic continuation of the spherical radial operator,
not a Cartesian grid with a fractional number of axes.
