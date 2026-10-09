# Dimension sweep with step 0.2

This run covers `d = 1.0, 1.2, ..., 5.0` (21 dimensions). Every dimension uses
its own finite-horizon `A_crit` estimate and one run at `+0.1%` above it. All
runs reached `|phi(0,t)| = 80`; alpha and beta are fitted from 20 sampled points
in `0.1 <= t_c - t <= 0.3`. Beta uses the Module 4 half-central-field radius
`R50`.

The five integer-dimension runs reuse the matching simulations and critical
searches from the parent `results/` folder. Each intermediate dimension has
its own collapse/no-collapse bracket and bisection; interpolation from the
integer thresholds only seeds the initial bracket (`+/-0.025`). The requested
bisection tolerance for new dimensions was `1e-5`; the reused integer brackets
are tighter. This keeps the physical and fit settings matched while avoiding
duplicate endpoint runs.

Alpha decreases overall from `0.944` at `d=1` to about `0.725` at `d=5`, with
a pronounced dip at `d=4`. Beta rises overall from `0.415` to `0.492`. The
sum `alpha + beta` ranges from about `1.098` at `d=4` to `1.360` at `d=1`;
it is not constrained to one. The `d=4` alpha dip is a feature to check with a
resolution study before assigning it physical significance. `t_c` is fitted
independently at every dimension and reaches about `12.532` near `d=2.6`.

The main figure's vertical bars are one-standard-error OLS slope uncertainties;
the sum's uncertainty includes the covariance of alpha and beta. These errors
describe fit residuals, not systematic uncertainty from resolution or fit
choices.

- [Alpha, beta, and sum versus dimension](alpha_beta_sum_vs_dimension.png)
- [Log-log fits for all 21 dimensions](loglog_fits_by_dimension.png)
- [Individually fitted collapse times](collapse_time_vs_dimension.png)
- [Dimension summary and fit statistics](dimension_summary.csv)
- [Critical brackets and trials](critical_summary.csv) · [full search history](critical_search.csv)
- [Trajectory samples](trajectory_samples.csv) · [fit points](scaling_fit_points.csv)
- [Settings](run_settings.csv)

Rerun from the project root with:

```bash
python3 General_dim/general_dim.py \
  --dimensions 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8 4.0 4.2 4.4 4.6 4.8 5.0 \
  --critical-tolerance 1e-5 \
  --reuse-from General_dim/results \
  --output-dir General_dim/results/d_step_0.2
```
