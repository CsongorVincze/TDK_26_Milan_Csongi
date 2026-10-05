# Half-central-field radius: collapse-time sensitivity

Every plot here is a separate supercritical amplitude, at `N=1000` only. For
each run the code truncates the trajectory at successively higher
`|phi(0)|` crossings (coarse levels 20, 30, …, 80; finer one-unit levels are
added around regions where `alpha+beta` approaches 1). It estimates `tc` at
each crossing by jointly fitting `tc` and the central-field tail exponent,
then freely refits the notes' exponents in
`0.1 <= tau=tc-t <= 0.3`:

$$
\log_{10}|\phi(0,t)|=-\alpha\log_{10}(t_c-t)+C_\phi,\qquad
\log_{10}R_{50}(t)=\beta\log_{10}(t_c-t)+C_R.
$$

Each plot shows `alpha`, `beta`, and `alpha+beta` versus the joint-fit `tc`,
connected across threshold choices. The `tc` estimate comes from the high-field
tail after the trajectory reaches the selected central-field crossing; no
fixed-`alpha=1` estimator is used. Vertical bars are OLS fit standard errors;
horizontal bars are fit-based `tc` uncertainties.
These do not include systematic uncertainty from the fit window, threshold,
resolution, or the amplitude scan.

The trajectories are `N=1000`; the adopted `Acrit` value (`1.202609700896`)
was established in the earlier `N=500` finite-horizon search, not re-bracketed
at `N=1000` for this scan. Thus the reported percentage offsets are relative
to that reference.

The dense amplitude scan found `alpha+beta` close to 1 near `+0.03125%` for
the joint-tail estimator. For example, at the matched `|phi(0)|=26` crossing,
`tc=12.51870`, `alpha=0.65358`, `beta=0.35489`, and `alpha+beta=1.00847 ±
0.00231` for this radius. The apparent agreement is exploratory: the search
adaptively looked near 1, and nearby threshold choices shift the sum. See the
CSV for every fit rather than treating this single point as a universal
exponent measurement.

The +0.02875% trajectory did not reach even `|phi(0)|=20` by `t=100`, so it has
no threshold-fit plot. Its full trajectory is retained at
[`../tc_precision_N1000_plus_0.02875pct.csv`](../tc_precision_N1000_plus_0.02875pct.csv),
and the stop status for all runs is in
[`tc_precision_run_manifest.csv`](tc_precision_run_manifest.csv). Historic
`N<1000` data remain untouched but are excluded from this search.

## Plots by amplitude

`tc_sensitivity_N1000_plus_<offset>pct.png`:

[+0.01%](tc_sensitivity_N1000_plus_0.01pct.png) ·
[+0.02%](tc_sensitivity_N1000_plus_0.02pct.png) ·
[+0.02125%](tc_sensitivity_N1000_plus_0.02125pct.png) ·
[+0.0225%](tc_sensitivity_N1000_plus_0.0225pct.png) ·
[+0.02375%](tc_sensitivity_N1000_plus_0.02375pct.png) ·
[+0.025%](tc_sensitivity_N1000_plus_0.025pct.png) ·
[+0.02625%](tc_sensitivity_N1000_plus_0.02625pct.png) ·
[+0.0275%](tc_sensitivity_N1000_plus_0.0275pct.png) ·
[+0.03%](tc_sensitivity_N1000_plus_0.03pct.png) ·
[+0.03125%](tc_sensitivity_N1000_plus_0.03125pct.png) ·
[+0.0325%](tc_sensitivity_N1000_plus_0.0325pct.png) ·
[+0.03375%](tc_sensitivity_N1000_plus_0.03375pct.png) ·
[+0.035%](tc_sensitivity_N1000_plus_0.035pct.png) ·
[+0.03625%](tc_sensitivity_N1000_plus_0.03625pct.png) ·
[+0.0375%](tc_sensitivity_N1000_plus_0.0375pct.png) ·
[+0.03875%](tc_sensitivity_N1000_plus_0.03875pct.png) ·
[+0.04%](tc_sensitivity_N1000_plus_0.04pct.png) ·
[+0.05%](tc_sensitivity_N1000_plus_0.05pct.png) ·
[+0.0625%](tc_sensitivity_N1000_plus_0.0625pct.png) ·
[+0.075%](tc_sensitivity_N1000_plus_0.075pct.png) ·
[+0.0875%](tc_sensitivity_N1000_plus_0.0875pct.png) ·
[+0.1%](tc_sensitivity_N1000_plus_0.1pct.png)

The complete fit table is [`tc_precision_search.csv`](tc_precision_search.csv);
raw N=1000 trajectories are in the parent [`results/`](..) folder. The Gaussian
definition is plotted independently in its
[own results gallery](../../../../gaussian_profile/results/alpha_beta_scaling/README.md).
Rerun the cached search from the project root with:

```bash
python3 module4/Scale_invariance/tc_precision_search.py
```
