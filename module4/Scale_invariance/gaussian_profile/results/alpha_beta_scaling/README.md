# Gaussian-profile radius: collapse-time sensitivity

Every plot here is a separate supercritical amplitude, at `N=1000` only. For
each run the trajectory is truncated at successive `|phi(0)|` crossings
(coarse levels 20, 30, …, 80, with one-unit refinement near regions where
`alpha+beta` approaches 1). At each crossing, `tc` and the tail central-field
exponent are estimated jointly from the high-field tail; the notes' exponents
are then freely fitted in `0.1 <= tau=tc-t <= 0.3`:

$$
\log_{10}|\phi(0,t)|=-\alpha\log_{10}(t_c-t)+C_\phi,\qquad
\log_{10}R_{50,G}(t)=\beta\log_{10}(t_c-t)+C_R.
$$

Each plot shows `alpha`, `beta`, and `alpha+beta` versus the joint-fit `tc`,
with threshold-choice points connected. There is no fixed-`alpha=1` estimator
in this search. Vertical error bars are OLS slope standard errors; horizontal
bars are fit-based `tc` uncertainties. They do not include systematic effects
from changing the fit window, cutoff, resolution, or amplitude.

All collapse trajectories use `N=1000`. The adopted `Acrit` value
(`1.202609700896`) is from the earlier `N=500` finite-horizon search and was
not re-bracketed at `N=1000` here; the percentage offsets are relative to that
reference.

Near `+0.03125%` the joint-tail estimates put the sum close to 1. At the
matched `|phi(0)|=26` crossing, `tc=12.51870`, `alpha=0.65358`,
`beta=0.34446`, and `alpha+beta=0.99804 ± 0.00243`. This is suggestive, not a
settled scaling law: the amplitude scan was explicitly refined around a
near-one result, and nearby thresholds and the radius definition shift the
sum. Compare the paired values in the half-radius gallery and use the complete
CSV for the sensitivity range.

The +0.02875% trajectory did not reach `|phi(0)|=20` by `t=100`, so no
threshold-based fit or plot exists for it. The raw trajectory is retained at
[`../tc_precision_N1000_plus_0.02875pct.csv`](../tc_precision_N1000_plus_0.02875pct.csv);
all crossing status is recorded in
[`tc_precision_run_manifest.csv`](tc_precision_run_manifest.csv). Historic
`N<1000` data remain in place but are not included in this search.

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
raw N=1000 trajectories are in the parent [`results/`](..). The paired script
loads both standalone radius solvers to guarantee a matched comparison; it
does not use the old beta-study folder. Rerun it from the project root with:

```bash
python3 module4/Scale_invariance/tc_precision_search.py
```
