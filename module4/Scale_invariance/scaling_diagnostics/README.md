# Scaling diagnostics

This folder audits two issues that matter when interpreting the current
`alpha + beta` fits: whether the N=1000 collapse outcome is monotone in the
initial amplitude, and how strongly the fitted exponents depend on the
remaining-time window and field cutoff.

Start with the [findings report](report.md). The main figures are:

- [Finite-horizon outcome map](results/n1000_finite_horizon_outcomes.png)
- [Fit-window and cutoff sensitivity](results/fit_window_sensitivity_N1000.png)

The corresponding tables are
[`n1000_finite_horizon_outcomes.csv`](results/n1000_finite_horizon_outcomes.csv)
and
[`fit_window_sensitivity_N1000.csv`](results/fit_window_sensitivity_N1000.csv).
The recorded-step local transition search is documented in
[`critical_resolution_check.py`](critical_resolution_check.py); its two output
CSVs are kept alongside the plots.

Both audit scripts use saved trajectories and do not rerun the PDE:

```bash
python3 module4/Scale_invariance/scaling_diagnostics/outcome_audit.py
python3 module4/Scale_invariance/scaling_diagnostics/fit_window_audit.py
```

The second script uses four remaining-time windows, all bounded by `tau=0.4`,
and the cutoff levels `|phi(0)|=20, 24, 40, 60, 80`. It retains only N=1000
data and compares the three original coarse amplitude offsets with the
search-selected `+0.03125%` case.

For the follow-up high-field audit of the Module 5 derivation, see
[module5_audit/report.md](module5_audit/report.md). It reruns the same solver
to `|phi(0)|=1280` across resolutions and parameter variations, and directly
compares the core PDE terms, `alpha + beta`, and gradient scaling.
