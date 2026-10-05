# Gaussian-profile core-width experiment

This is a standalone Module 4 simulation/analysis. It has its own radial-grid
solver, evolving Gaussian-width fit, near-critical supercritical runs, CSVs,
plots, and report; it imports neither the beta study nor the half-radius study.

From the project root, run:

```bash
python3 module4/Scale_invariance/gaussian_profile/study.py
```

By default it uses the previously determined N=500 finite-horizon reference
`Acrit = 1.202609700895846`, runs 0.01%, 0.05%, and 0.1% supercritical cases,
continues the closest case through `|phi(0)|=80`, and repeats that closest case
at N=1000. Add `--search-critical` to independently recompute the finite-time
critical bracket with this folder's own solver. Use `--analyze-existing` to
rebuild analysis products from this folder's CSVs without rerunning the PDE.
All products are written to [results](results/), independent of the current
working directory.

See [report.md](report.md) for the Gaussian ansatz, radius definition, fit
results, and its limitations.

At each sample the script fits the normalized inner profile over
`0.30 <= |phi(r)|/|phi(0)| <= 0.90` to
`exp(-r^2/(2 sigma^2))`. It reports `sigma` and the equivalent Gaussian
half-height radius `R50,G = sqrt(2 ln 2) sigma`; the fit-quality column and
profile overlays show where the Gaussian-shape assumption is credible.

For the current direct comparison with the half-height radius, run
`python3 module4/Scale_invariance/tc_precision_search.py`. The paired script
loads both standalone solvers, uses `N=1000` collapse trajectories, and writes
matched multi-amplitude, multi-threshold results into each study's own
`results/alpha_beta_scaling` folder. The one-amplitude shared-`tc`
baseline can still be regenerated with
`python3 module4/Scale_invariance/shared_tc_alpha.py`; its matched run is
`N=1000`, `+0.01%`. It uses only the joint high-field-tail fit for `tc` and
`alpha`; the retired fixed-`alpha=1` path is not regenerated. Radius
normalizations are labeled `C_R`.
