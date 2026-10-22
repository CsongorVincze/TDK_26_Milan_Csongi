# Project handoff

## Goal

Study oscillon collapse in the radial scalar-field model, especially near-critical
scaling with remaining time `tau = t_c - t`. The main observables are `alpha`
(central-field scaling) and `beta` (core-radius scaling); compare their values
and their sum without assuming the sum is one. The current focus is extending
the analysis over spatial dimension `d`, including non-integer values.

## Work completed

- Module 4 contains the spherical 3D collapse solver, critical-amplitude and
  supercritical studies, half-central-field (`R50`) and Gaussian core-radius
  scaling studies, plus diagnostics and CSV exports.
- `General_dim/` is a standalone general-dimension implementation beside
  `module4/`; its derivation is in `General_dim/derivation.tex`.
- General-dimension sweeps cover `d=1..5` at integer steps and at `0.2` steps.
  The dense outputs are in `General_dim/results/d_step_0.2/`; see its README,
  exponent plot, 21 log-log panels, summary CSV, and raw trajectory CSV.
- Dense-sweep settings: `N=500`, `R0=2`, critical event `|phi(0,t)|=12` by
  `t=100`, one run at `+0.1%` over each dimension's `A_crit`, continued to
  `|phi(0,t)|=80`. Each `t_c` is independently fitted from the last 40
  central-field samples with `|phi(0,t)| >= 4`. Alpha and beta are fitted over
  `0.1 <= t_c-t <= 0.3`; beta uses the first radius where the field magnitude
  falls to half its central value (`R50`). All 21 runs fit successfully.
- The dense scan shows a marked alpha dip near `d=4`; check its resolution
  dependence before interpreting it. Fit error bars are OLS residual errors,
  not systematic or resolution uncertainties.

## Important model distinction

`General_dim/` generalizes the Module 4 potential `V(phi)=phi^2/2-phi^4/4`,
retaining both mass and cubic-force terms, RK4, the stretched radial grid, and
the outer sponge. The original `milan_alapkod/milan_alapkod.py` uses a shifted
double-well potential and was used for long-time oscillon decay; do not conflate
that decay study with Module 4's finite-time-collapse scaling.

## Continue from here

- Read [`General_dim/results/d_step_0.2/README.md`](General_dim/results/d_step_0.2/README.md)
  and [`General_dim/README.md`](General_dim/README.md) first.
- Reproduce the dense run with the command in its README. Intermediate `A_crit`
  values were independently bracketed and bisected (the neighboring integer
  results only seed the initial bracket); integer endpoints reuse matching
  saved runs.
- A useful next validation is a resolution study of the `d=4` alpha dip and
  sensitivity to the finite-horizon critical threshold / fit window.

## Handoff cautions

- The user has explicitly required approval before any commit or push. Nothing
  in this work should be committed or pushed without that approval.
- At handoff, `git status` reports deletions of
  `milan_alapkod/fourier_decay_fit.py`, two Fourier/decay PNGs under
  `milan_alapkod/results/`, and
  `module4/results/module4_critical_central_field.png`. These were not changed
  as part of the general-dimension sweep; inspect with the user before
  restoring or committing them.
