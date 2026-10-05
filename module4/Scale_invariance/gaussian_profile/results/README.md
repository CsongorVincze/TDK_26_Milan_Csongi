# Gaussian-profile outputs

- `supercritical_plus_*.csv`: sampled time series, including the fitted
  Gaussian width, implied half-height radius, profile-fit `R^2`, and fit-point
  count.
- `tc_precision_N1000_plus_*.csv`: additional N=1000 trajectories for the
  multi-amplitude collapse-time sensitivity scan, including a run that did not
  reach `|phi(0)|=20` by the `t=100` horizon.
- `case_index.csv` and `run_settings.csv`: run amplitudes, cutoffs, grid
  settings, and the critical-amplitude value used.
- `gaussian_profile_snapshots.csv`: selected radial field profiles and their
  Gaussian fits at central-field levels 4, 8, 12, 20, 40, and 80.
- `beta_fit_sensitivity.csv`: log-log slopes for each `tc` method and fit
  window.
- `width_vs_central_field.csv`: direct `q=beta/gamma` fits and Gaussian-shape
  quality summaries.
- `shared_tc_alpha_fits.csv`: central-field `alpha`, Gaussian-radius `beta`,
  and radius normalization `C_R`, using the same per-run `tc` and fit window as
  the half-radius comparison. The focused comparison retains `N >= 1000` and
  `0.1 <= tc-t <= 0.3`; `alpha` is not the radius prefactor.
- PNG files: contraction and Gaussian profile-fit quality. Current paired α–β
  fit tables and one collapse-time sensitivity plot per amplitude are in this
  results folder's [`alpha_beta_scaling`](alpha_beta_scaling/README.md)
  gallery; matching half-radius figures use the same central-field data, `tc`,
  and fit window.

Regenerate the simulations and outputs with `study.py`; rerun only the fits
and figures with `study.py --analyze-existing`.
