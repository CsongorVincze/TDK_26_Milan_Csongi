# Half-radius results

Every trajectory CSV records `time`, central field `phi_center`, its magnitude,
central time derivative, `r50`, and `r50_squared`. Here `r50` is the first
interpolated radius where `abs(phi(r,t)) <= 0.5*abs(phi(0,t))`; it is not an
energy-weighted radius. The files retain the time-series points used for the
fits, so the collapse windows can be refit without rerunning the PDE.

- `run_settings.csv`, `critical_summary.csv`, and
  `critical_threshold_search.csv`: reproducibility settings, finite-horizon
  threshold, and its bisection trials.
- `critical_lower_boundary.csv` / `critical_upper_boundary.csv`: two bracket
  trajectories.
- `supercritical_plus_*.csv` and `supercritical_sweep_summary.csv`: three
  initial amplitudes above the bracket midpoint, through `|phi(0)|=20`.
- `continuation_*.csv`: the closest supercritical run continued to larger
  central-field cutoffs; includes an `N=1000` run through 80.
- `tc_precision_N1000_plus_*.csv`: additional `N=1000` trajectories used by the
  paired collapse-time sensitivity scan; includes the +0.02875% run that did
  not reach `|phi(0)|=20` by `t=100`.
- `resolution_*.csv` and `resolution_summary.csv`: fixed-amplitude comparison
  at `N=250, 500, 1000` through `|phi(0)|=20`.
- `beta_fit_sensitivity.csv`: `R50` log-log time fits over several remaining-
  time windows, using the joint high-field `(tc, alpha)` estimate; includes
  supercritical, continuation, and resolution runs.
- `width_vs_center_scaling.csv`: direct exponent in `R50 ~ |phi(0)|^(-q)`,
  which does not require estimating `tc`.
- `shared_tc_alpha_fits.csv`: central-field `alpha`, half-radius `beta`, and
  radius normalization `C_R`, all fitted with a shared per-run `tc` and
  matched remaining-time window. The focused comparison retains `N >= 1000`
  and `0.1 <= tc-t <= 0.3`; `alpha` is not the radius prefactor.

PNG figures visualize the threshold search, central field, half-radius over
time, continuation, and resolution comparison. Current paired α–β fit tables
and one plot per amplitude are in this results folder's
[`alpha_beta_scaling`](alpha_beta_scaling/README.md) gallery; matching
Gaussian-radius figures use the same central-field data, `tc`, and fit window.
