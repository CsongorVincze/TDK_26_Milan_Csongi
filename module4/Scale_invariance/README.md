# Module 4 scale-invariance studies

These are standalone collapse studies with two different core-radius
definitions:

- [Half-central-field radius](Fraction_of_central/half_radius/README.md),
  defined by the first point where the field falls to half its central value.
- [Gaussian-profile radius](gaussian_profile/README.md),
  inferred from a time-dependent fit to the inner field profile.

## Current threshold and collapse-time sensitivity search

The current paired search retains only spatial resolution `N=1000`; no
`N<1000` result enters its fit tables or plots. It scans supercritical
amplitudes and several central-field thresholds, estimating `tc` at each
threshold before fitting the notes' exponents in `0.1 <= tc-t <= 0.3`.
The same central-field data and collapse time are used for both radius
definitions, so the radius method is the intended comparison.

- [Half-radius results and plots](Fraction_of_central/half_radius/results/alpha_beta_scaling/README.md)
- [Gaussian-radius results and plots](gaussian_profile/results/alpha_beta_scaling/README.md)
- [Shared plot index](alpha_beta_scaling/README.md)
- [Search implementation](tc_precision_search.py)

Each amplitude gets one figure using only a joint fit of `tc` and the
central-field exponent on the high-field tail after each threshold crossing.
The plotted near-collapse `alpha` is then freely fitted in the stated time
window, while `beta` is fitted from the radius. Fit tables and run manifests
are saved beside the corresponding method's plots. Raw trajectories are kept
in that method's own `results/` folder so the studies remain independent.

The older [shared-`tc` baseline](shared_tc_alpha.py) is retained as a
single-run reference and now uses only the joint high-field-tail fit for
`tc` and `alpha`; the fixed-`alpha=1` estimator is retired. Broader N=500
exploratory data remain in their original folders but are excluded from this
N=1000 sensitivity search.
