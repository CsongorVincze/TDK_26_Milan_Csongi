# Module 5 scaling audit: what the simulations support

## Finding

The notes’ dimensional argument for `alpha + beta = 1` is incorrect, and the
simulations do not establish that sum as a universal collapse law. They do
show a narrow, retuned intermediate field window where a collapse-time-free
estimate is near one. At higher fields, independent resolution runs instead
give a stable sum near `3/2`. The mechanism behind the intermediate window
remains unresolved.

The issue in the notes is the gradient dimension. With `[phi] = T^-1` and
`[r] = T`, one has `[grad(phi)] = T^-2`, not `T^-1`. In the ansatz

```text
phi = tau^(-alpha) * Phi(r / tau^beta)
```

the gradient scales as `tau^(-alpha-beta)`; matching its dimension would give
a sum of 2, not 1. More generally, balancing the leading time, Laplacian,
and cubic terms of the simulated equation gives `alpha = beta = 1`, hence
sum 2. That dimensional argument therefore cannot derive the Module 5
relation.

## The strongest near-one evidence

To avoid choosing a collapse time `tc`, let `u = abs(phi(0,t))` and
`v = abs(du/dt)`. Define the measured slopes

```text
s = d(log(v)) / d(log(u))
q = -d(log(R50)) / d(log(u))
```

If these local power laws hold,

```text
u    ~ (tc - t)^(-alpha)
R50  ~ (tc - t)^beta
```

then

```text
alpha      = 1 / (s - 1)
beta       = alpha * q
alpha+beta = (1 + q) / (s - 1)
```

Thus sum one requires `s = 2 + q`. This estimate uses the evolved central
velocity and the measured radius at fixed central-field thresholds; it does
not fit `tc`.

| Run and field band | s | q | inferred alpha+beta |
|---|---:|---:|---:|
| N=1000, +0.030000%, 12–30 | 2.501 | 0.503 | 1.001 |
| N=2000, +0.030059063%, 12–30 | 2.524 | 0.502 | 0.986 |
| N=2000, +0.030054297%, 12–30 | 2.307 | 0.503 | 1.150 |

The first two runs agree across resolution; a smaller timestep repeats the
N=2000 value. Other nearby runs give 1.03–1.04. But the third row shows
strong amplitude sensitivity, so this is a narrow, run-specific intermediate
fit—not a precision exponent or proof of a universal law.

Collapse-time choice explains one earlier misleading fit. For the N=2000
+0.030059063% run, the `tau = 0.1–0.3` fit (fields 8.6–16.9) gives 1.008,
but the no-`tc` estimate on those same events is only 0.684. In the 12–30
band it is 0.986, and reconstructed local collapse times scatter less. So
the lower-field interval is still curved/transitional; its near-one result
depends on the selected `tc`. Fit `R^2` values are high, but are fit-quality
measures, not uncertainty estimates.

Changing the radius threshold from 50% to 25% or 75% gives sums 0.97–1.01
for the two near-one runs, with radius-regression `R^2 > 0.9999`. This makes
it unlikely that the result is solely an artifact of the exact half-height
choice, though the amplitude and interval sensitivity remains. See the
[collapse-time-free scaling plot](regime_scan/results/tc_free_validation/tc_free_scaling.png)
and [fit table](regime_scan/results/tc_free_validation/tc_free_fits.csv).

## Profile and late-stage checks

At `abs(phi_0) = 12, 20, 30`, three near-transition runs have almost
overlapping profiles after scaling radius by `R50` and field by `abs(phi_0)`.
The best fit to `1/(1+rho^m)` has `m = 2.005–2.105`; fixed
`1/(1+rho^2)` has RMSE 0.0005–0.009 over `0 <= rho <= 2`, compared with
0.063–0.071 for a half-height-normalized Gaussian. This is a useful
descriptive shape match, not independent radius/exponent evidence: the
normalization forces agreement at the center and at `rho = r/R50 = 1`, and
the interpolated radial samples are correlated. It also does not show that
the shape persists asymptotically. See the
[profile plot](regime_scan/results/profile_window_check/profile_evolution.png)
and [shape-fit data](regime_scan/results/profile_window_check/profile_shape_fits.csv).

At higher fields, shared-`tc` fits over `abs(phi_0) = 160–1280` give
`alpha ≈ 1`, `beta ≈ 0.5`, and sum 1.492–1.496 for N=500–4000. The late-stage
local ODE `phi_tt = phi^3` predicts `alpha = 1`; a smooth radial variation of
local blow-up time,

```text
T_loc(r) = tc + kappa*r^2 + ...
phi(r,t)/phi(0,t) ≈ 1 / (1 + kappa*r^2/(tc-t))
```

predicts `R50 proportional to (tc-t)^(1/2)` and sum 3/2. The high-field
profiles and PDE-term ratios support this late-stage interpretation.

The intermediate pair `(alpha, beta) ≈ (2/3, 1/3)` is not explained by the
same simple balance. Its time, Laplacian, and cubic terms scale as
`tau^(-8/3)`, `tau^(-4/3)`, and `tau^(-2)`, respectively; avoiding an
unmatched leading time term would require a singular profile in the
single-profile similarity ansatz. Numerically, the Laplacian is still
material in the 12–30 band (core norm ratio about 0.35), and the central
acceleration trend also departs from the constant-slope prediction. The
evidence points to a finite-window crossover or more complicated profile
evolution, not a resolved analytic mechanism.

## Scope and reproduction

“Transition” here is operational: reaching `abs(phi(0)) = 12` before
`t = 100`, not proof of finite-time blow-up or an exact critical amplitude.
At N=2000, +0.0300542224121% does not reach 12, while +0.030054296875%
reaches 1280; this extremely narrow finite-horizon bracket should not be
called a mathematical `A_critical`. The radius is `R50`, not the
energy-weighted `sqrt(<r^2>)`. Thresholds are interpolated between RK4
endpoints, and no formal confidence intervals are calculated.

Reproduce the profile-window check from the repository root with:

```bash
python3 module4/Scale_invariance/scaling_diagnostics/module5_audit/regime_scan.py --case 1000:0.03 --case 2000:0.0300590625 --case 2000:0.030054296875 --profile-level 12 --profile-level 20 --profile-level 30 --tag profile_window_check
```

The scripts, settings, event data, and fit tables are in this folder and its
`regime_scan/results/` subfolders. CSV outputs are locally available but
ignored by the project-wide Git rule.
