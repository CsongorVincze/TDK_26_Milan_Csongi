# Module 4 collapse and beta study

## Main finding

The core really does contract as the central field grows. The earlier signed-energy radius was the misleading diagnostic: the energy density changes sign, so positive and negative contributions nearly cancel in its moment. Its apparent expansion is not evidence that the field core expands.

The most reliable result so far is the direct width-versus-central-field scaling, `R_core ~ |phi(0,t)|^(-q)`, with `q = beta/gamma` close to 0.51–0.53. A time exponent `beta` cannot yet be quoted as a unique final value because the central-field exponent `gamma` and extrapolated collapse time `t_c` are still sensitive to how close to collapse we fit. In the extended runs, the two reasonable time-fit choices give `beta` about 0.45–0.49 by the `|phi(0)|=80` cutoff. A value near 0.5 is plausible, but it depends on the asymptotic central-field law and needs another convergence check before being treated as settled.

## What was tested

The new [runner](run_study.py) reuses the Module 4 radial RK4 solver and the saved N=500 critical-amplitude estimate, `Acrit = 1.2026097009`. Its main near-critical N=500 sweep uses amplitudes 0.01%, 0.05%, and 0.1% above that estimate. A second set of runs holds the smallest supercritical amplitude fixed and continues it to central-field cutoffs 30, 40, 60, and 80; one cutoff-80 continuation is repeated at N=1000. A fixed-amplitude resolution check also compares N=250, 500, and 1000 at cutoff 20.

The inherited setup is `Rmax=25`, initial width `R0=2`, grid stretch 5, core radius 8, and output spacing `dt_sample=0.01`. The time step is adaptively limited by the solver. The field equation and potential are unchanged; the energy density used for diagnostics is `rho = 0.5*pi^2 + 0.5*(dphi/dr)^2 + 0.5*phi^2 - 0.25*phi^4`.

We compare several definitions of core size:

- Signed energy RMS radius: `sqrt( integral(rho*r^4 dr) / integral(rho*r^2 dr) )`, evaluated with core cutoffs R=4, 8, and 12. This is the square root of the signed `<r^2>` observable; the raw `<r^2>` is also saved.
- Absolute- and positive-energy RMS radii, using `abs(rho)` and `max(rho, 0)` as weights inside R=8.
- Gradient-energy RMS radius, weighted by `(dphi/dr)^2` inside R=8.
- Field-profile widths `r50` and `r20`, where `abs(phi(r,t))` first falls to 50% or 20% of its central value.

For a radius-like quantity, the requested log-log fit is

```text
R(t) = C * (tc - t)^beta
log10(R) = beta * log10(tc - t) + constant
```

For raw `<r^2>`, the fitted log-log slope is twice the radius exponent, so the reported radius-equivalent beta is half that slope. We varied the remaining-time fit window (`tc-t` from 0.1–0.3, 0.2–0.6, 0.4–1.0, and 0.8–1.6) and compared multiple ways to estimate `tc`.

## Core contraction and why the old moment looked wrong

For the N=500 run only 0.01% above the saved critical estimate:

| Central field | Time | Signed RMS R8 | Absolute-energy RMS R8 | Gradient RMS R8 | Field r50 | Signed/absolute moment ratio |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 12.06 | 15.1500 | 2.9083 | 1.3474 | 0.9244 | 0.4356 | 0.2097 |
| 20 | 15.2032 | 3.0225 | 1.0029 | 0.7167 | 0.3382 | 0.1069 |
| 40 | 15.2411 | 3.2072 | 0.6675 | 0.5111 | 0.2397 | 0.0418 |
| 80 | 15.2594 | 3.3635 | 0.4473 | 0.3656 | 0.1698 | 0.0170 |

From central field 20 to 80, the profile half-width shrinks by about 50%, and the gradient and absolute-energy radii also shrink. In contrast, the signed RMS radius increases and is non-monotonic at higher cutoffs. At `|phi(0)|=80`, the signed moment is only about 1.7% of the corresponding absolute-weighted denominator. That severe cancellation makes its beta fit ill-conditioned.

The continued runs stop at the selected field threshold; they do not evolve past 80. There is no `<r>` statistic in this study: the field at the origin is `phi(0,t)`, while the radius observables are measured separately.

## Beta estimates and method dependence

All values below use the same near-critical N=500 trajectory and the common fit window `0.2 <= tc-t <= 0.6` unless a different cutoff is specified.

| Observable | `tc` method | Beta | Fit R-squared | Reading |
| --- | --- | ---: | ---: | --- |
| Signed-energy RMS R8 | `1/abs(phi0)` linear fit, last 40 tail samples (`gamma=1`) | 0.002 | 0.046 | Fails as a collapse-width fit |
| Absolute-energy RMS R8 | Same fixed-`gamma` `tc` | 0.438 | 0.987 | Contracts, but weight is sensitive to negative energy density |
| Gradient-energy RMS R8 | Same fixed-`gamma` `tc` | 0.519 | 0.998 | Strong radius-like scaling |
| Field `r50` | Same fixed-`gamma` `tc` | 0.518 | 1.000 | Strongest simple profile-width fit |
| Field `r20` | Same fixed-`gamma` `tc` | 0.524 | 1.000 | Consistent with `r50` in this window |
| Gradient-energy RMS R8 | Jointly fitted `tc` and `gamma`, 40 tail samples | 0.412 | 1.000 | Strong fit, with a different central-field model |
| Field `r50` | Jointly fitted `tc` and `gamma`, 40 tail samples | 0.430 | 1.000 | Strong fit, with a different central-field model |

The fixed-`gamma` method assumes `abs(phi0) ~ (tc-t)^(-1)` and fits `1/abs(phi0)` linearly in time. The alternative jointly fits

```text
abs(phi0) ~ (tc - t)^(-gamma)
```

without forcing `gamma=1`. At cutoff 20 it finds `gamma=0.753` and `tc=15.2560`, compared with `tc=15.3210` under the fixed-`gamma` method. The same `r50` time fit consequently gives beta 0.430 or 0.518. Both radius fits look straight over the chosen short window; the difference comes mainly from `tc`/`gamma`, not from a visible failure of the radius data to follow a line.

The independently fitted `tc` is also window-sensitive. On the cutoff-20 run, fitting `1/abs(phi0)` in the field band 2–4 gives `tc=15.672`, while the 4–8 band gives 15.369 and the last-40-tail rule gives 15.321. The highest bands contain only 6–7 sampled points; their impressive R-squared values do not remove their extrapolation uncertainty. This is why a high R-squared alone is not enough to establish beta.

The direct width-versus-field fit avoids `tc` altogether. Over the last 40 samples with `abs(phi0) >= 4`, the N=500 `r50` exponent `q=beta/gamma` is 0.5299, 0.5339, and 0.5307 for the three amplitudes. It is 0.5312, 0.5299, and 0.5307 for the N=250, 500, and 1000 runs at fixed amplitude. The fit R-squared is about 0.999. This ratio is much more stable than beta by itself.

Continuing the smallest supercritical run changes the joint central-field exponent from `gamma=0.753` at cutoff 20 to 0.853 at cutoff 80, while the direct `r50` ratio changes only from `q=0.530` to 0.511. The joint-fit time exponent for `r50` moves from beta 0.430 to 0.448; the fixed-`gamma=1` estimate moves from 0.518 to 0.487. At cutoff 80, N=500 and N=1000 agree closely: `r50=0.16982` versus 0.16917, with joint-fit beta 0.4482 for both. This supports a genuine, resolved contraction, but the drift in `gamma` means the asymptotic beta is not yet pinned down.

## Bottom line and caveats

1. Core contraction is visible in field widths, gradient energy, and positive/absolute energy diagnostics; it persists through `abs(phi0)=80`.
2. Do not use the existing signed-energy `<r^2>` moment as the core-width beta here. It is dominated by cancellations and can increase while the central profile contracts.
3. The robust measured scaling is `beta/gamma ~= 0.51–0.53`. If the central field ultimately has `gamma=1`, beta would be near 0.5. With the currently fitted `gamma~=0.85` at the highest tested cutoff, beta is nearer 0.45. The honest current answer is therefore a provisional beta range of about 0.45–0.50, not a definitive universal number.
4. The saved `Acrit` came from the existing N=500 threshold search. The N=250/500/1000 check holds that same amplitude fixed; it does not independently retune the critical amplitude at each resolution. The cutoff-80 convergence check is encouraging but only compares N=500 and 1000.
5. These simulations stop at a chosen finite central-field threshold. They show increasingly concentrated collapse-like dynamics, but do not prove the continuum PDE's singularity structure or establish behavior beyond `abs(phi0)=80`.

## Files and rerunning

Each trajectory CSV in [results](results/) contains the time, central field, signed `<r^2>` and RMS radii, positive/absolute/gradient radii, field widths, cancellation ratio, and energy diagnostics. The main analysis tables are `run_summary.csv`, `tc_sensitivity.csv`, `beta_fit_sensitivity.csv`, and `amplitude_scaling.csv`.

Key figures:

- [Core-observable evolution](results/observable_evolution.png)
- [Direct width versus central-field scaling](results/amplitude_scaling.png)
- [Continuation from cutoff 20 to 80](results/cutoff_continuation.png)
- [Explicit time-based log-log fits](results/representative_loglog_fits.png) and [free-gamma alternative](results/representative_loglog_fits_joint_tc.png)
- [Beta sensitivity to fit window](results/beta_fit_sensitivity.png), [collapse-time estimator](results/tc_sensitivity.png), and [resolution](results/resolution_sensitivity.png)

From the project root, rerun with:

```bash
MPLCONFIGDIR=/tmp/mplconfig_module4 python3 module4/beta_collapse_study/run_study.py
```
