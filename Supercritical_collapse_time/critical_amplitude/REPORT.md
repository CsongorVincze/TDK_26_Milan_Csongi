# Critical Gaussian amplitude: convergence and numerical error

## Result

For the positive-field collapse onset near A=1.20261 in the 3D Gaussian family,

**Ac = 1.202609733202 ± 5.0e-12**, dimensionless.

This is a conservative **model-based numerical uncertainty estimate**. It is
not a rigorous enclosure of the infinite-domain PDE threshold and not a
statistical confidence interval. The estimated band is
[1.202609733197000, 1.202609733207000]. It includes measured discretization and
sensitivity allowances, under the convergence assumptions described below.

The initial data are phi(r,0)=A exp(-r²/8), pi(r,0)=0 (R0=2). The PDE is
phi_tt=phi_rr+(2/r)phi_r-phi+phi³, with the original outer sponge tested
separately. The operational classification starts with reaching |phi(0,t)|=12
by t=100. Its stability is tested to t=200 and cutoff 96.

## An important correction to the earlier reference

The old N=500 reference, 1.2026097002647, was a **later negative-field
collapse-window edge**. It was not the first onset in its neighborhood.
The new amplitude map finds positive-field collapse already near
A=1.2026055 on that grid, followed by a non-crossing interval and then the
old negative-field branch. Thus the previous N=500/N=1000 comparison mixed
different boundaries and could not be used as a spatial error estimate.

Here, the positive-field onset is selected and tracked consistently. Its
N=500 second-order value is 1.202605393181419. Every selected
collapsing bracket endpoint is verified to reach positive field 12, and
branch-validation runs at raw offset 1e-7 reach positive field 48. The scan
establishes the onset within the mapped neighborhood [1.202602,1.202614]; it
does not prove that no other collapse islands exist at arbitrary lower A.
As there are multiple windows, Ac must always identify the chosen boundary.

## How the determination was calculated

The stretched radial map and inner nodes are fixed. Nested grids use
N=500,999,1997,3993, giving an exact factor of two in computational spacing
h=1/(N-1). The time step scales with h at CFL=0.4. Two spatial methods are
compared: the original second-order stencil and a fourth-order stencil
including an even-origin stencil. RK4, potential, Gaussian profile, and
boundary convention are otherwise shared.

The new stencils were first tested against the analytic radial Laplacian
of a Gaussian. Their observed interior convergence orders are approximately
2.0001 and 3.9993. The accelerated original stencil also reproduces the
original Python solver's event times within floating-point roundoff on two
independent comparison runs. The fourth-order stencil falls back to second
order next to the outer boundary; actual threshold convergence and the
extended-domain tests verify that this does not dominate this calculation.

Each threshold has opposite classified bracket endpoints. Bisection search
half-widths are kept separate from spatial error:

Second-order spatial discretization:

| N | Critical midpoint | Search half-width |
| ---: | ---: | ---: |
| 500 | 1.202605393181419 | 5.722e-13 |
| 999 | 1.202608648163605 | 7.629e-13 |
| 1997 | 1.202609461940147 | 6.020e-13 |
| 3993 | 1.202609665387076 | 9.314e-13 |

Fourth-order spatial discretization:

| N | Critical midpoint | Search half-width |
| ---: | ---: | ---: |
| 500 | 1.202609733763905 | 1.910e-14 |
| 999 | 1.202609733237439 | 1.898e-14 |
| 1997 | 1.202609733204517 | 1.898e-14 |
| 3993 | 1.202609733202442 | 2.442e-14 |


## Richardson extrapolation and spatial uncertainty

For one tracked boundary, the tested leading-error model is

    Ac(h) = Ac(0) + C h^p + higher-order terms.

For three grids h,h/2,h/4, let d0=Ac(h)-Ac(h/2) and
d1=Ac(h/2)-Ac(h/4). Then

    p = log2(d0/d1)
    Ac(0)_estimate = Ac(h/4) - d1²/(d0-d1).

Observed orders from successive triplets:

- Second-order scheme: 1.999945, 1.999980.
- Fourth-order scheme: 3.999251, 3.987675.

The latest second-order extrapolation is 1.202609733203961.
The latest fourth-order extrapolation is 1.202609733202302.
Their difference is 1.659e-12.
The fourth-order estimate is the central reported value before rounding.

The standard absolute Grid Convergence Index allowance on the fine grid is

    U_GCI = 1.25 |Ac(0)_estimate - Ac(h_fine)|.

The safety factor 1.25 and observed-order approach follow the
[NASA/NPARC convergence guidance](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html).
This allowance is an engineering numerical estimate, not a theorem proving
that the exact PDE threshold lies inside it.

For the reported spatial allowance, take the maximum of the fourth-order
GCI, the discrepancy between independent-scheme extrapolations (including
their propagated search brackets), and the change between successive
fourth-order extrapolations (also including search brackets). This prevents
reporting a small GCI when the independent checks disagree.

Bracket propagation uses outward-rounded interval arithmetic through the
Richardson expression. It bounds the effect of the supplied search intervals
**within that algebraic extrapolation model**. It does not bound omitted
higher-order terms or convert a floating-point solver into a validated PDE
solver. Those are assessed empirically by the extra grid and independent
stencil.

## Other error sources and the combined estimate

At the N=999 fourth-order grid, independent local searches tested CFL=0.2
and 0.1; horizon 200; cutoffs 48 and 96; extended outer radii at least 50
and 75; sponge strength 3.5 and zero; and extended floating-point precision
(long double). A joint test also combines radius at least 75, horizon 200,
cutoff 96, and CFL=0.2. Domain extension adds outer grid nodes, preserving every
inner grid point and the inner time-step scale. The sponge is moved outward
for larger domains, not simply rescaled on an altered core grid.

No systematic shift is resolved beyond the search intervals in these tests.
For each setting, its allowance is the absolute midpoint shift plus both
search half-widths. For each category, use the largest allowance among its
tests and retain that N=999 value without reducing it for finer grids.
The deterministic contributions are added, not combined as independent
random standard deviations:

| Contribution | Absolute amplitude allowance |
| --- | ---: |
| Spatial convergence and independent-method comparison | 3.966e-12 |
| Time-step sensitivity | 5.007e-14 |
| Outer domain and sponge | 5.806e-14 |
| Horizon and collapse cutoff | 5.806e-14 |
| Arithmetic precision | 5.007e-14 |
| Propagated search brackets | 3.086e-14 |
| Sum before reporting/rounding | 4.213e-12 |
| Reported allowance | 5.0e-12 |


The allowance is rounded upward for presentation, including rounding of the
reported center. Extended arithmetic agrees at the tested grid, but this
is not a formal roundoff bound. The budget also assumes the observed spatial
convergence persists and that the tested domain/horizon/cutoff variations
are representative. A rigorous mathematical guarantee would require a
validated PDE error analysis or interval-in-time/space solver; that has not
been performed here.

## How close can a supercritical run be labeled reliably?

For a fixed known starting amplitude A and an estimated threshold band
[Ac-U,Ac+U], the corresponding distance band is

    delta_min = A/(Ac+U)-1
    delta_max = A/(Ac-U)-1.

Thus an accurate Ac must be much more accurate than the imposed amplitude
excess. The following intervals use the reported threshold band and include
floating-point amplitude construction in the nominal delta. The final column
is a **fraction**, not a percentage:

| Target raw delta | Starting amplitude | Estimated interval for delta | Relative uncertainty in distance |
| ---: | ---: | ---: | ---: |
| 1e-04 | 1.202729994175320 | [9.9999996e-05, 1.0000000e-04] | 4.158e-08 |
| 1e-06 | 1.202610935811733 | [9.9999584e-07, 1.0000042e-06] | 4.158e-06 |
| 1e-08 | 1.202609745228097 | [9.9958424e-09, 1.0004158e-08] | 4.158e-04 |
| 1e-10 | 1.202609733322261 | [9.5842223e-11, 1.0415757e-10] | 4.158e-02 |


To limit uncertainty to 0.01 of the stated distance, the minimum raw offset
implied by this threshold allowance is approximately
4.158e-10. This concerns the
threshold reference only; it is not a guarantee that the evolution or its
blow-up time is resolved at that distance.

In particular, do not use the new continuum Ac with the old N=500
second-order solver and assume delta=1e-8 is realized accurately. That grid's
same-branch threshold differs from the continuum estimate by about 4.34e-6.
Even N=3993 with the second-order stencil has a threshold shift about
6.78e-8. The fourth-order solver's grid thresholds are far closer to the
continuum estimate. A future near-critical collapse study should use the
fourth-order solver, compare fixed-A trajectories on multiple grids, and
separately verify the late collapsing core and blow-up-time estimator.

## Check of one specified supercritical amplitude

The same amplitude A=1.2026097452283997, nominal raw delta=1e-8 relative to the unrounded continuum estimate, was run with the fourth-order stencil:

| N | CFL | Time to 12 | Time to 24 | Time to 96 |
| ---: | ---: | ---: | ---: | ---: |
| 999 | 0.4 | 23.5562983104 | 23.6363892031 | 23.6857471754 |
| 1997 | 0.4 | 23.5558113793 | 23.6358940212 | 23.6852509625 |
| 3993 | 0.4 | 23.5557810454 | 23.6358631341 | 23.6852209505 |
| 1997 | 0.2 | 23.5558114148 | 23.6358940456 | 23.6852519538 |

All four runs reach positive field 96. These are cutoff times, **not** extrapolated singularity times. This test checks the same initial condition across grids; it does not retune the amplitude to each grid threshold. A blow-up-time estimate still needs a resolved late-field fit and its own convergence study.

## Reproduction and artifacts

From the project root:

```bash
MPLCONFIGDIR=/tmp/tdk-critical-mpl python3 Supercritical_collapse_time/critical_convergence/run.py all
```

Requires Python, NumPy, Matplotlib, and g++. Simulations are cached by source
hash and all numerical settings in trials.jsonl. Delete or move the results
folder only if a fresh uncached rerun is desired; ordinary reruns reuse the
saved computations. No commit or push is performed. All figures are English
PNG exports following TDK_plot_skill; no PDFs are generated.

Key outputs:

- determination.json: full result, observed orders, and uncertainty components.
- grid_convergence.csv, error_budget.csv, sensitivity_differences.csv.
- *_trials.csv and trials.jsonl: every actual classification.
- critical_convergence.png: independent convergence and extrapolation checks.
- branch_tracking.png: the collapse windows and the old branch mismatch.
- branch_validation.csv and branch_*.csv: same-branch trajectory checks.
- fixed_amplitude_validation.csv and fixed_delta1e-8_*.csv.
- usable_offsets.csv: raw-offset uncertainty propagation.
- manufactured_stencil_validation.csv and evolution_validation.csv.

The main uncertainty calculation and all empirical assumptions are exposed
in analyze.py; it does not attach an unsupported statistical confidence level.
