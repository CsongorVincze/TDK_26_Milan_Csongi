# Critical-amplitude precision and evolution through t = 100

All new offset axes use the raw dimensionless ratio `delta=(A-Ac)/Ac`, on a
logarithmic axis. Times are first crossings of `|phi(0,t)|=12` (`t12`), not
fits to a singularity time `tc`. Every run ends at t=100 or its first cutoff.
No PDFs are generated.

## Precision of Ac

The reference in the original 502-amplitude sweep was
`Ac_old=1.2026097002526517`, with bracket width `9.895306796e-11`.
Its half-width corresponds to relative offset uncertainty `4.114e-11`.
Consequently, at delta=1e-9 the bracket uncertainty alone is about 0.041
of the imposed excess; the nominal offset should not be treated as exact.
The new N=500 midpoint is only `1.2052e-11` above that old reference.

Fresh local bisections, with all classifications evolved through t=100:

| Grid N | CFL | Bisection-selected midpoint | Bracket width |
| ---: | ---: | ---: | ---: |
| 500 | 0.4 | 1.2026097002647038 | 5.9952e-15 |
| 500 | 0.2 | 1.2026097002648428 | 6.2172e-15 |
| 500 | 0.1 | 1.2026097002648548 | 5.9952e-15 |
| 1000 | 0.4 | 1.2026086500147362 | 7.9162e-10 |

For N=500, CFL=0.4, the refined endpoints are:

- Non-crossing: `1.2026097002647007`.
- Crossing: `1.2026097002647067`.

This tiny bracket is **search precision within a deterministic double-precision
implementation**, not an error bound on Ac. The largest tested time-step shift
is `1.51e-13` in A, or `1.26e-13` as a raw offset, much larger than the refined
bracket half-width. Offsets comparable to that scale are not time-step robust.
Roundoff and changes in arithmetic ordering are not bounded by the bracket.
An independent check with the original Python solver classifies amplitudes
`Ac_refined +/- 1e-12` consistently (lower does not cross, upper does).

Increasing N to 1000 selects a transition about `1.05025e-6` lower in A,
a relative shift of about `8.73e-7`. This is a resolution-sensitivity diagnostic,
not a continuum-error estimate. These are bisection-selected local transitions;
non-monotone outcomes forbid assuming a unique global threshold. **Six-decimal
continuum accuracy has not been established.** The much tighter N=500 bracket
cannot remove that spatial-discretization uncertainty.

## Windows and the new close-offset scan

The previous percentage axis must be divided by 100 to obtain raw delta.
For example, the former tick 1e-3 meant raw delta=1e-5.
The first previously sampled non-crossing band is near raw delta=3.2e-5;
the second is near 2.4e-3. Switching to raw units does not move the data
relative to each other; it only rescales their coordinate labels.

202 fresh simulations use the refined N=500 midpoint and requested offsets
zero plus 201 logarithmically spaced values from 1e-13 to 1e-3. The CSV saves
both requested offsets and actual floating-point `(A-Ac)/Ac` values.
One new positive-offset run, near delta=3.16228e-5, misses the cutoff by t=100.
The zero-offset midpoint happens to cross at about t=24.72; the first positive
sample crosses near t=23.00. Classification at the last few floating-point
amplitude increments should not be overinterpreted.

The updated near-critical plot also includes existing amplitudes in this
range, re-expressed relative to the refined midpoint. Thus it preserves the
densely sampled window edges rather than relying on the new logarithmic
sampling to resolve each narrow window. The old CSVs retain their original
reference for provenance; their coordinates are recomputed when overlaid.
Below about 3e-5 the sampled curve remains on one negative-field collapse
branch. This finite scan does not prove that arbitrarily narrow windows are
absent, and it does not establish continuum behavior. The visible branch
jumps are accompanied by changes in the final field sign.

## Evolution through t = 100

All 87 original runs with end time greater than 20 were freshly rerun:

- 42 reached the field cutoff after t=20.
- 45 remained below the cutoff through t=100.

The rerun classifications and cutoff times agree with their saved originals.
The surviving central fields have substantial early transients followed by
smaller oscillations, not delayed cutoff crossings before t=100. Symmetric-log
field scales make these late oscillations visible without clipping the early
transients. This alone does not prove a persistent localized oscillon.
Late-time structure is measured in the existing Rmax=25 domain and sponge;
its dependence on the outer domain and boundary has not been tested here.

Six representative simulations additionally save the radial field every
0.2 time units, with a final snapshot at t=100 or the cutoff. The radius array
is saved separately. Space-time panels display r<=8 and use a symmetric-log
color normalization to expose weak late fields. Panels are blank after a
cutoff because post-cutoff evolution was not simulated.

## Figures and raw data

- `full_evolution_t100.png`: all original trajectories, plus t=20..100 detail.
- `surviving_oscillons_t100.png`: three individual non-crossing central fields.
- `spatial_evolution_t100.png`: six representative radial field histories.
- `nearcritical_precision.png`: refined close-offset scan and precision scales.
- `precision_summary.csv`, `search_N*_cfl*.csv`: every precision-search trial.
- `python_reference_check.csv`: original Python solver comparisons.
- `survivors_t100.csv`: rerun outcomes and trajectory paths.
- `nearcritical_sweep.csv`: the 202 new near-critical runs.
- `spatial_runs.csv`, `spatial_radius.csv`, `fields_*.csv`: spatial outputs.
- `trajectories/*.csv`: sampled central fields for all new runs.

From the project root, rerun the complete audit with:

```bash
MPLCONFIGDIR=/tmp/supercritical-mpl python3 Supercritical_collapse_time/precision_longtime.py
```

Regenerate its figures without simulation:

```bash
MPLCONFIGDIR=/tmp/supercritical-mpl python3 Supercritical_collapse_time/precision_longtime.py --plot-only
```
