# Module 5 scaling audit

The [short findings report](report.md) checks the proposed
`alpha + beta = 1` relation against the same radial Klein–Gordon solver used
by the Module 4 half-radius study. It compares resolution, timestep, CFL, and
initial-amplitude choices.

Run the default experiment matrix from the repository root with:

```bash
python3 module4/Scale_invariance/scaling_diagnostics/module5_audit/run_audit.py
```

The default matrix reaches `|phi(0)|=1280` at N=500, 1000, 2000, and 4000;
checks a smaller CFL and field-dependent timestep at N=1000; and compares
`+0.01%`, `+0.05%`, and `+0.1%` amplitudes relative to the old N=500 finite-
horizon reference. This is not a fresh critical-amplitude search.

Results are written under `results/`: `event_diagnostics.csv` contains the
measured core/PDE quantities, `scaling_fits.csv` contains shared-`t_c` fits,
`term_ratio_fits.csv` summarizes PDE-term trends, `central_ode_fits.csv` checks
the linearized `1/|phi(0)|` relation and the evolved-velocity ratio, and
`run_settings.csv` records each run.
To additionally compare the near-core profile against rational
local-ODE and Gaussian shapes, run:

```bash
python3 module4/Scale_invariance/scaling_diagnostics/module5_audit/run_audit.py --profiles-only
```

That mode saves rescaled profile samples, model-fit errors, and a comparison
figure without replacing the main summary tables. The four main outputs
visualize PDE balance, core width, direct gradient scaling, and profile shape.
The output CSVs are ignored by the project-wide Git rule. Selected summaries
remain in the results folders locally; they would need deliberate force-adding
if you later approve including them in Git.

To scan how the exponents drift with collapse stage and amplitude, use the
separate regime-scan driver:

```bash
python3 module4/Scale_invariance/scaling_diagnostics/module5_audit/regime_scan.py
```

The default scans several N=1000 amplitudes. Add repeatable cases in the form
`N:offset_percent[:CFL[:field_step]]` to study a different resolution, and use
`--tag` to keep that batch in its own results subfolder. For example:

```bash
python3 module4/Scale_invariance/scaling_diagnostics/module5_audit/regime_scan.py --case 2000:0.0300590625 --tag n2000_finer2
```

To add profile snapshots at selected central-field thresholds, repeat
`--profile-level`; for example:

```bash
python3 module4/Scale_invariance/scaling_diagnostics/module5_audit/regime_scan.py --case 1000:0.03 --case 2000:0.0300590625 --profile-level 12 --profile-level 20 --profile-level 30 --tag profile_window_check
```

This writes `profile_snapshots.csv`, `profile_shape_fits.csv`, and
`profile_evolution.png` alongside the scaling outputs. The current comparison
shows nearly overlapping normalized profiles in that field band, closer to
`1/(1+(r/R50)^2)` than to a Gaussian; treat this as an intermediate-window
shape check, not an asymptotic result.

This scan fits one collapse time from central-field events 160–1280, then
uses it for both the central-field and R50 exponents. It also records a
collapse-time-free check. If `u = abs(phi(0))`,
`s = d(log(abs(u_dot))) / d(log(u))`, and
`q = -d(log(R50)) / d(log(u))`, a local power law implies
`alpha + beta = (1 + q) / (s - 1)`. The resulting `tc_free_fits.csv` and
`tc_free_scaling.png` show how that check compares across fixed field bands.
The fit table also compares radii at 25%, 50%, and 75% of the central-field
magnitude, to check sensitivity to the core-radius convention.
The scan records these alongside `tau = tc - t` fits, event data, fit-band
sensitivity, and plots. Its outputs stay under
`regime_scan/results/` (or the tagged subfolder) and do not overwrite the
main audit results.
