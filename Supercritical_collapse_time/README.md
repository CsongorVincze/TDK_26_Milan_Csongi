# Gaussian collapse windows

<!-- survey-summary -->
**Current result: 6 non-crossing windows and 7 collapse bands**, from 3783 primary amplitudes. All counts are finite-horizon numerical observations. Validated window interiors: 6/6.
<!-- /survey-summary -->

This folder contains the current amplitude-window search and the validated
critical-amplitude reference. Figures use English labels and raw
`delta = (A - Ac) / Ac` ratios. PNG only.

## Start here

- [Window search report](results/REPORT.md): counts, definitions, and limitations.
- [Window map](results/window_map.png): cutoff time versus log(raw amplitude excess).
- [Full amplitude range](results/amplitude_overview.png).
- [Individual gaps](results/gap_details.png).
- [Representative trajectories](results/representative_trajectories.png).
- [Window intervals](results/windows.csv) and [all primary samples](results/amplitude_samples.csv).
- [Resolution checks](results/resolution_checks.csv).

The new survey spans **A = 0.2 to 4 Ac**, with dense positive raw offsets
**1e-8 to 3**, subcritical controls, and adaptive sampling of changes in
collapse sign, zero-crossing count, and no-crossing gaps. Runs reach field 48
or t=100. Plotted t12 values preserve the original field-12 cutoff-time
convention; they are not extrapolated singularity times.

## Critical-amplitude reference

`Ac = 1.202609733202 +/- 5e-12` is the prior continuum estimate for the
positive-field collapse onset near A=1.20261. Its uncertainty is a numerical
estimate, not a rigorous PDE bound. Read [the reference report](critical_amplitude/REPORT.md)
and [its folder guide](critical_amplitude/README.md). The broad survey uses
N=999 with the fourth-order stencil; representative windows are checked on
N=1997 and N=3993. Window-edge positions remain resolution dependent.

## Folder layout

- `results/`: current figures, readable report, and compact CSV tables.
- `critical_amplitude/`: the accepted critical-amplitude report and summary.
- `src/`: the active solver, search/plot script, and TDK plot style.
- `data/`: complete simulation cache and compressed representative trajectories.
- `archive/`: one verified archive of all superseded work, with checksums.
- `.cache/`: generated executable and progress logs; not needed for inspecting results.

The old studies were archived only after every original file passed a
SHA-256 content check. They no longer occupy the active results folder.
The archive also preserves the earlier negative-field window reference and
its old trajectories; those should not be mixed with the current onset.

## Reproduce or refresh

Requires g++, Python, NumPy, and Matplotlib. From the project root:

```bash
python3 Supercritical_collapse_time/src/window_search.py
```

The search resumes from `data/runs.jsonl`. For figures and tables only:

```bash
python3 Supercritical_collapse_time/src/window_search.py --stage plots
```

Individual stages are `survey`, `refine`, and `validate`. Validation uses
fixed starting amplitudes across resolutions; non-crossing representatives
also get a larger-domain t=200 check. No PDFs, commits, or pushes are created.
