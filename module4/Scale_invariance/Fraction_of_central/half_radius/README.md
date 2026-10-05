# Half-central-field core radius

This study is self-contained: `study.py` includes its own radial-grid solver,
finite-time critical-amplitude search, supercritical runs, diagnostics, fits,
and plots. It imports nothing from the sibling `beta_collapse_study` or the older
Module 4 script.

The core radius is the first radius where the field magnitude has fallen to
half its central magnitude. See [report.md](report.md) for the results and
definitions. Generated CSVs and plots are written to [results](results/).

From the project root, run:

```bash
python3 module4/Scale_invariance/Fraction_of_central/half_radius/study.py
```

Useful options include `--critical-tmax`, `--critical-tolerance`, `--t-max`,
and `--sample-interval`. Use `--analyze-existing` to recalculate fits and
refresh simulation-diagnostic plots from the saved CSVs without rerunning the
PDE. Output paths are resolved from the script location, not the current
working directory.

The paired threshold-sensitivity search (only `N=1000` fits) loads both
standalone radius solvers and writes matched results to each method's own
`results/alpha_beta_scaling` folder:

```bash
python3 module4/Scale_invariance/tc_precision_search.py
```

The single-amplitude shared-`tc` baseline can be regenerated with:

```bash
python3 module4/Scale_invariance/shared_tc_alpha.py
```

That baseline now uses only the joint high-field-tail fit for `tc` and `alpha`;
the retired fixed-`alpha=1` path is not regenerated. Current multi-amplitude
plots and fit tables are documented in
[results/alpha_beta_scaling](results/alpha_beta_scaling/README.md). Both use the
same central-field data and `tc`; the radius normalization is called `C_R`, not
the notes' `alpha`.
