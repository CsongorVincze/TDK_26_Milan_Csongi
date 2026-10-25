# Critical-amplitude reference

The accepted reference is **Ac = 1.202609733202 +/- 5e-12**, a continuum
extrapolation for the positive-field collapse onset in this Gaussian family.
The uncertainty is a conservative model-based numerical estimate.

- `REPORT.md`: original detailed derivation and assumptions.
- `reference.json`: complete determination, grid results, and error budget.
- `error_budget.csv`, `grid_convergence.csv`: compact numerical tables.
- `convergence.png`: independent-method convergence figure.
- `previous_branch_map.png`: why the old N=500 window edge was a different boundary.

The original report describes a historical `critical_convergence/` source
layout. Its full scripts, raw trials and reproducibility data are preserved
inside `../archive/previous_studies_2026-10-25.tar.gz` under that path. The
current active solver is `../src/solver.cpp`; it adds zero-crossing counts
without changing the evolved equations or integration scheme.
