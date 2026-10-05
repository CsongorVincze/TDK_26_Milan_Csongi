# Module 4 collapse study

This folder contains a reproducible collapse-diagnostics run, its CSV data, plots, and the findings report.

- Read the conclusions in [report.md](report.md).
- Run `python3 module4/beta_collapse_study/run_study.py` from the project root to regenerate the simulations and outputs. Matplotlib uses a non-interactive backend, so no display is needed.
- Simulation time series and analysis products are in [results](results/).

The script reuses the Module 4 solver and saved critical-amplitude estimate. It does not alter the field equation or the main simulation program.
