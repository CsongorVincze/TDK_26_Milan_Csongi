#!/usr/bin/env python3
"""Audit one nearby finite-horizon outcome transition at N=1000."""

import argparse
import csv
import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
HALF_STUDY_PATH = HERE.parent / "Fraction_of_central" / "half_radius" / "study.py"
N = 1000
A_CRIT_N500 = 1.202609700895846
T_MAX, DT_SAMPLE, PHI_CUTOFF = 100.0, 0.01, 12.0
INITIAL_OFFSETS_PERCENT = (0.02875, 0.03)


def load_study():
    spec = importlib.util.spec_from_file_location("half_radius_critical_audit", HALF_STUDY_PATH)
    study = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = study
    spec.loader.exec_module(study)
    return study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tolerance", type=float, default=1e-8)
    args = parser.parse_args()
    if args.tolerance <= 0:
        parser.error("Tolerance must be positive.")

    study = load_study()
    trials = []

    def classify(amplitude):
        rows, crossed, failed, event_time, steps = study.evolve(
            amplitude, N, T_MAX, DT_SAMPLE, PHI_CUTOFF, record=True,
        )
        if failed:
            raise RuntimeError(f"Numerical failure at A={amplitude:.12f}")
        row = {
            "N": N, "amplitude": amplitude,
            "offset_from_N500_Acrit_percent": 100.0*(amplitude/A_CRIT_N500-1.0),
            "horizon": T_MAX, "collapse_cutoff": PHI_CUTOFF,
            "collapsed_by_horizon": int(crossed),
            "event_time": event_time if crossed else "",
            "max_abs_phi_center": max(row["abs_phi_center"] for row in rows),
            "steps": steps,
        }
        trials.append(row)
        print(f"A={amplitude:.12f}: "
              f"{'cutoff reached' if crossed else 'no cutoff by horizon'}"
              + (f" at t={event_time:.6f}" if crossed else ""), flush=True)
        return bool(crossed)

    lo = A_CRIT_N500*(1.0+INITIAL_OFFSETS_PERCENT[0]/100.0)
    hi = A_CRIT_N500*(1.0+INITIAL_OFFSETS_PERCENT[1]/100.0)
    if classify(lo) or not classify(hi):
        raise RuntimeError("The saved-run bracket did not straddle the N=1000 cutoff.")

    while hi-lo > args.tolerance:
        mid = 0.5*(lo+hi)
        if classify(mid):
            hi = mid
        else:
            lo = mid

    midpoint = 0.5*(lo+hi)
    summary = [{
        "N": N, "horizon": T_MAX, "collapse_cutoff": PHI_CUTOFF,
        "tolerance": args.tolerance, "transition_lower_noncollapsing": lo,
        "transition_upper_collapsing": hi, "transition_midpoint": midpoint,
        "bracket_width": hi-lo,
        "midpoint_offset_from_N500_reference_percent": 100.0*(midpoint/A_CRIT_N500-1.0),
        "interpretation": "local transition only; other nearby amplitudes also collapse",
        "cfl": 0.4, "r_max": 25.0, "stretch": 5.0,
        "sponge_fraction": 0.8, "sponge_strength": 7.0,
        "sponge_power": 3.0, "R0": 2.0, "initial_velocity": 0.0,
    }]
    RESULTS.mkdir(parents=True, exist_ok=True)
    for path, rows in ((RESULTS/"local_transition_search_N1000.csv", trials),
                       (RESULTS/"local_transition_summary_N1000.csv", summary)):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    print(f"N=1000 local finite-horizon transition: [{lo:.12f}, {hi:.12f}] "
          f"(midpoint {midpoint:.12f})")


if __name__ == "__main__":
    main()
