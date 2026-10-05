#!/usr/bin/env python3
"""Summarize N=1000 collapse outcomes across the saved amplitude scan."""

import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MANIFEST = (ROOT / "Fraction_of_central" / "half_radius" / "results" /
            "alpha_beta_scaling" / "tc_precision_run_manifest.csv")
RESULTS = HERE / "results"
LOCAL_TRANSITIONS = RESULTS / "local_transition_search_N1000.csv"
HORIZON = 100.0


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def first_crossing(rows, field, level):
    """Linearly interpolate the first crossing of an absolute central field."""
    for index, row in enumerate(rows):
        value = abs(float(row[field]))
        if value >= level:
            if index == 0:
                return float(row["time"])
            before = rows[index - 1]
            left = abs(float(before[field]))
            fraction = (level - left) / (value - left)
            return float(before["time"]) + fraction * (
                float(row["time"]) - float(before["time"]))
    return None


def main():
    rows = []
    manifest_rows = read_csv(MANIFEST)
    trajectory_dir = MANIFEST.parent.parent
    for item in manifest_rows:
        trajectory_path = trajectory_dir / item["half_trajectory_csv"]
        trajectory = read_csv(trajectory_path)
        values = [abs(float(row["abs_phi_center"])) for row in trajectory]
        rows.append({
            "N": 1000,
            "offset_from_N500_Acrit_percent": float(item["offset_percent"]),
            "amplitude": float(item["amplitude"]),
            "collapsed_by_t100_phi12": int(max(values) >= 12.0),
            "first_phi12_time": first_crossing(trajectory, "abs_phi_center", 12.0) or "",
            "first_phi20_time": first_crossing(trajectory, "abs_phi_center", 20.0) or "",
            "first_phi80_time": first_crossing(trajectory, "abs_phi_center", 80.0) or "",
            "max_abs_phi_center": max(values),
            "last_time": float(trajectory[-1]["time"]),
            "source": "saved amplitude scan",
        })

    if LOCAL_TRANSITIONS.exists():
        for item in read_csv(LOCAL_TRANSITIONS):
            offset = float(item["offset_from_N500_Acrit_percent"])
            amplitude = float(item["amplitude"])
            event_time = item["event_time"]
            replacement = {
                "N": 1000,
                "offset_from_N500_Acrit_percent": offset,
                "amplitude": amplitude,
                "collapsed_by_t100_phi12": int(item["collapsed_by_horizon"]),
                "first_phi12_time": event_time,
                "first_phi20_time": "",
                "first_phi80_time": "",
                "max_abs_phi_center": float(item["max_abs_phi_center"]),
                "last_time": HORIZON if not event_time else float(event_time),
                "source": "local recorded-step transition bisection",
            }
            match = next((index for index, old in enumerate(rows)
                          if math.isclose(old["amplitude"], amplitude, abs_tol=1e-13)), None)
            if match is None:
                rows.append(replacement)
            else:
                rows[match] = replacement

    rows.sort(key=lambda row: row["amplitude"])
    RESULTS.mkdir(parents=True, exist_ok=True)
    table_path = RESULTS / "n1000_finite_horizon_outcomes.csv"
    with table_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    sweep = [row for row in rows if row["source"] == "saved amplitude scan"]
    transition = [row for row in rows if row["source"] != "saved amplitude scan"]
    collapsed = [row for row in sweep if row["collapsed_by_t100_phi12"]]
    noncollapsed = [row for row in sweep if not row["collapsed_by_t100_phi12"]]
    t_collapsed = [row for row in transition if row["collapsed_by_t100_phi12"]]
    t_noncollapsed = [row for row in transition if not row["collapsed_by_t100_phi12"]]
    fig, (overview, zoom) = plt.subplots(1, 2, figsize=(12.4, 5.7),
                                         gridspec_kw={"width_ratios": [1.2, 1]})
    for ax, xlim in ((overview, (0.007, 0.103)), (zoom, (0.027, 0.032))):
        ax.scatter([r["offset_from_N500_Acrit_percent"] for r in collapsed],
                   [float(r["first_phi12_time"]) for r in collapsed],
                   color="#2563eb", marker="o", s=42,
                   label="Saved scan: crossed $|\\phi_0|=12$")
        ax.scatter([r["offset_from_N500_Acrit_percent"] for r in noncollapsed],
                   [r["last_time"] for r in noncollapsed],
                   color="#dc2626", marker="x", s=60, linewidth=2,
                   label="Saved scan: no crossing by $t=100$")
        if t_collapsed:
            ax.scatter([r["offset_from_N500_Acrit_percent"] for r in t_collapsed],
                       [float(r["first_phi12_time"]) for r in t_collapsed],
                       facecolors="none", edgecolors="#111827", marker="D", s=55,
                       label="Local bisection: crossing")
        if t_noncollapsed:
            ax.scatter([r["offset_from_N500_Acrit_percent"] for r in t_noncollapsed],
                       [r["last_time"] for r in t_noncollapsed],
                       color="#111827", marker="+", s=58,
                       label="Local bisection: no crossing")
        ax.set_xlim(*xlim)
        ax.set_ylim(0, 106)
        ax.set_xlabel("Excess over N=500 $A_{crit}$ reference (%)")
        ax.grid(alpha=0.25)
    overview.set_ylabel("First $|\\phi(0,t)|=12$ time (or horizon if not reached)")
    overview.set_title("Amplitude scan")
    zoom.set_title("Local transition detail")
    overview.legend(fontsize=7.5, loc="best")
    fig.suptitle("N=1000 finite-horizon outcomes are not monotone in amplitude")
    fig.tight_layout()
    fig.savefig(RESULTS / "n1000_finite_horizon_outcomes.png", dpi=180)
    plt.close(fig)

    failed = [row for row in rows if not row["collapsed_by_t100_phi12"]]
    print(f"Audited {len(rows)} unique N=1000 amplitudes; "
          f"{len(failed)} did not reach |phi0|=12 by t={HORIZON:g}.")
    print(f"Wrote {table_path} and {RESULTS / 'n1000_finite_horizon_outcomes.png'}")


if __name__ == "__main__":
    main()
