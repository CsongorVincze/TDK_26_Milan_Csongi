#!/usr/bin/env python3
"""Compare alpha+beta across tau windows using cached N=1000 trajectories."""

import csv
import importlib.util
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SEARCH_PATH = ROOT / "tc_precision_search.py"
RESULTS = HERE / "results"
OFFSETS = (0.01, 0.03125, 0.05, 0.1)
WINDOWS = ((0.05, 0.1), (0.1, 0.2), (0.1, 0.3), (0.2, 0.4))
CUTOFFS = (20, 24, 40, 60, 80)


def load_search():
    spec = importlib.util.spec_from_file_location("tc_precision_search_audit", SEARCH_PATH)
    search = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = search
    spec.loader.exec_module(search)
    return search


def main():
    search = load_search()
    output = []
    for offset in OFFSETS:
        half_path, gaussian_path = search.trajectory_paths(offset)
        half_rows = search.read_csv(half_path)
        gaussian_rows = search.read_csv(gaussian_path)
        for cutoff in CUTOFFS:
            half_cut, crossing_time = search.interpolate_crossing(half_rows, cutoff)
            gaussian_cut, _ = search.interpolate_crossing(gaussian_rows, cutoff)
            if half_cut is None or gaussian_cut is None:
                continue
            tc_fit = search.threshold_joint_tc(half_cut)
            if tc_fit is None:
                continue
            for window in WINDOWS:
                search.TAU_WINDOW = window
                for method, radius_rows, radius_col in (
                        ("half-radius", half_cut, "r50"),
                        ("Gaussian-radius", gaussian_cut, "r50_gaussian")):
                    fit = search.fit_exponents(
                        half_cut, radius_rows, radius_col, tc_fit["tc"], cutoff, offset)
                    if fit is None:
                        continue
                    output.append({
                        "N": 1000, "offset_percent_from_N500_Acrit": offset,
                        "amplitude": search.A_CRIT_REFERENCE_N500*(1+offset/100),
                        "search_selected_amplitude": int(offset == 0.03125),
                        "radius_method": method, "cutoff_phi_center": cutoff,
                        "cutoff_crossing_time": crossing_time,
                        "tc_method": "threshold_joint_tail",
                        "tc": tc_fit["tc"], "tc_se": tc_fit["tc_se"],
                        "tau_min": window[0], "tau_max": window[1],
                        "alpha": fit["alpha"], "alpha_se": fit["alpha_se"],
                        "beta": fit["beta"], "beta_se": fit["beta_se"],
                        "alpha_plus_beta": fit["alpha_plus_beta"],
                        "alpha_plus_beta_se": fit["alpha_plus_beta_se"],
                        "alpha_r2": fit["alpha_r2"], "beta_r2": fit["beta_r2"],
                        "fit_points": fit["fit_points"],
                    })

    RESULTS.mkdir(parents=True, exist_ok=True)
    table_path = RESULTS / "fit_window_sensitivity_N1000.csv"
    if output:
        with table_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(output[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(output)

    window_labels = [f"{lo:g}–{hi:g}" for lo, hi in WINDOWS]
    colors = {0.01: "#1d4ed8", 0.03125: "#dc2626",
              0.05: "#059669", 0.1: "#7c3aed"}
    markers = {20: "o", 24: "P", 40: "s", 60: "^", 80: "D"}
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.6), sharey=True)
    for ax, method in zip(axes, ("half-radius", "Gaussian-radius")):
        for offset in OFFSETS:
            for cutoff in CUTOFFS:
                selected = [row for row in output
                            if row["radius_method"] == method
                            and row["offset_percent_from_N500_Acrit"] == offset
                            and row["cutoff_phi_center"] == cutoff]
                selected.sort(key=lambda row: row["tau_min"])
                if not selected:
                    continue
                xs = [window_labels.index(
                    f"{row['tau_min']:g}–{row['tau_max']:g}") for row in selected]
                ax.errorbar(xs,
                            [row["alpha_plus_beta"] for row in selected],
                            yerr=[row["alpha_plus_beta_se"] for row in selected],
                            fmt=markers[cutoff], color=colors[offset], ms=5,
                            capsize=1.5, alpha=0.78, linestyle="none")
        ax.axhline(1.0, color="#111827", linestyle="--", linewidth=1)
        ax.set_xticks(range(len(window_labels)), window_labels)
        ax.set_xlabel(r"Remaining-time fit window $\tau=t_c-t$")
        ax.set_title(method)
        ax.grid(alpha=0.2)
    axes[0].set_ylabel(r"Freely fitted $\alpha+\beta$")
    fig.suptitle("Exponent-sum sensitivity to fit window, cutoff, and amplitude")
    handles = [plt.Line2D([0], [0], color=colors[offset], marker="o", linestyle="none",
                          label=f"+{offset:g}%" + (" (search-selected)" if offset == 0.03125 else ""))
               for offset in OFFSETS]
    handles.extend(plt.Line2D([0], [0], color="#555", marker=markers[cutoff], linestyle="none",
                              label=f"$|\\phi_0|={cutoff}$") for cutoff in CUTOFFS)
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8,
               bbox_to_anchor=(0.5, -0.015))
    fig.tight_layout(rect=(0, 0.10, 1, 0.94))
    plot_path = RESULTS / "fit_window_sensitivity_N1000.png"
    fig.savefig(plot_path, dpi=180, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {len(output)} cached-trajectory fits to {table_path}")
    print(f"Wrote {plot_path}")


if __name__ == "__main__":
    main()
