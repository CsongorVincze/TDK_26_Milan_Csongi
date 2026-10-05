#!/usr/bin/env python3
"""Fit central-field alpha and two core-radius betas using a shared t_c."""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


HERE = Path(__file__).resolve().parent
FRACTION = HERE/"Fraction_of_central"/"half_radius"/"results"
GAUSSIAN = HERE/"gaussian_profile"/"results"
HALF_PLOT_DIR = FRACTION/"alpha_beta_scaling"
GAUSSIAN_PLOT_DIR = GAUSSIAN/"alpha_beta_scaling"
# Keep the near-collapse bracket only; the remaining-time upper bound is < 0.4.
# Focused near-collapse bracket; discard intervals whose upper tau bound > 0.4.
WINDOWS = ((0.1, 0.3),)
# Keep only matched runs with N >= 1000; currently only this case is available.
RUNS = (
    ("plus_0.01pct_N1000", FRACTION/"continuation_N1000_phi80.csv",
     GAUSSIAN/"supercritical_plus_0.01pct_N1000.csv", 1000, 0.01),
)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = []
        for row in csv.DictReader(handle):
            for key, value in row.items():
                try:
                    row[key] = float(value)
                except (TypeError, ValueError):
                    pass
            rows.append(row)
        return rows


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fit_line(x, y):
    if len(x) < 5 or np.ptp(x) == 0:
        return np.nan, np.nan, np.nan
    slope, intercept = np.polyfit(x, y, 1)
    residual = y-(slope*x+intercept)
    total = np.sum((y-y.mean())**2)
    r2 = 1.0-np.sum(residual**2)/total if total > 0 else np.nan
    return float(slope), float(intercept), float(r2)


def shared_tc(central_rows):
    """Estimate tc and alpha jointly from the saved high-field tail."""
    time = np.asarray([row["time"] for row in central_rows])
    center = np.asarray([row["abs_phi_center"] for row in central_rows])
    ids = np.flatnonzero(center >= 4.0)[-40:]
    if len(ids) < 8:
        raise ValueError("Need at least eight tail samples with |phi(0)| >= 4")

    fit_time, log_center = time[ids], np.log(center[ids])
    left, right = time[-1]+1e-5, time[-1]+1.0
    for _ in range(4):
        trial_tc = np.linspace(left, right, 1201)
        log_tau = np.log(trial_tc[:, None]-fit_time[None, :])
        centered_x = log_tau-log_tau.mean(axis=1, keepdims=True)
        centered_y = log_center-log_center.mean()
        slopes = np.sum(centered_x*centered_y, axis=1)/np.sum(centered_x**2, axis=1)
        offsets = log_center.mean()-slopes*log_tau.mean(axis=1)
        rss = np.sum((log_center[None, :]-(slopes[:, None]*log_tau+offsets[:, None]))**2,
                     axis=1)
        rss[slopes >= 0] = np.inf
        index = int(np.argmin(rss))
        step = (right-left)/(len(trial_tc)-1)
        tc, alpha, best_rss = (float(trial_tc[index]), float(-slopes[index]),
                               float(rss[index]))
        left, right = max(time[-1]+1e-6, tc-5*step), tc+5*step
    total = np.sum((log_center-log_center.mean())**2)
    joint_r2 = 1.0-best_rss/total if total > 0 else np.nan
    joint = (tc, alpha, float(joint_r2), len(ids))
    return {"joint_tc_alpha": joint}


def fit_observable(rows, radius_column, case, n, offset, method, tc_fit,
                   tc_source, central_rows, windows=WINDOWS):
    time = np.asarray([row["time"] for row in rows])
    center = np.asarray([row["abs_phi_center"] for row in rows])
    radius = np.asarray([row[radius_column] for row in rows])
    central_time = np.asarray([row["time"] for row in central_rows])
    central_field = np.asarray([row["abs_phi_center"] for row in central_rows])
    tc, alpha_tail, tc_r2, tc_n = tc_fit
    tau = tc-time
    central_tau = tc-central_time
    fits = []
    for lower, upper in windows:
        use = ((tau >= lower) & (tau <= upper) & (center >= 4.0)
               & np.isfinite(radius) & (radius > 0.0))
        alpha_use = ((central_tau >= lower) & (central_tau <= upper)
                     & (central_field >= 4.0) & np.isfinite(central_field))
        beta, log_c_r, r2 = fit_line(np.log10(tau[use]), np.log10(radius[use]))
        central_slope, log_c_phi, alpha_r2 = fit_line(
            np.log10(central_tau[alpha_use]), np.log10(central_field[alpha_use]))
        alpha = -central_slope
        fits.append({
            "case": case, "N": n, "offset_percent": offset,
            "observable": radius_column, "tc_method": method,
            "tc_source": tc_source, "tc": tc,
            "alpha_central": alpha,
            "log10_C_phi": log_c_phi,
            "C_phi": float(10.0**log_c_phi) if np.isfinite(log_c_phi) else np.nan,
            "alpha_fit_r2": alpha_r2,
            "alpha_fit_points": int(alpha_use.sum()),
            "alpha_tail_joint": alpha_tail,
            "tc_fit_r2": tc_r2, "tc_fit_points": tc_n,
            "tau_min": lower, "tau_max": upper, "fit_points": int(use.sum()),
            "beta": beta, "log10_C_R": log_c_r,
            "C_R": float(10.0**log_c_r) if np.isfinite(log_c_r) else np.nan,
            "beta_fit_r2": r2,
        })
    return fits


def pair_fits(rows):
    by_key = {(r["case"], r["tc_method"], r["tau_min"], r["tau_max"],
               r["N"], r["observable"]): r
              for r in rows}
    comparison = []
    cases = sorted({r["case"] for r in rows})
    for case in cases:
        for method, lower, upper, n in sorted({
                (r["tc_method"], r["tau_min"], r["tau_max"], r["N"])
                for r in rows if r["case"] == case}):
            half = by_key.get((case, method, lower, upper, n, "r50"))
            gauss = by_key.get((case, method, lower, upper, n, "r50_gaussian"))
            if half is None or gauss is None:
                continue
            comparison.append({
                "case": case, "N": n, "tc_method": method,
                "tc": half["tc"],
                "tau_min": lower, "tau_max": upper,
                "alpha_central": half["alpha_central"],
                "alpha_fit_r2": half["alpha_fit_r2"],
                "alpha_fit_points": half["alpha_fit_points"],
                "alpha_tail_joint": half["alpha_tail_joint"],
                "beta_half_radius": half["beta"],
                "C_R_half_radius": half["C_R"],
                "log10_C_R_half_radius": half["log10_C_R"],
                "alpha_plus_beta_half_radius": half["alpha_central"]+half["beta"],
                "beta_gaussian": gauss["beta"],
                "C_R_gaussian": gauss["C_R"],
                "log10_C_R_gaussian": gauss["log10_C_R"],
                "alpha_plus_beta_gaussian": half["alpha_central"]+gauss["beta"],
                "C_R_gaussian_over_half": (
                    gauss["C_R"]/half["C_R"] if half["C_R"] else np.nan),
                "r2_half_radius": half["beta_fit_r2"],
                "r2_gaussian": gauss["beta_fit_r2"],
                "n_half_radius": half["fit_points"],
                "n_gaussian": gauss["fit_points"],
            })
    return comparison


def plot_alpha_beta_window(rows, central_rows, radius_column, radius_label,
                           n, offset, tc, method, window, output):
    """Plot alpha and beta fits together for one remaining-time window."""
    center_time = np.asarray([row["time"] for row in central_rows])
    center = np.asarray([row["abs_phi_center"] for row in central_rows])
    time = np.asarray([row["time"] for row in rows])
    field = np.asarray([row["abs_phi_center"] for row in rows])
    radius = np.asarray([row[radius_column] for row in rows])
    lower, upper = window
    fig, ax_alpha = plt.subplots(figsize=(8.5, 5.5))
    ax_beta = ax_alpha.twinx()
    tc_label = r"$t_c$ from joint high-field central-field tail fit"
    fig.suptitle(
        rf"$N={n}$, $A=A_{{crit}}(1+{offset:g}\%)$, $\tau\in[{lower:g},{upper:g}]$"
        + "\n" + tc_label,
        fontsize=13)

    central_tau = tc-center_time
    alpha_use = ((central_tau >= lower) & (central_tau <= upper)
                 & (center >= 4.0) & (center > 0.0))
    x_alpha, y_alpha = (np.log10(central_tau[alpha_use]),
                        np.log10(center[alpha_use]))
    alpha_slope, alpha_intercept, alpha_r2 = fit_line(x_alpha, y_alpha)
    alpha = -alpha_slope
    ax_alpha.plot(x_alpha, y_alpha, ".", ms=5, color="#7c3aed")

    tau = tc-time
    beta_use = ((tau >= lower) & (tau <= upper) & (field >= 4.0)
                & np.isfinite(radius) & (radius > 0.0))
    x_beta, y_beta = np.log10(tau[beta_use]), np.log10(radius[beta_use])
    beta, beta_intercept, beta_r2 = fit_line(x_beta, y_beta)
    ax_beta.plot(x_beta, y_beta, ".", ms=5, color="#2563eb")

    x_fit = np.linspace(np.log10(lower), np.log10(upper), 80)
    alpha_label = (rf"Central field: $\alpha={alpha:.3f}$, $R^2={alpha_r2:.4f}$, "
                   rf"$n={int(alpha_use.sum())}$"
                   if np.isfinite(alpha) else
                   f"Central field: no fit (n={int(alpha_use.sum())})")
    beta_label = (rf"${radius_label}$: $\beta={beta:.3f}$, $R^2={beta_r2:.4f}$, "
                  rf"$n={int(beta_use.sum())}$"
                  if np.isfinite(beta) else
                  rf"${radius_label}$: no fit (n={int(beta_use.sum())})")
    if np.isfinite(alpha_slope):
        ax_alpha.plot(x_fit, alpha_slope*x_fit+alpha_intercept,
                      "--", lw=1.7, color="#7c3aed")
    if np.isfinite(beta):
        ax_beta.plot(x_fit, beta*x_fit+beta_intercept,
                     "--", lw=1.7, color="#2563eb")

    ax_alpha.set_xlabel(r"Shared time coordinate: $\log_{10}(t_c-t)$")
    ax_alpha.set_ylabel(r"Central field: $\log_{10}|\phi(0,t)|$", color="#7c3aed")
    ax_beta.set_ylabel(rf"Core radius: $\log_{{10}}{radius_label}$", color="#2563eb")
    ax_alpha.tick_params(axis="y", colors="#7c3aed")
    ax_beta.tick_params(axis="y", colors="#2563eb")
    ax_alpha.set_xlim(np.log10(lower), np.log10(upper))
    ax_alpha.grid(alpha=0.25)

    handles = [
        Line2D([0], [0], color="#7c3aed", marker=".", linestyle="--",
               label=alpha_label),
        Line2D([0], [0], color="#2563eb", marker=".", linestyle="--",
               label=beta_label),
    ]
    ax_alpha.legend(handles=handles, loc="best", fontsize=9, framealpha=0.9)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(output, dpi=170)
    plt.close(fig)


def plot_all_alpha_beta():
    """Write matching central-field/radius fit figures to each study's results."""
    HALF_PLOT_DIR.mkdir(parents=True, exist_ok=True)
    GAUSSIAN_PLOT_DIR.mkdir(parents=True, exist_ok=True)
    # Remove only obsolete figures produced by the retired fixed-alpha method.
    for plot_dir in (HALF_PLOT_DIR, GAUSSIAN_PLOT_DIR):
        for stale_plot in plot_dir.glob("*tc_alpha1*.png"):
            stale_plot.unlink()
    count = 0
    for case, half_path, gaussian_path, n, offset in RUNS:
        central_rows, gaussian_rows = read_csv(half_path), read_csv(gaussian_path)
        for method, tc_fit in shared_tc(central_rows).items():
            if tc_fit is None:
                continue
            tc = tc_fit[0]
            slug = case.replace("plus_", "").replace("pct_", "pct_")
            suffix = "tc_joint_tail"
            for lower, upper in WINDOWS:
                filename = (f"alpha_beta_{slug}_{suffix}_tau_"
                            f"{lower:g}_{upper:g}.png")
                plot_alpha_beta_window(
                    central_rows, central_rows, "r50", r"R_{50}",
                    n, offset, tc, method, (lower, upper),
                    HALF_PLOT_DIR/filename)
                plot_alpha_beta_window(
                    gaussian_rows, central_rows, "r50_gaussian", r"R_{50,G}",
                    n, offset, tc, method, (lower, upper),
                    GAUSSIAN_PLOT_DIR/filename)
                count += 2
    return count


def main():
    half_fits, gaussian_fits = [], []
    for case, half_path, gaussian_path, n, offset in RUNS:
        half_rows, gaussian_rows = read_csv(half_path), read_csv(gaussian_path)
        estimators = shared_tc(half_rows)
        for method, tc_fit in estimators.items():
            if tc_fit is None:
                continue
            half_fits.extend(fit_observable(
                half_rows, "r50", case, n, offset, method, tc_fit,
                half_path.name, half_rows,
            ))
            gaussian_fits.extend(fit_observable(
                gaussian_rows, "r50_gaussian", case, n, offset, method,
                tc_fit, half_path.name, half_rows,
            ))
    comparison = pair_fits(half_fits+gaussian_fits)
    write_csv(FRACTION/"shared_tc_alpha_fits.csv", half_fits)
    write_csv(GAUSSIAN/"shared_tc_alpha_fits.csv", gaussian_fits)
    write_csv(HERE/"shared_tc_alpha_comparison.csv", comparison)
    count = plot_all_alpha_beta()
    print(f"Wrote {count} individual alpha-beta figures ({count//2} per results folder)")


if __name__ == "__main__":
    main()
