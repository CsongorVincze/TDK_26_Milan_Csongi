#!/usr/bin/env python3
"""Scan collapse-time estimates from progressively later field thresholds."""

import argparse
import csv
import importlib.util
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
HALF_DIR = HERE/"Fraction_of_central"/"half_radius"
GAUSSIAN_DIR = HERE/"gaussian_profile"
HALF_RESULTS = HALF_DIR/"results"
GAUSSIAN_RESULTS = GAUSSIAN_DIR/"results"
HALF_PLOTS = HALF_RESULTS/"alpha_beta_scaling"
GAUSSIAN_PLOTS = GAUSSIAN_RESULTS/"alpha_beta_scaling"

N = 1000
A_CRIT_REFERENCE_N500 = 1.202609700895846
T_MAX, SAMPLE_INTERVAL, STOP_LEVEL = 100.0, 0.01, 80.0
TAU_WINDOW = (0.1, 0.3)
BASE_OFFSETS_PERCENT = (0.01, 0.05, 0.1)
COARSE_CUTOFFS = (20, 30, 40, 50, 60, 70, 80)
MIN_FIT_POINTS = 5


def load_study(name, path):
    """Load a sibling standalone solver without making either study depend on this file."""
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


HALF_STUDY = load_study("half_radius_tc_search_solver", HALF_DIR/"study.py")
GAUSSIAN_STUDY = load_study("gaussian_tc_search_solver", GAUSSIAN_DIR/"study.py")


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key, value in row.items():
            try:
                row[key] = float(value)
            except (TypeError, ValueError):
                pass
    return rows


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def offset_slug(offset_percent):
    return f"{offset_percent:g}pct"


def trajectory_paths(offset_percent):
    slug = offset_slug(offset_percent)
    if math.isclose(offset_percent, 0.01, abs_tol=1e-10):
        return (HALF_RESULTS/"continuation_N1000_phi80.csv",
                GAUSSIAN_RESULTS/"supercritical_plus_0.01pct_N1000.csv")
    return (HALF_RESULTS/f"tc_precision_N1000_plus_{slug}.csv",
            GAUSSIAN_RESULTS/f"tc_precision_N1000_plus_{slug}.csv")


def simulate_case(offset_percent, force=False):
    half_path, gaussian_path = trajectory_paths(offset_percent)
    if half_path.exists() and gaussian_path.exists() and not force:
        return read_csv(half_path), read_csv(gaussian_path), half_path, gaussian_path

    amplitude = A_CRIT_REFERENCE_N500*(1.0+offset_percent/100.0)
    case = f"tc_precision_N{N}_plus_{offset_slug(offset_percent)}"
    print(f"Simulating N={N}, +{offset_percent:g}%: A={amplitude:.12f}, "
          f"stop at |phi(0)|={STOP_LEVEL:g}", flush=True)

    half_rows, h_crossed, h_failed, h_event, h_steps = HALF_STUDY.evolve(
        amplitude, N, T_MAX, SAMPLE_INTERVAL, STOP_LEVEL,
    )
    gaussian_rows, _, g_crossed, g_failed, g_event, g_steps = GAUSSIAN_STUDY.evolve(
        amplitude, N, T_MAX, STOP_LEVEL, SAMPLE_INTERVAL, case,
        capture_profiles=False, fit_radius=8.0,
    )
    if h_failed or g_failed:
        raise RuntimeError(
            f"N={N}, +{offset_percent:g}% encountered a numerical failure "
            f"(half={h_failed}, Gaussian={g_failed})"
        )
    write_csv(half_path, half_rows)
    write_csv(gaussian_path, gaussian_rows)
    h_status = f"reached {h_event:.6f}" if h_crossed else f"not reached by t={T_MAX:g}"
    g_status = f"reached {g_event:.6f}" if g_crossed else f"not reached by t={T_MAX:g}"
    print(f"  saved trajectories; |phi0|={STOP_LEVEL:g}: half {h_status}, "
          f"Gaussian {g_status}; steps {h_steps}, {g_steps}", flush=True)
    return half_rows, gaussian_rows, half_path, gaussian_path


def interpolate_crossing(rows, level):
    """Return rows truncated at the first upward crossing of |phi(0)|=level."""
    center = lambda row: abs(float(row["abs_phi_center"]))
    if not rows or center(rows[-1]) < level:
        return None, None
    if center(rows[0]) >= level:
        return [rows[0]], float(rows[0]["time"])

    for index in range(1, len(rows)):
        before, after = rows[index-1], rows[index]
        left, right = center(before), center(after)
        if left < level <= right:
            fraction = (level-left)/(right-left)
            crossing = {}
            for key, value in before.items():
                try:
                    crossing[key] = float(value)+fraction*(float(after[key])-float(value))
                except (TypeError, ValueError):
                    crossing[key] = value
            crossing["abs_phi_center"] = float(level)
            crossing["time"] = float(before["time"]+fraction*(after["time"]-before["time"]))
            result = rows[:index]+[crossing]
            return result, crossing["time"]
    return None, None


def line_fit(x, y):
    """OLS line and conventional one-standard-error slope/intercept uncertainties."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    n = len(x)
    if n < MIN_FIT_POINTS or np.ptp(x) == 0:
        return None
    xbar, ybar = float(x.mean()), float(y.mean())
    sxx = float(np.sum((x-xbar)**2))
    slope = float(np.sum((x-xbar)*(y-ybar))/sxx)
    intercept = ybar-slope*xbar
    residuals = y-(slope*x+intercept)
    rss = float(np.sum(residuals**2))
    total = float(np.sum((y-ybar)**2))
    dof = n-2
    variance = rss/dof
    slope_variance = variance/sxx
    intercept_variance = variance*(1.0/n+xbar*xbar/sxx)
    slope_intercept_covariance = -xbar*variance/sxx
    return {
        "slope": slope, "intercept": float(intercept),
        "r2": 1.0-rss/total if total > 0 else float("nan"),
        "slope_se": math.sqrt(max(0.0, slope_variance)),
        "intercept_se": math.sqrt(max(0.0, intercept_variance)),
        "slope_intercept_covariance": slope_intercept_covariance,
        "slope_variance": slope_variance,
        "residuals": residuals, "sxx": sxx, "dof": dof, "n": n,
    }


def threshold_joint_tc(truncated_rows):
    """Jointly fit tc and alpha to the large-field tail after each cutoff."""
    tail = [row for row in truncated_rows if float(row["abs_phi_center"]) >= 4.0][-40:]
    if len(tail) < MIN_FIT_POINTS:
        return None
    time = np.asarray([float(row["time"]) for row in tail])
    log_center = np.log10([float(row["abs_phi_center"]) for row in tail])
    left, right = time[-1]+1e-5, time[-1]+1.0
    best = None
    for _ in range(4):
        trial_tc = np.linspace(left, right, 1201)
        log_tau = np.log10(trial_tc[:, None]-time[None, :])
        centered_x = log_tau-log_tau.mean(axis=1, keepdims=True)
        centered_y = log_center-log_center.mean()
        slopes = np.sum(centered_x*centered_y, axis=1)/np.sum(centered_x**2, axis=1)
        offsets = log_center.mean()-slopes*log_tau.mean(axis=1)
        rss = np.sum((log_center[None, :]-(slopes[:, None]*log_tau+offsets[:, None]))**2,
                     axis=1)
        rss[slopes >= 0] = np.inf
        index = int(np.argmin(rss))
        step = (right-left)/(len(trial_tc)-1)
        best = (float(trial_tc[index]), float(-slopes[index]), float(rss[index]))
        left, right = max(time[-1]+1e-6, best[0]-5*step), best[0]+5*step

    tc, alpha_tail, rss = best
    tau = tc-time
    x = np.log10(tau)
    fit = line_fit(x, log_center)
    if fit is None or tc <= float(truncated_rows[-1]["time"]):
        return None
    jacobian = np.column_stack((
        np.ones_like(tau), -x, -alpha_tail/(np.log(10.0)*tau),
    ))
    dof = len(time)-3
    if dof <= 0:
        return None
    try:
        covariance = (rss/dof)*np.linalg.inv(jacobian.T@jacobian)
    except np.linalg.LinAlgError:
        return None
    return {
        "tc": tc, "tc_se": math.sqrt(max(0.0, float(covariance[2, 2]))),
        "tc_fit_r2": fit["r2"], "tc_fit_points": len(tail),
        "alpha_tail": alpha_tail,
    }


def matched_field(radius_rows, central_rows):
    """Map common central-field samples onto a radius trajectory by sample time."""
    return {round(float(row["time"]), 8): float(row["abs_phi_center"])
            for row in central_rows}


def fit_exponents(central_rows, radius_rows, radius_column, tc, cutoff, offset_percent):
    center_by_time = matched_field(radius_rows, central_rows)
    lower, upper = TAU_WINDOW
    xs, fields, radii = [], [], []
    for row in radius_rows:
        time = float(row["time"])
        field = center_by_time.get(round(time, 8))
        radius = float(row[radius_column])
        tau = tc-time
        if (field is not None and lower <= tau <= upper and field >= 4.0
                and np.isfinite(radius) and radius > 0.0):
            xs.append(math.log10(tau))
            fields.append(math.log10(field))
            radii.append(math.log10(radius))
    alpha_fit = line_fit(xs, fields)
    beta_fit = line_fit(xs, radii)
    if alpha_fit is None or beta_fit is None:
        return None

    alpha = -alpha_fit["slope"]
    beta = beta_fit["slope"]
    # The same x samples are used for both slopes; include their residual covariance.
    if (alpha_fit["n"] == beta_fit["n"] and
            math.isclose(alpha_fit["sxx"], beta_fit["sxx"], rel_tol=1e-10)):
        slope_cov = float(np.dot(alpha_fit["residuals"], beta_fit["residuals"]))
        slope_cov /= alpha_fit["dof"]*alpha_fit["sxx"]
        sum_variance = (alpha_fit["slope_variance"]+beta_fit["slope_variance"]
                        -2.0*slope_cov)
    else:
        sum_variance = alpha_fit["slope_variance"]+beta_fit["slope_variance"]

    return {
        "N": N, "offset_percent": offset_percent,
        "cutoff_phi_center": cutoff,
        "tau_min": lower, "tau_max": upper, "tc": tc,
        "alpha": alpha, "alpha_se": alpha_fit["slope_se"],
        "alpha_r2": alpha_fit["r2"],
        "beta": beta, "beta_se": beta_fit["slope_se"],
        "beta_r2": beta_fit["r2"],
        "alpha_plus_beta": alpha+beta,
        "alpha_plus_beta_se": math.sqrt(max(0.0, sum_variance)),
        "fit_points": alpha_fit["n"],
    }


def evaluate_cutoffs(half_rows, gaussian_rows, offset_percent, cutoffs):
    half_results, gaussian_results = [], []
    for cutoff in sorted(set(float(value) for value in cutoffs)):
        half_cut, crossing_time = interpolate_crossing(half_rows, cutoff)
        gaussian_cut, _ = interpolate_crossing(gaussian_rows, cutoff)
        if half_cut is None or gaussian_cut is None:
            continue
        tc = threshold_joint_tc(half_cut)
        if tc is None:
            continue
        shared = {
            "tc_estimator": "threshold_joint_tail",
            "cutoff_crossing_time": crossing_time,
            "tc": tc["tc"], "tc_se": tc["tc_se"],
            "tc_fit_r2": tc["tc_fit_r2"], "tc_fit_points": tc["tc_fit_points"],
            "alpha_tail": tc["alpha_tail"],
        }
        half_fit = fit_exponents(
            half_cut, half_cut, "r50", tc["tc"], cutoff, offset_percent)
        gaussian_fit = fit_exponents(
            half_cut, gaussian_cut, "r50_gaussian", tc["tc"], cutoff, offset_percent)
        if half_fit is not None:
            half_results.append({**shared, **half_fit})
        if gaussian_fit is not None:
            gaussian_results.append({**shared, **gaussian_fit})
    return half_results, gaussian_results


def best_search_point(results_by_method):
    choices = [row for rows in results_by_method for row in rows]
    return min(choices, key=lambda row: abs(row["alpha_plus_beta"]-1.0)) if choices else None


def refine_cutoffs(half_rows, gaussian_rows, offset_percent):
    coarse_half, coarse_gaussian = evaluate_cutoffs(
        half_rows, gaussian_rows, offset_percent, COARSE_CUTOFFS)
    coarse_rows = coarse_half+coarse_gaussian
    if not coarse_rows:
        return list(COARSE_CUTOFFS), coarse_half, coarse_gaussian

    # Resolve every coarse threshold region already close to the target, not just
    # whichever one happens to be the single closest point.
    centers = {row["cutoff_phi_center"] for row in coarse_rows
               if abs(row["alpha_plus_beta"]-1.0) <= 0.20}
    if not centers:
        centers = {best_search_point((coarse_half, coarse_gaussian))["cutoff_phi_center"]}
    fine = set()
    for center in centers:
        fine.update(round(value, 4) for value in np.arange(
            max(20.0, center-5.0), min(STOP_LEVEL, center+5.0)+0.5, 1.0))
    cutoffs = sorted(set(COARSE_CUTOFFS).union(fine))
    return cutoffs, *evaluate_cutoffs(half_rows, gaussian_rows, offset_percent, cutoffs)


def refine_offsets(base_offsets, coarse_results):
    choices = []
    for offset, methods in coarse_results.items():
        row = best_search_point(methods)
        if row is not None:
            choices.append((abs(row["alpha_plus_beta"]-1.0), offset))
    if not choices:
        return list(base_offsets), None
    best_offset = min(choices)[1]
    ordered = sorted(base_offsets)
    index = ordered.index(best_offset)
    intervals = []
    if index > 0:
        intervals.append((ordered[index-1], ordered[index]))
    if index+1 < len(ordered):
        intervals.append((ordered[index], ordered[index+1]))
    refined = []
    for lower, upper in intervals:
        refined.extend(lower+(upper-lower)*step/4.0 for step in (1, 2, 3))
    return sorted(set(ordered+refined)), best_offset


def refine_crossing_offsets(offsets, coarse_results, subdivisions=8):
    """Densify amplitude intervals that bracket alpha+beta=1 at matching settings."""
    additions = set()
    ordered = sorted(offsets)
    for lower, upper in zip(ordered, ordered[1:]):
        for method_index in (0, 1):
            left_rows = coarse_results[lower][method_index]
            right_rows = coarse_results[upper][method_index]
            left = {(row["tc_estimator"], row["cutoff_phi_center"]): row["alpha_plus_beta"]
                    for row in left_rows}
            right = {(row["tc_estimator"], row["cutoff_phi_center"]): row["alpha_plus_beta"]
                     for row in right_rows}
            if any((left[key]-1.0)*(right[key]-1.0) < 0
                   for key in left.keys() & right.keys()):
                additions.update(lower+(upper-lower)*step/subdivisions
                                 for step in range(1, subdivisions))
    return sorted(additions-set(ordered))


def plot_case(rows, offset_percent, radius_label, output):
    if not rows:
        return
    fig, ax = plt.subplots(figsize=(8.8, 6.0))
    series = (
        ("alpha", "alpha_se", r"$\alpha$", "#7c3aed", "o"),
        ("beta", "beta_se", r"$\beta$", "#2563eb", "s"),
        ("alpha_plus_beta", "alpha_plus_beta_se", r"$\alpha+\beta$", "#15803d", "^"),
    )
    points = sorted(rows, key=lambda row: row["tc"])
    for value_key, error_key, label, color, marker in series:
        ax.errorbar(
            [row["tc"] for row in points],
            [row[value_key] for row in points],
            yerr=[row[error_key] for row in points],
            xerr=[row["tc_se"] for row in points],
            fmt=f"{marker}-", color=color, capsize=2, lw=1.4, ms=4,
            label=label,
        )
    ax.axhline(1.0, color="#555555", linestyle="--", linewidth=1.0,
               label=r"target $\alpha+\beta=1$")
    ax.set_xlabel(r"Joint high-field-tail estimate of collapse time $t_c$")
    ax.set_ylabel("Fitted exponent")
    ax.set_title(r"Joint central-field tail fit for $(t_c,\alpha)$")
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=9)
    fig.suptitle(
        rf"$N={N}$, $A=A_{{crit}}^{{(N=500)}}(1+{offset_percent:g}\%)$, "
        + rf"$\tau\in[{TAU_WINDOW[0]},{TAU_WINDOW[1]}]$"
        + rf"; $\beta$ from {radius_label}", fontsize=13)
    fig.text(0.5, 0.015,
             "Vertical bars: OLS slope SE; horizontal bars: joint-fit tc SE. "
             "Threshold levels are listed in the CSV.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(output, dpi=180)
    plt.close(fig)


def write_outputs(all_rows, manifests):
    HALF_PLOTS.mkdir(parents=True, exist_ok=True)
    GAUSSIAN_PLOTS.mkdir(parents=True, exist_ok=True)
    half_rows = [row for offset_rows in all_rows.values() for row in offset_rows[0]]
    gaussian_rows = [row for offset_rows in all_rows.values() for row in offset_rows[1]]
    write_csv(HALF_PLOTS/"tc_precision_search.csv", half_rows)
    write_csv(GAUSSIAN_PLOTS/"tc_precision_search.csv", gaussian_rows)
    write_csv(HALF_PLOTS/"tc_precision_run_manifest.csv", manifests)
    write_csv(GAUSSIAN_PLOTS/"tc_precision_run_manifest.csv", manifests)
    for offset, (half_fit_rows, gaussian_fit_rows) in all_rows.items():
        slug = offset_slug(offset)
        plot_case(half_fit_rows, offset, r"$R_{50}$",
                  HALF_PLOTS/f"tc_sensitivity_N1000_plus_{slug}.png")
        plot_case(gaussian_fit_rows, offset, r"$R_{50,G}$",
                  GAUSSIAN_PLOTS/f"tc_sensitivity_N1000_plus_{slug}.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offset-percent", type=float, nargs="+",
                        default=BASE_OFFSETS_PERCENT,
                        help="coarse supercritical offsets in percent (default: 0.01 0.05 0.1)")
    parser.add_argument("--no-dense-amplitudes", action="store_true",
                        help="skip the data-driven amplitude refinement")
    parser.add_argument("--force-runs", action="store_true",
                        help="rerun N=1000 trajectories for non-reference offsets")
    args = parser.parse_args()
    base_offsets = sorted(set(args.offset_percent))
    if not base_offsets or any(value <= 0 for value in base_offsets):
        parser.error("All supercritical offsets must be positive percentages.")

    trajectories, manifests = {}, []
    coarse_results = {}

    def ensure_coarse(offsets, force=False):
        for offset in offsets:
            if offset in coarse_results:
                continue
            half_rows, gaussian_rows, half_path, gaussian_path = simulate_case(
                offset, force=force)
            trajectories[offset] = (half_rows, gaussian_rows)
            half_coarse, gaussian_coarse = evaluate_cutoffs(
                half_rows, gaussian_rows, offset, COARSE_CUTOFFS)
            coarse_results[offset] = (half_coarse, gaussian_coarse)
            manifests.append({
                "N": N, "offset_percent": offset,
                "amplitude": A_CRIT_REFERENCE_N500*(1+offset/100.0),
                "Acrit_reference_N500": A_CRIT_REFERENCE_N500,
                "stop_level": STOP_LEVEL, "sample_interval": SAMPLE_INTERVAL,
                "half_reached_stop": int(half_rows[-1]["abs_phi_center"] >= STOP_LEVEL-1e-8),
                "gaussian_reached_stop": int(gaussian_rows[-1]["abs_phi_center"] >= STOP_LEVEL-1e-8),
                "half_last_time": half_rows[-1]["time"],
                "gaussian_last_time": gaussian_rows[-1]["time"],
                "half_trajectory_csv": half_path.name,
                "gaussian_trajectory_csv": gaussian_path.name,
            })
            print(f"Scanned {len(half_coarse)} half-radius and "
                  f"{len(gaussian_coarse)} Gaussian coarse threshold estimates "
                  f"for +{offset:g}%", flush=True)

    ensure_coarse(base_offsets, force=args.force_runs)

    best_coarse = None
    if not args.no_dense_amplitudes:
        offsets, best_coarse = refine_offsets(base_offsets, coarse_results)
    else:
        offsets = base_offsets
    if best_coarse is not None:
        print(f"Closest coarse sum to 1 occurred at +{best_coarse:g}%; "
              "densifying adjacent amplitude interval(s).", flush=True)

    ensure_coarse(offsets, force=args.force_runs)
    if not args.no_dense_amplitudes:
        crossing_offsets = refine_crossing_offsets(offsets, coarse_results)
        if crossing_offsets:
            intervals = [(lower, upper) for lower, upper in
                         zip(sorted(offsets), sorted(offsets)[1:])
                         if any(lower < value < upper for value in crossing_offsets)]
            print(f"Found alpha+beta=1 bracket(s) in {intervals}; "
                  f"adding {len(crossing_offsets)} denser amplitude points.", flush=True)
            offsets = sorted(set(offsets+crossing_offsets))
            ensure_coarse(crossing_offsets, force=args.force_runs)

    final_rows = {}
    for offset in offsets:
        half_rows, gaussian_rows = trajectories[offset]
        cutoffs, half_fits, gaussian_fits = refine_cutoffs(
            half_rows, gaussian_rows, offset)
        final_rows[offset] = (half_fits, gaussian_fits)
        best = best_search_point((half_fits, gaussian_fits))
        if best:
            print(f"  +{offset:g}%: {len(cutoffs)} threshold levels; closest "
                  f"sum={best['alpha_plus_beta']:.4f} at |phi0|={best['cutoff_phi_center']:g}",
                  flush=True)
    write_outputs(final_rows, manifests)
    half_plots = sum(bool(rows[0]) for rows in final_rows.values())
    gaussian_plots = sum(bool(rows[1]) for rows in final_rows.values())
    print(f"Wrote {half_plots} half-radius and {gaussian_plots} Gaussian-radius "
          "supercritical-value plots; summary tables are in both "
          "alpha_beta_scaling folders.")


if __name__ == "__main__":
    main()
