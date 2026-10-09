#!/usr/bin/env python3
"""Map exponent-sum drift across amplitude, field range, tc fit, and timestep.

Runs the established Module 4 solver and reuses the audit's R50/PDE diagnostics.
All outputs live in this experiment's own results folder.
"""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import run_audit as audit


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "regime_scan" / "results"
TC_BANDS = ((4, 12), (8, 20), (12, 30), (20, 40), (20, 80),
            (40, 80), (40, 160), (80, 320), (80, 1280),
            (160, 1280), (320, 1280))
FIELD_BANDS = ((4, 12), (12, 30), (30, 80), (80, 200),
               (200, 500), (500, 1280))
TAU_BANDS = ((0.1, 0.3), (0.05, 0.1), (0.01, 0.03), (0.003, 0.01))
PARAMETRIC_FIELD_BANDS = ((6, 12), (8.6, 17), (12, 30), (20, 40),
                          (30, 80), (80, 200), (200, 500), (500, 1280))


def line_fit(x, y):
    x, y = np.asarray(x), np.asarray(y)
    if len(x) < 6 or np.ptp(x) == 0:
        return None
    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (slope*x + intercept)
    total = np.sum((y-y.mean())**2)
    return (float(slope), float(1-np.sum(residual**2)/total)
            if total else float("nan"))


def fit_local_exponents(events, case, tc, tc_band):
    rows = []
    for low, high in FIELD_BANDS:
        selected = [row for row in events if low <= row["level"] <= high
                    and tc > row["time"] and row["r50"] > 0]
        if len(selected) < 6:
            continue
        tau = tc - np.asarray([row["time"] for row in selected])
        log_tau = np.log(tau)
        alpha_fit = line_fit(log_tau, -np.log(
            [row["abs_phi_center"] for row in selected]))
        beta_fit = line_fit(log_tau, np.log([row["r50"] for row in selected]))
        gamma_fit = line_fit(log_tau, -np.log([row["gradient_r50"] for row in selected]))
        if not alpha_fit or not beta_fit or not gamma_fit:
            continue
        alpha, alpha_r2 = alpha_fit
        beta, beta_r2 = beta_fit
        gamma, gamma_r2 = gamma_fit
        rows.append({
            "case": case, "tc_fit_level_min": tc_band[0],
            "tc_fit_level_max": tc_band[1], "tc": tc,
            "field_min": low, "field_max": high, "points": len(selected),
            "tau_min": float(np.min(tau)), "tau_max": float(np.max(tau)),
            "field_geomean": float(np.sqrt(low*high)),
            "alpha": alpha, "alpha_r2": alpha_r2,
            "beta_r50": beta, "beta_r2": beta_r2,
            "alpha_plus_beta": alpha+beta,
            "gamma_gradient": gamma, "gamma_r2": gamma_r2,
            "gamma_minus_alpha_plus_beta": gamma-alpha-beta,
            "median_cells_per_r50": float(np.median(
                [row["cells_per_r50"] for row in selected])),
        })
    return rows


def fit_tau_windows(events, case, tc, tc_band):
    rows = []
    for lower, upper in TAU_BANDS:
        selected = [row for row in events if row["abs_phi_center"] >= 4
                    and lower <= tc-row["time"] <= upper and row["r50"] > 0]
        if len(selected) < 6:
            continue
        tau = tc - np.asarray([row["time"] for row in selected])
        log_tau = np.log(tau)
        alpha_fit = line_fit(log_tau, -np.log(
            [row["abs_phi_center"] for row in selected]))
        beta_fit = line_fit(log_tau, np.log([row["r50"] for row in selected]))
        gamma_fit = line_fit(log_tau, -np.log([row["gradient_r50"] for row in selected]))
        if not alpha_fit or not beta_fit or not gamma_fit:
            continue
        alpha, alpha_r2 = alpha_fit
        beta, beta_r2 = beta_fit
        gamma, gamma_r2 = gamma_fit
        rows.append({
            "case": case, "tc_fit_level_min": tc_band[0],
            "tc_fit_level_max": tc_band[1], "tc": tc,
            "tau_window_min": lower, "tau_window_max": upper,
            "tau_observed_min": float(np.min(tau)),
            "tau_observed_max": float(np.max(tau)), "points": len(selected),
            "phi_min": float(min(row["abs_phi_center"] for row in selected)),
            "phi_max": float(max(row["abs_phi_center"] for row in selected)),
            "alpha": alpha, "alpha_r2": alpha_r2,
            "beta_r50": beta, "beta_r2": beta_r2,
            "alpha_plus_beta": alpha+beta,
            "gamma_gradient": gamma, "gamma_r2": gamma_r2,
            "gamma_minus_alpha_plus_beta": gamma-alpha-beta,
        })
    return rows


def fit_tc_free(events, case):
    """Estimate alpha+beta without estimating tc.

    If u=|phi(0)|~tau^-alpha, then |du/dt|~u^(1+1/alpha).
    If R50~tau^beta, then R50~u^(-beta/alpha). These two slopes
    determine alpha and beta/alpha without a fitted collapse time.
    """
    rows = []
    for low, high in PARAMETRIC_FIELD_BANDS:
        selected = [row for row in events if low <= row["abs_phi_center"] <= high
                    and row.get("phi_center") is not None
                    and row["r50"] > 0 and row["phi_center"]*row["pi_center"] > 0]
        if len(selected) < 6:
            continue
        log_u = np.log([row["abs_phi_center"] for row in selected])
        log_speed = np.log([abs(row["pi_center"]) for row in selected])
        log_accel_ratio = np.log([
            max(row["rhs_over_cubic_center"], np.finfo(float).tiny)
            for row in selected
        ])
        speed_fit = line_fit(log_u, log_speed)
        accel_fit = line_fit(log_u, log_accel_ratio)
        if not speed_fit or not accel_fit:
            continue
        speed_slope, speed_r2 = speed_fit
        accel_slope, accel_r2 = accel_fit
        radius_fits = {}
        for radius_key in ("r25_fraction", "r50", "r75_fraction"):
            radius_rows = [row for row in selected
                           if row.get(radius_key) is not None
                           and row[radius_key] > 0]
            radius_fit = line_fit(
                np.log([row["abs_phi_center"] for row in radius_rows]),
                np.log([row[radius_key] for row in radius_rows]),
            ) if len(radius_rows) >= 6 else None
            if radius_fit:
                q_radius = -radius_fit[0]
                radius_fits[radius_key] = {
                    "q": q_radius, "r2": radius_fit[1],
                    "sum": float("nan"),
                }
            else:
                radius_fits[radius_key] = {
                    "q": float("nan"), "r2": float("nan"),
                    "sum": float("nan"),
                }
        q = radius_fits["r50"]["q"]
        alpha = 1/(speed_slope-1) if speed_slope > 1 else float("nan")
        beta = alpha*q
        for radius_result in radius_fits.values():
            radius_result["sum"] = alpha*(1+radius_result["q"])
        if np.isfinite(alpha):
            tc_local = np.asarray([
                row["time"] + alpha*row["abs_phi_center"]/abs(row["pi_center"])
                for row in selected
            ])
            tc_local_mean = float(np.mean(tc_local))
            tc_local_std = float(np.std(tc_local))
            tc_local_span = float(np.ptp(tc_local))
        else:
            tc_local_mean = tc_local_std = tc_local_span = float("nan")
        rows.append({
            "case": case["name"], "N": case["n"],
            "offset_percent": case["offset"], "field_step": case["field_step"],
            "phi_min": float(min(row["abs_phi_center"] for row in selected)),
            "phi_max": float(max(row["abs_phi_center"] for row in selected)),
            "points": len(selected), "speed_slope_s": speed_slope,
            "speed_r2": speed_r2, "q_beta_over_alpha": q,
            "radius_r2": radius_fits["r50"]["r2"],
            "alpha_from_velocity": alpha,
            "beta_from_radius": beta, "alpha_plus_beta_tc_free": alpha+beta,
            "target_residual_s_minus_2_plus_q": speed_slope-(2+q),
            "tc_local_mean_from_powerlaw": tc_local_mean,
            "tc_local_std_from_powerlaw": tc_local_std,
            "tc_local_span_from_powerlaw": tc_local_span,
            "accel_over_cubic_slope": accel_slope,
            "accel_over_cubic_r2": accel_r2,
            "slope_predicted_from_speed_power": 2*speed_slope-4,
            "q_r25_fraction": radius_fits["r25_fraction"]["q"],
            "r2_r25_fraction": radius_fits["r25_fraction"]["r2"],
            "sum_r25_fraction": radius_fits["r25_fraction"]["sum"],
            "q_r75_fraction": radius_fits["r75_fraction"]["q"],
            "r2_r75_fraction": radius_fits["r75_fraction"]["r2"],
            "sum_r75_fraction": radius_fits["r75_fraction"]["sum"],
        })
    return rows


def plot_regimes(fits, cases, output):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharex=True)
    variables = (("alpha", r"local central exponent $\alpha$"),
                 ("beta_r50", r"local width exponent $\beta_{50}$"),
                 ("alpha_plus_beta", r"local sum $\alpha+\beta_{50}$"))
    colors = plt.cm.viridis(np.linspace(0.05, 0.9, len(cases)))
    for case, color in zip(cases, colors):
        selected = sorted((row for row in fits if row["case"] == case["name"]
                           and row["tc_fit_level_min"] == 160),
                          key=lambda row: row["field_geomean"])
        label = f"+{case['offset']:g}%" + (
            " (dt cap .05)" if case["field_step"] == 0.05 else "")
        for ax, (key, _) in zip(axes, variables):
            ax.plot([row["field_geomean"] for row in selected],
                    [row[key] for row in selected], "o-", color=color,
                    label=label)
    axes[2].axhline(1.0, color="black", ls="--", lw=1,
                    label=r"notes target $1$")
    axes[2].axhline(1.5, color="black", ls=":", lw=1,
                    label=r"local-ODE profile $3/2$")
    for ax, (_, ylabel) in zip(axes, variables):
        ax.set_xscale("log")
        ax.set_xlabel(r"geometric-mean central field in fit band")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.25)
    axes[0].legend(fontsize=7, loc="best")
    axes[2].legend(fontsize=7, loc="best")
    fig.suptitle(r"Local exponents drift across collapse stages "
                 r"(shared $t_c$ fit from $|\phi_0|=160$--$1280$)")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_tau_windows(fits, cases, output):
    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    colors = plt.cm.viridis(np.linspace(0.05, 0.9, len(cases)))
    for case, color in zip(cases, colors):
        selected = sorted((row for row in fits if row["case"] == case["name"]
                           and row["tc_fit_level_min"] == 160),
                          key=lambda row: row["tau_window_min"], reverse=True)
        if not selected:
            continue
        x = np.asarray([np.sqrt(row["tau_window_min"]*row["tau_window_max"])
                        for row in selected])
        y = np.asarray([row["alpha_plus_beta"] for row in selected])
        label = f"+{case['offset']:g}%" + (
            " (dt cap .05)" if case["field_step"] == 0.05 else "")
        ax.plot(x, y, "o-", color=color, label=label)
    ax.axhline(1.0, color="black", ls="--", lw=1, label=r"notes target $1$")
    ax.axhline(1.5, color="black", ls=":", lw=1,
               label=r"high-field trend $3/2$")
    ax.set_xscale("log")
    ax.set_xlabel(r"geometric center of fit window in $\tau=t_c-t$")
    ax.set_ylabel(r"joint-window fit $\alpha+\beta_{50}$")
    ax.set_title(r"The fitted sum changes with collapse stage "
                 r"(same high-field $t_c$ for both exponents)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_tc_drift(tc_rows, cases, output):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    colors = plt.cm.viridis(np.linspace(0.05, 0.9, len(cases)))
    for case, color in zip(cases, colors):
        rows = sorted((row for row in tc_rows if row["case"] == case["name"]),
                      key=lambda row: np.sqrt(row["fit_level_min"]*
                                              row["fit_level_max"]))
        reference = next((row["tc"] for row in rows
                          if row["fit_level_min"] == 320
                          and row["fit_level_max"] == 1280), None)
        if not rows or reference is None:
            continue
        x = [np.sqrt(row["fit_level_min"]*row["fit_level_max"])
             for row in rows]
        label = f"+{case['offset']:g}%" + (
            " (dt cap .05)" if case["field_step"] == 0.05 else "")
        axes[0].semilogx(x, [row["tc"]-reference for row in rows], "o-",
                         color=color, label=label)
        axes[1].semilogx(x, [row["alpha_tail"] for row in rows], "o-",
                         color=color, label=label)
    axes[0].axhline(0, color="black", ls="--", lw=1)
    axes[0].set(xlabel=r"geometric-mean field in $t_c,\alpha$ fit band",
                ylabel=r"$t_c(\mathrm{band})-t_c(320\text{--}1280)$",
                title="Collapse-time estimate drifts with fit band")
    axes[1].set(xlabel=r"geometric-mean field in fit band",
                ylabel=r"free-fit tail exponent $\alpha_{\rm fit}$",
                title="The apparent exponent drifts too")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_tc_free(fits, cases, output):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharex=True)
    variables = (
        ("alpha_from_velocity", r"$\alpha$ from $|\partial_t\phi_0|$ vs $|\phi_0|$"),
        ("q_beta_over_alpha", r"$q=\beta/\alpha$ from $R_{50}$ vs $|\phi_0|$"),
        ("alpha_plus_beta_tc_free", r"$\alpha+\beta$ without a fitted $t_c$"),
    )
    colors = plt.cm.viridis(np.linspace(0.05, 0.9, len(cases)))
    for case, color in zip(cases, colors):
        selected = sorted((row for row in fits if row["case"] == case["name"]),
                          key=lambda row: np.sqrt(row["phi_min"]*row["phi_max"]))
        if not selected:
            continue
        x = [np.sqrt(row["phi_min"]*row["phi_max"]) for row in selected]
        label = f"N={case['n']}, +{case['offset']:g}%"
        if case["field_step"] != 0.15:
            label += f", dt cap={case['field_step']:g}"
        for ax, (key, _) in zip(axes, variables):
            ax.plot(x, [row[key] for row in selected], "o-", color=color,
                    label=label)
    axes[2].axhline(1.0, color="black", ls="--", lw=1, label="target 1")
    axes[2].axhline(1.5, color="black", ls=":", lw=1, label="late trend 3/2")
    for ax, (_, ylabel) in zip(axes, variables):
        ax.set_xscale("log")
        ax.set_xlabel(r"geometric-mean central field in fixed fit band")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.25)
    axes[0].legend(fontsize=7, loc="best")
    axes[2].legend(fontsize=7, loc="best")
    fig.suptitle("Parametric scaling check that does not fit the collapse time")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_profiles(profiles, cases, output):
    fig, axes = plt.subplots(1, len(cases), figsize=(5.2*len(cases), 4.6),
                             squeeze=False)
    levels = sorted({row["level"] for row in profiles})
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(levels)))
    for ax, case in zip(axes[0], cases):
        for level, color in zip(levels, colors):
            rows = sorted((row for row in profiles
                           if row["case"] == case["name"]
                           and row["level"] == level),
                          key=lambda row: row["rho_r_over_r50"])
            if rows:
                ax.plot([row["rho_r_over_r50"] for row in rows],
                        [row["normalized_abs_phi"] for row in rows],
                        color=color, label=fr"$|\phi_0|={level:g}$")
        rho = np.linspace(0, 2, 101)
        ax.plot(rho, 1/(1+rho**2), "k--", label=r"$1/(1+\rho^2)$")
        ax.plot(rho, np.exp(-np.log(2)*rho**2), "k:", label="Gaussian")
        ax.set(xlabel=r"$r/R_{50}$", ylabel=r"$|\phi|/|\phi_0|$",
               title=f"N={case['n']}, A offset +{case['offset']:g}%",
               xlim=(0, 2), ylim=(0, 1.03))
        ax.grid(True, alpha=0.25)
    axes[0][0].legend(fontsize=7)
    fig.suptitle("Near-transition core profiles, normalized by R50")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=audit.parse_case,
                        help="N:offset_percent[:CFL[:field_step]], repeatable")
    parser.add_argument("--max-level", type=float, default=1280.0)
    parser.add_argument("--tag", default="",
                        help="optional subfolder for this independent scan batch")
    parser.add_argument("--profile-level", action="append", type=float,
                        help="also save normalized radial profiles at this center field")
    args = parser.parse_args()
    cases = args.case or [audit.parse_case(value) for value in (
        "1000:0.01", "1000:0.029990845", "1000:0.03",
        "1000:0.03125", "1000:0.0325", "1000:0.05",
        "1000:0.03125:0.4:0.05",
    )]
    if args.max_level < 1280:
        parser.error("max-level must be at least 1280 for the declared fit bands")

    profile_levels = args.profile_level or []
    if any(level <= 0 or level > args.max_level for level in profile_levels):
        parser.error("profile levels must be positive and no larger than max-level")
    audit.LEVELS = tuple(np.unique(np.r_[
        np.geomspace(4.0, args.max_level, 121),
        [x for x in (10, 12, 20, 30, 40, 80, 160, 320, 640, 1280)
         if x <= args.max_level],
        profile_levels,
    ]))
    solver = audit.load_solver()
    all_events, all_profiles, settings, tc_rows, fit_rows, tau_rows, tc_free_rows = (
        [], [], [], [], [], [], []
    )
    for case in cases:
        events, steps, final_time, profiles = audit.simulate(
            solver, case, args.max_level, capture_profile_levels=profile_levels,
            verbose=False,
        )
        all_events.extend(events)
        all_profiles.extend(profiles)
        settings.append({**case, "steps": steps, "final_time": final_time,
                         "highest_level_reached": max(
                             (row["level"] for row in events), default=""),
                         "event_count": len(events)})
        print(f"{case['name']}: {len(events)} thresholds through "
              f"{max((row['level'] for row in events), default=0):g}; "
              f"{steps} steps, t={final_time:.6f}", flush=True)
        tc_free_rows.extend(fit_tc_free(events, case))
        for tc_band in TC_BANDS:
            fit = audit.fit_powerlaw(events, *tc_band, "abs_phi_center")
            if not fit:
                continue
            tc, alpha_tail, r2, n = fit
            tc_rows.append({"case": case["name"], "N": case["n"],
                            "offset_percent": case["offset"],
                            "field_step": case["field_step"],
                            "fit_level_min": tc_band[0],
                            "fit_level_max": tc_band[1], "tc": tc,
                            "alpha_tail": alpha_tail, "r2": r2, "points": n})
            fit_rows.extend(fit_local_exponents(events, case["name"], tc, tc_band))
            tau_rows.extend(fit_tau_windows(events, case["name"], tc, tc_band))

    output = RESULTS / args.tag if args.tag else RESULTS
    output.mkdir(parents=True, exist_ok=True)
    audit.write_csv(output/"regime_events.csv", all_events)
    audit.write_csv(output/"tc_sensitivity.csv", tc_rows)
    audit.write_csv(output/"local_exponent_fits.csv", fit_rows)
    audit.write_csv(output/"tau_window_fits.csv", tau_rows)
    audit.write_csv(output/"tc_free_fits.csv", tc_free_rows)
    audit.write_csv(output/"run_settings.csv", settings)
    if all_profiles:
        audit.write_csv(output/"profile_snapshots.csv", all_profiles)
        audit.write_csv(output/"profile_shape_fits.csv",
                        audit.fit_profile_shapes(all_profiles))
        plot_profiles(all_profiles, cases, output/"profile_evolution.png")
    if tc_free_rows:
        plot_tc_free(tc_free_rows, cases, output/"tc_free_scaling.png")
    if any(row["tc_fit_level_min"] == 160 for row in fit_rows):
        plot_regimes(fit_rows, cases, output/"exponent_drift.png")
        plot_tau_windows(tau_rows, cases, output/"sum_vs_tau_window.png")
    if any(row["fit_level_min"] == 320 and row["fit_level_max"] == 1280
           for row in tc_rows):
        plot_tc_drift(tc_rows, cases, output/"tc_fit_band_drift.png")
    if not tc_rows:
        print("No field thresholds were reached; wrote only run settings.", flush=True)
    elif not any(row["tc_fit_level_min"] == 160 for row in fit_rows):
        print("No fit spans |phi0|=160--1280; skipped high-field plots.", flush=True)
    print(f"Wrote regime scan to {output}", flush=True)


if __name__ == "__main__":
    main()
