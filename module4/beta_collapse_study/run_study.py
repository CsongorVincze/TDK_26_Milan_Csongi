#!/usr/bin/env python3
"""Compare Module 4 collapse-width observables and beta-fit sensitivities."""

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MODULE4 = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(MODULE4))
from milan_alapkod_module4 import Config, make_grid, rk4  # noqa: E402


OUT = Path(__file__).resolve().parent / "results"
CRITICAL_SUMMARY = MODULE4 / "results" / "module4_critical_summary.json"
LEVEL = 20.0
SAMPLE_INTERVAL = 0.01
RADIUS_CUTOFFS = (4.0, 8.0, 12.0)
BASE_OBSERVABLES = (
    "rms_signed_R4", "rms_signed_R8", "rms_signed_R12",
    "rms_absrho_R8", "rms_positive_rho_R8", "rms_gradient_R8",
    "field_r50", "field_r20",
)
WINDOWS = ((0.10, 0.30), (0.20, 0.60), (0.40, 1.00), (0.80, 1.60))
TC_BANDS = ((2.0, 4.0), (4.0, 8.0), (8.0, 12.0),
            (12.0, LEVEL), (20.0, 30.0), (30.0, 40.0),
            (40.0, 60.0), (60.0, 80.0))
CONTINUATION_LEVELS = (30.0, 40.0, 60.0, 80.0)


def write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def weighted_rms(weight, radius):
    den = np.trapezoid(weight*radius**2, x=radius)
    scale = np.trapezoid(np.abs(weight)*radius**2, x=radius)
    if abs(den) <= 1e-12*max(scale, 1.0):
        return float("nan"), float("nan"), float(den)
    mean_r2 = np.trapezoid(weight*radius**4, x=radius)/den
    rms = np.sqrt(mean_r2) if np.isfinite(mean_r2) and mean_r2 > 0 else float("nan")
    return float(mean_r2), float(rms), float(den)


def field_width(phi, radius, fraction, cutoff):
    use = radius <= cutoff
    r, value = radius[use], np.abs(phi[use])
    if len(r) < 2 or value[0] <= 0:
        return float("nan")
    crossings = np.flatnonzero(value[1:] <= fraction*value[0])
    if not len(crossings):
        return float("nan")
    i = int(crossings[0]+1)
    y0, y1 = value[i-1], value[i]
    if y0 == y1:
        return float(r[i])
    f = np.clip((fraction*value[0]-y0)/(y1-y0), 0.0, 1.0)
    return float(r[i-1]+f*(r[i]-r[i-1]))


def measure(t, phi, pi, config, grid):
    dphi = np.gradient(phi, grid.x, edge_order=2)/grid.rx
    dphi[0] = 0.0
    rho = 0.5*pi**2 + 0.5*dphi**2 + 0.5*phi**2 - 0.25*phi**4
    row = {"time": float(t), "phi_center": float(phi[0]),
           "abs_phi_center": float(abs(phi[0]))}
    moment_values = {}
    for cutoff in RADIUS_CUTOFFS:
        mask = grid.r <= cutoff
        r = grid.r[mask]
        mean_r2, rms, den = weighted_rms(rho[mask], r)
        moment_values[cutoff] = (mean_r2, rms, den)
        row[f"mean_r2_signed_R{cutoff:g}"] = mean_r2
        row[f"rms_signed_R{cutoff:g}"] = rms
        row[f"signed_energy_den_R{cutoff:g}"] = den
    mask = grid.r <= config.core_radius
    r = grid.r[mask]
    _, row["rms_absrho_R8"], abs_den = weighted_rms(np.abs(rho[mask]), r)
    _, row["rms_positive_rho_R8"], _ = weighted_rms(np.maximum(rho[mask], 0), r)
    _, row["rms_gradient_R8"], _ = weighted_rms(dphi[mask]**2, r)
    _, _, signed_den = moment_values[config.core_radius]
    row["signed_cancellation_ratio_R8"] = (
        signed_den/abs_den if abs_den > 0 else float("nan")
    )
    row["field_r50"] = field_width(phi, grid.r, 0.50, config.core_radius)
    row["field_r20"] = field_width(phi, grid.r, 0.20, config.core_radius)
    row["total_energy"] = float(4*np.pi*np.trapezoid(rho*grid.r**2, x=grid.r))
    return row


def evolve(amplitude, n, t_max, sample_interval, level=LEVEL):
    config = Config(n=n, collapse_level=level, sample_interval=sample_interval)
    grid = make_grid(config)
    phi = amplitude*np.exp(-grid.r**2/(2*config.r0**2))
    pi = np.zeros_like(phi)
    phi[-1] = pi[-1] = 0.0
    t, steps, next_sample = 0.0, 0, sample_interval
    rows = [measure(t, phi, pi, config, grid)]
    cutoff_reached = abs(phi[0]) >= level
    numerical_failure = False
    while t < t_max and not cutoff_reached:
        peak = max(float(np.max(np.abs(phi))), 1.0)
        dt = min(grid.dt, 0.15/peak, t_max-t, max(next_sample-t, 1e-12))
        next_phi, next_pi = rk4(phi, pi, dt, grid)
        if not np.all(np.isfinite(next_phi)) or not np.all(np.isfinite(next_pi)):
            numerical_failure = True
            break
        cutoff_reached = abs(next_phi[0]) >= level
        if cutoff_reached:
            rise = abs(next_phi[0])-abs(phi[0])
            frac = np.clip((level-abs(phi[0]))/rise, 0.0, 1.0)
            t += dt*frac
            phi += frac*(next_phi-phi)
            pi += frac*(next_pi-pi)
            phi[0] = np.copysign(level, next_phi[0])
        else:
            t += dt
            phi, pi = next_phi, next_pi
        steps += 1
        if t >= next_sample-1e-11 or cutoff_reached:
            rows.append(measure(t, phi, pi, config, grid))
            while next_sample <= t+1e-11:
                next_sample += sample_interval
    if rows[-1]["time"] < t-1e-11:
        rows.append(measure(t, phi, pi, config, grid))
    return rows, config, grid, steps, cutoff_reached, numerical_failure, t


def line_fit(x, y):
    if len(x) < 5 or np.ptp(x) == 0:
        return float("nan"), float("nan"), float("nan")
    slope, intercept = np.polyfit(x, y, 1)
    total = np.sum((y-y.mean())**2)
    residual = np.sum((y-(slope*x+intercept))**2)
    r2 = 1-residual/total if total > 0 else float("nan")
    return float(slope), float(intercept), float(r2)


def estimate_tc(rows, band):
    time = np.asarray([r["time"] for r in rows])
    center = np.asarray([r["abs_phi_center"] for r in rows])
    lo, hi = band
    ids = np.flatnonzero((center >= lo) & (center <= hi))[-40:]
    if len(ids) < 6:
        return None
    slope, intercept, r2 = line_fit(time[ids], 1/center[ids])
    tc = -intercept/slope if slope != 0 else float("nan")
    if not np.isfinite(tc) or slope >= 0 or tc <= time[ids[-1]]:
        return None
    return {"tc": tc, "r2": r2, "n": len(ids),
            "band_low": lo, "band_high": hi, "gamma": 1.0}


def estimate_tc_powerlaw(rows, count):
    """Jointly fit |phi(0,t)| ~ (tc-t)^(-gamma), without fixing gamma=1."""
    time = np.asarray([r["time"] for r in rows])
    center = np.asarray([r["abs_phi_center"] for r in rows])
    ids = np.flatnonzero(center >= 4.0)[-count:]
    if len(ids) < 8:
        return None
    fit_time, log_center = time[ids], np.log(center[ids])
    left, right = time[-1]+1e-5, time[-1]+1.0
    best = None
    for _ in range(4):
        trial_tc = np.linspace(left, right, 1201)
        log_tau = np.log(trial_tc[:, None]-fit_time[None, :])
        centered_x = log_tau-log_tau.mean(axis=1, keepdims=True)
        centered_y = log_center-log_center.mean()
        denom = np.sum(centered_x**2, axis=1)
        slopes = np.sum(centered_x*centered_y, axis=1)/denom
        intercepts = log_center.mean()-slopes*log_tau.mean(axis=1)
        residual = log_center[None, :] - (
            slopes[:, None]*log_tau + intercepts[:, None]
        )
        rss = np.sum(residual**2, axis=1)
        rss[slopes >= 0] = np.inf
        index = int(np.argmin(rss))
        step = (right-left)/(len(trial_tc)-1)
        best = (float(trial_tc[index]), float(slopes[index]),
                float(intercepts[index]), float(rss[index]), step)
        left = max(time[-1]+1e-6, best[0]-5*step)
        right = best[0]+5*step
    tc, slope, intercept, rss, _ = best
    total = np.sum((log_center-log_center.mean())**2)
    return {"tc": tc, "r2": 1-rss/total if total > 0 else float("nan"),
            "n": len(ids), "band_low": 4.0, "band_high": float(center[-1]),
            "gamma": -slope}


def tc_candidates(rows):
    final_level = rows[-1]["abs_phi_center"]
    candidates = []
    for low, high in TC_BANDS:
        use_high = min(high, final_level+1e-8)
        fit = estimate_tc(rows, (low, use_high)) if use_high > low else None
        if fit:
            fit["source"] = f"field_{low:g}_{high:g}"
            candidates.append(fit)
    use = np.flatnonzero(
        np.asarray([r["abs_phi_center"] for r in rows]) >= max(2.0, 0.2*LEVEL)
    )[-40:]
    if len(use) >= 6:
        time = np.asarray([r["time"] for r in rows])
        center = np.asarray([r["abs_phi_center"] for r in rows])
        slope, intercept, r2 = line_fit(time[use], 1/center[use])
        tc = -intercept/slope if slope != 0 else float("nan")
        if np.isfinite(tc) and slope < 0 and tc > time[use[-1]]:
            candidates.append({"tc": tc, "r2": r2, "n": len(use),
                               "band_low": max(2.0, 0.2*LEVEL),
                               "band_high": final_level,
                               "source": "current_last40_rule", "gamma": 1.0})
    for count in (20, 40, 80):
        fit = estimate_tc_powerlaw(rows, count)
        if fit:
            fit["source"] = f"joint_powerlaw_{count}"
            candidates.append(fit)
    return candidates


def beta_fits(case, tc_fits):
    rows = case["rows"]
    time = np.asarray([r["time"] for r in rows])
    observables = BASE_OBSERVABLES + ("mean_r2_signed_R8",)
    result = []
    for tc_fit in tc_fits:
        tau = tc_fit["tc"]-time
        for lower, upper in WINDOWS:
            temporal = (tau >= lower) & (tau <= upper)
            for name in observables:
                value = np.asarray([r[name] for r in rows])
                use = temporal & np.isfinite(value) & (value > 0)
                x, y = np.log10(tau[use]), np.log10(value[use])
                slope, intercept, r2 = line_fit(x, y)
                result.append({
                    "case": case["case"], "N": case["N"],
                    "amplitude": case["amplitude"], "observable": name,
                    "tc_source": tc_fit["source"], "tc": tc_fit["tc"],
                    "tc_fit_r2": tc_fit["r2"], "tau_min": lower,
                    "tau_max": upper, "n": int(use.sum()),
                    "gamma_central_fit": tc_fit.get("gamma", ""),
                    "slope_log_observable": slope,
                    "beta_radius_equivalent": slope/2 if name == "mean_r2_signed_R8" else slope,
                    "fit_r2": r2,
                })
    return result


def amplitude_scaling_fits(case):
    """Fit width ~ |phi(0)|^(-beta/gamma) on a common central-field tail."""
    rows = case["rows"]
    center = np.asarray([r["abs_phi_center"] for r in rows])
    ids = np.flatnonzero(center >= 4.0)[-40:]
    if len(ids) < 8:
        return []
    tc_fit = estimate_tc_powerlaw(rows, 40)
    gamma = tc_fit["gamma"] if tc_fit else float("nan")
    result = []
    for name in BASE_OBSERVABLES + ("mean_r2_signed_R8",):
        value = np.asarray([r[name] for r in rows])
        use = np.isfinite(value[ids]) & (value[ids] > 0)
        x, y = np.log10(center[ids][use]), np.log10(value[ids][use])
        slope, intercept, r2 = line_fit(x, y)
        result.append({
            "case": case["case"], "N": case["N"],
            "amplitude": case["amplitude"], "observable": name,
            "central_field_min": float(center[ids][0]),
            "central_field_max": float(center[ids[-1]]),
            "n": int(use.sum()), "central_gamma_joint_fit": gamma,
            "slope_log_observable_vs_log_phi": slope,
            "beta_over_gamma": (-slope/2 if name == "mean_r2_signed_R8" else -slope),
            "beta_if_joint_gamma": ((-slope/2 if name == "mean_r2_signed_R8" else -slope)*gamma
                                    if np.isfinite(gamma) else float("nan")),
            "fit_r2": r2,
        })
    return result


def save_run_csv(case, out):
    columns = list(case["rows"][0].keys())
    write_csv(out, columns, ([row.get(key, "") for key in columns]
                             for row in case["rows"]))


def plot_observable_evolution(cases, out):
    focus = [case for case in cases if case["N"] == 500 and "offset" in case
             and "_phi" not in case["case"]]
    fig, axes = plt.subplots(len(focus), 1, figsize=(10, 3.4*len(focus)),
                             sharex=False, squeeze=False)
    axes = axes[:, 0]
    colors = plt.get_cmap("tab10").colors
    for ax, case in zip(axes, focus):
        rows = case["rows"]
        time = np.asarray([r["time"] for r in rows])
        center = np.asarray([r["abs_phi_center"] for r in rows])
        start_index = int(np.flatnonzero(center >= 2.0)[0]) if np.any(center >= 2) else 0
        event_time = case["event_time"] or time[-1]
        use = (time >= max(time[start_index], event_time-2.0)) & (time <= event_time)
        for color, name in zip(colors, BASE_OBSERVABLES):
            values = np.asarray([r[name] for r in rows])
            base_indices = np.flatnonzero(use & np.isfinite(values) & (values > 0))
            if not len(base_indices):
                continue
            base = values[base_indices[0]]
            ax.plot(time[use], values[use]/base, color=color, lw=1.4, label=name)
        ax.axhline(1, color="black", lw=0.8, alpha=0.5)
        ax.set(ylabel="width / width at window start",
               title=f"N=500; A0={case['amplitude']:.9f} "
                     f"(+{100*case['offset']:g}% over saved Acrit)")
        ax.grid(alpha=0.25)
        ax.legend(ncol=2, fontsize=7.5, loc="best")
    axes[-1].set_xlabel("time t")
    fig.suptitle("Different core-size diagnostics near the phi=20 cutoff", y=1.0)
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_beta_sensitivity(fits, out):
    names = list(BASE_OBSERVABLES)
    fig, axes = plt.subplots(4, 2, figsize=(12, 13), sharex=True)
    axes = axes.flat
    cases = sorted({row["case"] for row in fits
                    if row["tc_source"] == "current_last40_rule" and row["N"] == 500
                    and "_phi" not in row["case"]})
    for ax, name in zip(axes, names):
        for case_name in cases:
            vals = [row for row in fits if row["case"] == case_name
                    and row["observable"] == name
                    and row["tc_source"] == "current_last40_rule"]
            vals.sort(key=lambda row: row["tau_min"])
            if vals:
                label = f"{vals[0]['amplitude']:.6f}"
                ax.plot(range(len(vals)), [row["beta_radius_equivalent"] for row in vals],
                        marker="o", ms=3.5, lw=1, label=label)
        ax.axhline(0, color="black", lw=0.8, alpha=0.6)
        ax.set_title(name, fontsize=9)
        ax.set_xticks(range(len(WINDOWS)),
                      labels=[f"{lo:g}–{hi:g}" for lo, hi in WINDOWS])
        ax.grid(alpha=0.25)
    axes[0].legend(title="A0", fontsize=7, title_fontsize=8)
    fig.supxlabel(r"fit window in remaining time $\tau=t_c-t$")
    fig.supylabel(r"radius-equivalent slope $\beta$ (for raw $\langle r^2\rangle$, slope/2)")
    fig.suptitle("Beta diagnostics depend on observable and fit window")
    fig.tight_layout(rect=(0.05, 0.03, 1, 0.97))
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_tc_sensitivity(tc_rows, out):
    fig, ax = plt.subplots(figsize=(9, 4.8))
    runs = sorted({row["case"] for row in tc_rows})
    source_order = {f"field_{lo:g}_{hi:g}": i
                    for i, (lo, hi) in enumerate(TC_BANDS)}
    source_order["current_last40_rule"] = len(source_order)
    for count in (20, 40, 80):
        source_order[f"joint_powerlaw_{count}"] = len(source_order)+count
    for case in runs:
        values = [row for row in tc_rows if row["case"] == case]
        values.sort(key=lambda row: source_order.get(row["source"], 99))
        if values:
            if "_plus_" in case:
                suffix = (f", cutoff {case.rsplit('_phi', 1)[1]}"
                          if "_phi" in case else "")
                label = (f"N{values[0]['N']} +{100*values[0]['offset']:g}%"
                         f"{suffix}")
            else:
                label = f"N{values[0]['N']}, fixed A0"
            ax.plot([row["source"] for row in values], [row["tc"] for row in values],
                    marker="o", ms=4, lw=1, label=label)
    ax.set(ylabel=r"estimated collapse time $t_c$",
           title=r"$t_c$ sensitivity to estimator and field-amplitude window")
    ax.tick_params(axis="x", rotation=20)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_representative_loglog(case, tc_fit, out):
    rows = case["rows"]
    time = np.asarray([row["time"] for row in rows])
    tau = tc_fit["tc"]-time
    axes_names = ("rms_signed_R8", "rms_absrho_R8",
                  "rms_gradient_R8", "field_r50")
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for ax, name in zip(axes.flat, axes_names):
        values = np.asarray([row[name] for row in rows])
        use = (tau >= 0.2) & (tau <= 0.6) & np.isfinite(values) & (values > 0)
        x, y = np.log10(tau[use]), np.log10(values[use])
        slope, intercept, r2 = line_fit(x, y)
        ax.plot(x, y, "o", ms=3.5, color="#2563eb", label="samples")
        if np.isfinite(slope):
            line_x = np.linspace(x.min(), x.max(), 80)
            ax.plot(line_x, slope*line_x+intercept, "--", color="#dc2626",
                    label="OLS fit")
        ax.set(title=f"{name}; slope={slope:.3f}, R²={r2:.3f}, n={len(x)}",
               xlabel=r"$\log_{10}(t_c-t)$", ylabel=rf"$\log_{{10}}({name})$")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle(
        f"Explicit log-log fits: {case['case']}; {tc_fit['source']} "
        + rf"($t_c={tc_fit['tc']:.5f}$), $\tau=0.2$–$0.6$"
    )
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_beta_tc_methods(beta_rows, out):
    focus = [row for row in beta_rows
             if row["case"] == "N500_plus_0.01pct"
             and row["observable"] in ("rms_signed_R8", "field_r50")
             and row["tau_min"] == 0.2 and row["tau_max"] == 0.6
             and int(row["n"]) >= 8 and np.isfinite(row["beta_radius_equivalent"])]
    sources = list(dict.fromkeys(row["tc_source"] for row in focus))
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    for ax, name in zip(axes, ("rms_signed_R8", "field_r50")):
        rows = [row for row in focus if row["observable"] == name]
        rows.sort(key=lambda row: sources.index(row["tc_source"]))
        ax.plot([row["tc_source"] for row in rows],
                [row["beta_radius_equivalent"] for row in rows],
                "o-", color="#2563eb")
        ax.axhline(0, color="black", lw=0.8, alpha=0.6)
        ax.set(ylabel=r"radius-equivalent $\beta$", title=name)
        ax.grid(alpha=0.25)
    axes[-1].tick_params(axis="x", rotation=20)
    axes[-1].set_xlabel(r"collapse-time estimator; common fit window $0.2\leq t_c-t\leq0.6$")
    fig.suptitle("Beta sensitivity to collapse-time estimation (N=500, +0.01%)")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_amplitude_scaling(case, out):
    rows = case["rows"]
    center = np.asarray([row["abs_phi_center"] for row in rows])
    displays = (
        ("rms_signed_R8", "signed-energy RMS radius"),
        ("rms_absrho_R8", "absolute-energy RMS radius"),
        ("rms_gradient_R8", "gradient-energy RMS radius"),
        ("field_r50", "field half-width"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for ax, (name, display) in zip(axes.flat, displays):
        values = np.asarray([row[name] for row in rows])
        candidates = np.flatnonzero(center >= 4.0)[-40:]
        use = np.zeros(len(center), dtype=bool)
        use[candidates] = True
        use &= np.isfinite(values) & (values > 0)
        x, y = np.log10(center[use]), np.log10(values[use])
        slope, intercept, r2 = line_fit(x, y)
        ax.plot(x, y, "o", ms=3.5, color="#2563eb", label="samples")
        if np.isfinite(slope):
            line_x = np.linspace(x.min(), x.max(), 80)
            ax.plot(line_x, slope*line_x+intercept, "--", color="#dc2626",
                    label="OLS fit")
        q = -slope/2 if name == "mean_r2_signed_R8" else -slope
        ax.set_title(f"{display}\nq={q:.3f}, R²={r2:.3f}, n={len(x)}",
                     fontsize=9)
        ax.set(xlabel=r"$\log_{10}|\phi(0,t)|$",
               ylabel=rf"$\log_{{10}}({name})$")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle("Direct width-vs-field scaling (no collapse-time estimate)")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_cutoff_continuation(cases, out):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    colors = {20: "#2563eb", 30: "#dc2626", 40: "#059669",
              60: "#9333ea", 80: "#ea580c"}
    line_styles = {500: "-", 1000: "--"}
    displays = (("rms_signed_R8", "signed-energy RMS radius"),
                ("field_r50", "field half-width"))
    for ax, (name, title) in zip(axes, displays):
        for case in cases:
            rows = case["rows"]
            center = np.asarray([row["abs_phi_center"] for row in rows])
            value = np.asarray([row[name] for row in rows])
            use = (center >= 4.0) & np.isfinite(value) & (value > 0)
            ax.plot(np.log10(center[use]), np.log10(value[use]),
                    marker=".", ls=line_styles.get(case["N"], ":"),
                    ms=3, lw=1, color=colors[case["collapse_level"]],
                    label=f"N={case['N']}; stop at |phi(0)|={case['collapse_level']}")
        ax.axvline(np.log10(LEVEL), color="black", ls=":", lw=0.9,
                   label="previous cutoff 20")
        ax.set(title=title, xlabel=r"$\log_{10}|\phi(0,t)|$",
               ylabel=rf"$\log_{{10}}({name})$")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7.5)
    fig.suptitle("Same near-critical run continued to progressively higher field cutoffs")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_resolution(cases, out):
    focus = [case for case in cases if case["resolution_test"]
             and case.get("collapse_level", LEVEL) == LEVEL]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    colors = {250: "#2563eb", 500: "#dc2626", 1000: "#059669"}
    for case in focus:
        rows = case["rows"]
        time = np.asarray([r["time"] for r in rows])
        center = np.asarray([r["abs_phi_center"] for r in rows])
        axes[0].plot(time, center, color=colors.get(case["N"]),
                     label=f"N={case['N']} ({case['status']})")
        radius = np.asarray([r["rms_signed_R8"] for r in rows])
        good = np.isfinite(radius) & (radius > 0)
        if np.any(good):
            base = radius[np.flatnonzero(good)[0]]
            axes[1].plot(time[good], radius[good]/base,
                         color=colors.get(case["N"]), label=f"N={case['N']}")
    axes[0].axhline(LEVEL, color="black", ls=":", lw=0.9, label="stop level")
    axes[0].set(xlabel="time t", ylabel=r"$|\phi(0,t)|$",
                title="Same A0 across spatial resolutions")
    axes[1].set(xlabel="time t", ylabel="signed RMS radius / initial value",
                title=r"Existing signed-energy $R_{core}=8$ observable")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t-max", type=float, default=100.0)
    parser.add_argument("--sample-interval", type=float, default=SAMPLE_INTERVAL)
    parser.add_argument("--resolutions", type=int, nargs="+", default=(250, 500, 1000))
    args = parser.parse_args()
    if not CRITICAL_SUMMARY.exists():
        parser.error(f"Missing saved critical-amplitude summary: {CRITICAL_SUMMARY}")
    if args.t_max <= 0 or args.sample_interval <= 0 or any(n < 20 for n in args.resolutions):
        parser.error("Use positive times and resolutions N >= 20.")

    acrit = float(json.loads(CRITICAL_SUMMARY.read_text(encoding="utf-8"))["Acrit_estimate"])
    offsets = (1e-4, 5e-4, 1e-3)
    study_cases = []
    baseline = None
    for offset in offsets:
        amplitude = acrit*(1+offset)
        name = f"N500_plus_{100*offset:g}pct"
        print(f"Running {name}: A0={amplitude:.12f}", flush=True)
        rows, config, grid, steps, crossed, failed, end_time = evolve(
            amplitude, 500, args.t_max, args.sample_interval,
        )
        case = {"case": name, "N": 500, "amplitude": amplitude,
                "offset": offset, "rows": rows, "steps": steps,
                "event_time": end_time if crossed else None,
                "status": ("cutoff20" if crossed else "no cutoff by t_max"),
                "numerical_failure": failed, "resolution_test": False,
                "collapse_level": int(LEVEL)}
        save_run_csv(case, OUT/f"{name}.csv")
        study_cases.append(case)
        if np.isclose(offset, 1e-4):
            baseline = case

    continuation_cases = []
    continuation_specs = [(500, level) for level in CONTINUATION_LEVELS]
    continuation_specs.append((1000, 80.0))
    for n, level in continuation_specs:
        name = f"N{n}_plus_0.01pct_phi{int(level)}"
        print(f"Continuation {name}: same A0, stop at |phi(0)|={level:g}",
              flush=True)
        rows, config, grid, steps, crossed, failed, end_time = evolve(
            baseline["amplitude"], n, args.t_max, args.sample_interval, level,
        )
        case = {"case": name, "N": n, "amplitude": baseline["amplitude"],
                "offset": 1e-4, "rows": rows, "steps": steps,
                "event_time": end_time if crossed else None,
                "status": (f"cutoff{level:g}" if crossed else "no cutoff by t_max"),
                "numerical_failure": failed, "resolution_test": (n != 500),
                "collapse_level": int(level)}
        save_run_csv(case, OUT/f"{name}.csv")
        study_cases.append(case)
        continuation_cases.append(case)

    fixed_amplitude = acrit*(1+1e-4)
    for n in args.resolutions:
        if n == 500:
            case = dict(baseline)
            case["case"] = "N500_fixed_A_resolution_reference"
            case["resolution_reference_duplicate"] = True
        else:
            name = f"N{n}_fixed_A{fixed_amplitude:.9f}"
            print(f"Resolution check {name}", flush=True)
            rows, config, grid, steps, crossed, failed, end_time = evolve(
                fixed_amplitude, n, args.t_max, args.sample_interval,
            )
            case = {"case": name, "N": n, "amplitude": fixed_amplitude,
                    "offset": 1e-4, "rows": rows, "steps": steps,
                    "event_time": end_time if crossed else None,
                    "status": ("cutoff20" if crossed else "no cutoff by t_max"),
                    "numerical_failure": failed, "resolution_test": True}
            save_run_csv(case, OUT/f"{name}.csv")
        case["resolution_test"] = True
        study_cases.append(case)

    all_tc, all_beta, all_amplitude_scaling, run_summary = [], [], [], []
    analysis_cases = [case for case in study_cases
                      if not case.get("resolution_reference_duplicate")]
    for case in analysis_cases:
        fits = tc_candidates(case["rows"])
        for fit in fits:
            all_tc.append({"case": case["case"], "N": case["N"],
                           "amplitude": case["amplitude"],
                           "offset": case["offset"], **fit})
        beta = beta_fits(case, fits)
        all_beta.extend(beta)
        all_amplitude_scaling.extend(amplitude_scaling_fits(case))
        max_center_row = max(case["rows"], key=lambda row: row["abs_phi_center"])
        event_row = case["rows"][-1]
        run_summary.append({
            "case": case["case"], "N": case["N"],
            "amplitude": case["amplitude"],
            "collapse_level": case.get("collapse_level", LEVEL),
            "offset_percent_vs_N500_Acrit": 100*(case["amplitude"]/acrit-1),
            "status": case["status"], "event_or_end_time": event_row["time"],
            "max_abs_phi_center": max_center_row["abs_phi_center"],
            "max_abs_phi_time": max_center_row["time"],
            "rms_signed_R8_at_end": event_row["rms_signed_R8"],
            "mean_r2_signed_R8_at_end": event_row["mean_r2_signed_R8"],
            "steps": case["steps"], "numerical_failure": case["numerical_failure"],
        })
    write_csv(OUT/"run_summary.csv", list(run_summary[0]),
              ([row.get(key, "") for key in run_summary[0]] for row in run_summary))
    if all_tc:
        write_csv(OUT/"tc_sensitivity.csv", list(all_tc[0]),
                  ([row.get(key, "") for key in all_tc[0]] for row in all_tc))
    if all_beta:
        write_csv(OUT/"beta_fit_sensitivity.csv", list(all_beta[0]),
                  ([row.get(key, "") for key in all_beta[0]] for row in all_beta))
    if all_amplitude_scaling:
        write_csv(OUT/"amplitude_scaling.csv", list(all_amplitude_scaling[0]),
                  ([row.get(key, "") for key in all_amplitude_scaling[0]]
                   for row in all_amplitude_scaling))
    plot_observable_evolution(study_cases[:len(offsets)], OUT/"observable_evolution.png")
    plot_beta_sensitivity(all_beta, OUT/"beta_fit_sensitivity.png")
    plot_tc_sensitivity(all_tc, OUT/"tc_sensitivity.png")
    plot_resolution(study_cases, OUT/"resolution_sensitivity.png")
    base_case = next(case for case in analysis_cases
                     if case["case"] == "N500_plus_0.01pct")
    base_tc = next(fit for fit in tc_candidates(base_case["rows"])
                   if fit["source"] == "current_last40_rule")
    plot_representative_loglog(
        base_case, base_tc, OUT/"representative_loglog_fits.png",
    )
    base_joint_tc = next(fit for fit in tc_candidates(base_case["rows"])
                         if fit["source"] == "joint_powerlaw_40")
    plot_representative_loglog(
        base_case, base_joint_tc, OUT/"representative_loglog_fits_joint_tc.png",
    )
    plot_beta_tc_methods(all_beta, OUT/"beta_tc_method_sensitivity.png")
    plot_amplitude_scaling(
        base_case, OUT/"amplitude_scaling.png",
    )
    plot_cutoff_continuation(
        [baseline] + continuation_cases,
        OUT/"cutoff_continuation.png",
    )

    metadata = {
        "Acrit_from_saved_N500_search": acrit, "base_collapse_cutoff": LEVEL,
        "continuation_runs": [{"N": n, "cutoff": int(level)}
                              for n, level in continuation_specs],
        "fit_formula": "log10(radius) = beta*log10(tc - t) + constant",
        "sample_interval": args.sample_interval, "horizon": args.t_max,
        "offsets_for_N500": offsets, "resolutions": list(args.resolutions),
        "radius_observables": list(BASE_OBSERVABLES),
        "fit_tau_windows": WINDOWS, "tc_field_bands": TC_BANDS,
        "outputs": ["run_summary.csv", "tc_sensitivity.csv",
                    "beta_fit_sensitivity.csv", "observable_evolution.png",
                    "amplitude_scaling.csv", "amplitude_scaling.png",
                    "cutoff_continuation.png",
                    "beta_fit_sensitivity.png", "tc_sensitivity.png",
                    "resolution_sensitivity.png", "representative_loglog_fits.png",
                    "representative_loglog_fits_joint_tc.png",
                    "beta_tc_method_sensitivity.png"],
    }
    (OUT/"study_settings.json").write_text(
        json.dumps(metadata, indent=2)+"\n", encoding="utf-8",
    )
    print(f"Saved study data and figures under {OUT}")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    main()
