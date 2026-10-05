#!/usr/bin/env python3
"""Fourier-filter the alapkod central-field decay and fit its envelope."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


FOLDER = Path(__file__).resolve().parent
SOURCE = FOLDER / "results" / "milan_alapkod_data_20260930_105152.csv"
OUT = FOLDER / "results"
BASELINE_WINDOW = (4000.0, 4700.0)
SEARCH_START = 4700.0
ONSET_FRACTION, END_FRACTION = 0.75, 0.10
MARKER_RMS_WINDOW = 100.0
ENVELOPE_CUTOFF = 0.05  # cycles per simulation-time unit


def moving_rms(signal, width):
    samples = max(3, int(round(width)))
    if samples % 2 == 0:
        samples += 1
    return np.sqrt(np.convolve(signal**2, np.ones(samples)/samples, mode="same"))


def fft_lowpass(signal, dt, cutoff):
    """Brick-wall FFT low-pass, with reflection padding to limit edge wrap."""
    pad = len(signal)//2
    extended = np.pad(signal, (pad, pad), mode="reflect")
    frequencies = np.fft.rfftfreq(len(extended), dt)
    spectrum = np.fft.rfft(extended)
    spectrum[frequencies > cutoff] = 0.0
    filtered = np.fft.irfft(spectrum, n=len(extended))
    return filtered[pad:pad+len(signal)]


def decay_markers(time, signal, dt):
    smooth = moving_rms(signal, MARKER_RMS_WINDOW/dt)
    baseline = (time >= BASELINE_WINDOW[0]) & (time < BASELINE_WINDOW[1])
    baseline_level = float(np.median(smooth[baseline]))
    starts = np.flatnonzero((time >= SEARCH_START)
                            & (smooth <= ONSET_FRACTION*baseline_level))
    if not len(starts):
        raise RuntimeError("Could not find the start of the central-field decay.")
    start = int(starts[0])
    ends = np.flatnonzero((np.arange(len(time)) >= start)
                          & (smooth <= END_FRACTION*baseline_level))
    if not len(ends):
        raise RuntimeError("Could not find the end of the selected decay interval.")
    return smooth, baseline_level, start, int(ends[0])


def dominant_frequency(time, signal, start, end, dt):
    """Find the strongest nonzero Fourier peak around the marked collapse."""
    use = (time >= start-400.0) & (time <= end+150.0)
    values = signal[use]
    values = (values-values.mean())*np.hanning(len(values))
    frequencies = np.fft.rfftfreq(len(values), dt)
    power = np.abs(np.fft.rfft(values))**2
    candidates = (frequencies >= 0.05) & (frequencies <= 0.5)
    peak = int(np.flatnonzero(candidates)[np.argmax(power[candidates])])
    return float(frequencies[peak]), frequencies, power


def exponential_fit(time, envelope, start, end):
    """Fit A*exp(-lambda*(t-start)) by least squares in envelope units."""
    use = (time >= start) & (time <= end) & np.isfinite(envelope) & (envelope > 0)
    elapsed, values = time[use]-start, envelope[use]
    low, high = 0.001, 0.03
    best = None
    for _ in range(4):
        rates = np.linspace(low, high, 5001)
        basis = np.exp(-rates[:, None]*elapsed[None, :])
        amplitudes = (basis @ values)/np.sum(basis**2, axis=1)
        residuals = values[None, :]-amplitudes[:, None]*basis
        sse = np.sum(residuals**2, axis=1)
        index = int(np.argmin(sse))
        step = (high-low)/(len(rates)-1)
        best = (float(rates[index]), float(amplitudes[index]),
                float(sse[index]), step)
        low = max(0.0, best[0]-4*step)
        high = best[0]+4*step
    rate, amplitude, sse, _ = best
    total = float(np.sum((values-values.mean())**2))
    return {"rate": rate, "amplitude_at_start": amplitude,
            "r_squared_amplitude_space": 1-sse/total if total else float("nan"),
            "sample_count": int(use.sum()), "mask": use}


def log_linear_fit(time, envelope, start, end):
    use = (time >= start) & (time <= end) & np.isfinite(envelope) & (envelope > 0)
    x, y = time[use]-start, np.log(envelope[use])
    slope, intercept = np.polyfit(x, y, 1)
    prediction = intercept+slope*x
    total = np.sum((y-y.mean())**2)
    return {"rate": float(-slope), "amplitude_at_start": float(np.exp(intercept)),
            "r_squared_log_space": (float(1-np.sum((y-prediction)**2)/total)
                                    if total else float("nan"))}


def save_plot(path, time, raw, filtered, envelope, marker_rms,
              frequencies, power, peak_frequency, signal_cutoff,
              start_time, end_time, fit):
    visible = (time >= start_time-400.0) & (time <= end_time+250.0)
    fit_time = time[fit["mask"]]
    fit_curve = fit["amplitude_at_start"]*np.exp(
        -fit["rate"]*(fit_time-start_time)
    )
    fig, axes = plt.subplots(3, 1, figsize=(11, 10),
                             gridspec_kw={"height_ratios": (1, 1.5, 1.3)})

    spectrum_db = 10*np.log10(np.maximum(power/power.max(), 1e-14))
    axes[0].plot(frequencies, spectrum_db, color="#334155", lw=1)
    axes[0].axvline(peak_frequency, color="#059669", ls="--",
                    label=f"dominant peak {peak_frequency:.3f} cycles/time")
    axes[0].axvline(signal_cutoff, color="#dc2626", ls=":",
                    label=f"retained below {signal_cutoff:.3f} cycles/time")
    axes[0].set(xlim=(0, min(0.65, frequencies[-1]),), ylim=(-80, 3),
                ylabel="relative power (dB)",
                title="Fourier spectrum around the marked collapse")
    axes[0].grid(alpha=0.2)
    axes[0].legend(fontsize=8)

    axes[1].plot(time[visible], raw[visible], color="#64748b", lw=0.75,
                 alpha=0.6, label="raw central field")
    axes[1].plot(time[visible], filtered[visible], color="#2563eb", lw=1,
                 label="FFT low-pass central field")
    axes[1].plot(time[visible], envelope[visible], color="#059669", lw=2,
                 label="Fourier-smoothed RMS envelope")
    axes[1].plot(fit_time, fit_curve, color="#dc2626", ls="--", lw=2,
                 label=f"exponential fit: lambda={fit['rate']:.5f}")
    axes[1].axvspan(start_time, end_time, color="#f59e0b", alpha=0.10,
                    label="fit interval")
    axes[1].set(xlabel="time t", ylabel="central field / envelope",
                title="Filtered signal and fitted decay envelope")
    axes[1].grid(alpha=0.2)
    axes[1].legend(fontsize=8, ncol=2)

    fit_envelope = envelope[fit["mask"]]
    axes[2].semilogy(fit_time, fit_envelope, ".", ms=3,
                     color="#059669", label="envelope samples")
    axes[2].semilogy(fit_time, fit_curve, "--", lw=2, color="#dc2626",
                     label="least-squares exponential")
    axes[2].set(xlabel="time t", ylabel="envelope (log scale)",
                title=(f"Amplitude-space R-squared={fit['r_squared_amplitude_space']:.3f}; "
                       f"e-folding time={1/fit['rate']:.1f}"))
    axes[2].grid(alpha=0.2, which="both")
    axes[2].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=SOURCE)
    parser.add_argument("--output-dir", type=Path, default=OUT)
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f"Input CSV not found: {args.input}")

    data = np.genfromtxt(args.input, delimiter=",", names=True)
    time, raw = data["time"], data["phi_center"]
    dt = float(np.median(np.diff(time)))
    nyquist = 0.5/dt
    marker_rms, baseline, i0, i1 = decay_markers(time, raw, dt)
    start_time, end_time = float(time[i0]), float(time[i1])
    peak_frequency, frequencies, power = dominant_frequency(
        time, raw, start_time, end_time, dt,
    )
    signal_cutoff = 1.33*peak_frequency
    filtered = fft_lowpass(raw, dt, signal_cutoff)
    envelope_power = fft_lowpass(filtered**2, dt, ENVELOPE_CUTOFF)
    envelope = np.sqrt(np.maximum(envelope_power, 0.0))
    fit = exponential_fit(time, envelope, start_time, end_time)
    log_fit = log_linear_fit(time, envelope, start_time, end_time)

    # Cutoff dependence and shorter windows check whether one rate describes
    # the full fall or only one phase of it.
    cutoff_sensitivity = []
    for cutoff in sorted(set((0.25, signal_cutoff, 0.35, 0.40))):
        filtered_test = fft_lowpass(raw, dt, cutoff)
        env_test = np.sqrt(np.maximum(
            fft_lowpass(filtered_test**2, dt, ENVELOPE_CUTOFF), 0.0,
        ))
        f = exponential_fit(time, env_test, start_time, end_time)
        cutoff_sensitivity.append({"signal_cutoff": cutoff,
                                   "lambda": f["rate"],
                                   "r_squared": f["r_squared_amplitude_space"]})
    windows = ((start_time, end_time), (start_time, min(end_time, 5100.0)),
               (max(start_time, 5000.0), end_time),
               (max(start_time, 5050.0), end_time))
    window_sensitivity = []
    for a, b in windows:
        f = exponential_fit(time, envelope, a, b)
        window_sensitivity.append({"start": a, "end": b, "lambda": f["rate"],
                                   "e_folding_time": 1/f["rate"],
                                   "r_squared": f["r_squared_amplitude_space"],
                                   "sample_count": f["sample_count"]})

    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = args.output_dir/"fourier_decay_fit.png"
    save_plot(plot_path, time, raw, filtered, envelope, marker_rms,
              frequencies, power, peak_frequency, signal_cutoff,
              start_time, end_time, fit)

    data_path = args.output_dir/"fourier_decay_envelope.csv"
    view = (time >= start_time-400.0) & (time <= end_time+250.0)
    fitted = np.full_like(time, np.nan, dtype=float)
    fitted[fit["mask"]] = fit["amplitude_at_start"]*np.exp(
        -fit["rate"]*(time[fit["mask"]]-start_time)
    )
    with data_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("time", "phi_center_raw", "phi_center_lowpass",
                         "fourier_rms_envelope", "moving_rms_100",
                         "exponential_fit", "in_fit_interval"))
        writer.writerows(zip(time[view], raw[view], filtered[view], envelope[view],
                             marker_rms[view], fitted[view], fit["mask"][view]))

    summary = {
        "source_csv": str(args.input.resolve()),
        "sample_interval": dt, "nyquist_cycles_per_time": nyquist,
        "decay_window": {"start": start_time, "end": end_time,
                         "duration": end_time-start_time,
                         "marker_method": "100-unit moving RMS falls from 75% to 10% of median over 4000<=t<4700"},
        "spectrum": {"dominant_frequency_cycles_per_time": peak_frequency,
                     "dominant_period": 1/peak_frequency,
                     "signal_lowpass_cutoff": signal_cutoff,
                     "envelope_power_lowpass_cutoff": ENVELOPE_CUTOFF},
        "primary_exponential_fit": {
            "model": "E(t)=E_start*exp(-lambda*(t-start))",
            "method": "nonlinear least squares in envelope amplitude",
            "lambda_per_time_unit": fit["rate"],
            "e_folding_time": 1/fit["rate"],
            "E_start": fit["amplitude_at_start"],
            "r_squared_amplitude_space": fit["r_squared_amplitude_space"],
            "sample_count": fit["sample_count"],
        },
        "log_linear_fit_alternative": {
            "lambda_per_time_unit": log_fit["rate"],
            "e_folding_time": 1/log_fit["rate"],
            "E_start": log_fit["amplitude_at_start"],
            "r_squared_log_space": log_fit["r_squared_log_space"],
        },
        "signal_cutoff_sensitivity": cutoff_sensitivity,
        "fit_window_sensitivity": window_sensitivity,
        "plot": str(plot_path.resolve()), "filtered_data": str(data_path.resolve()),
    }
    summary_path = args.output_dir/"fourier_decay_fit_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(f"Decay interval: {start_time:.3f} to {end_time:.3f}")
    print(f"lambda={fit['rate']:.6f} per time unit; e-fold={1/fit['rate']:.2f}; "
          f"amplitude-space R2={fit['r_squared_amplitude_space']:.4f}")
    print(f"Log-linear alternative lambda={log_fit['rate']:.6f}; "
          f"plot: {plot_path.resolve()}")


if __name__ == "__main__":
    main()
