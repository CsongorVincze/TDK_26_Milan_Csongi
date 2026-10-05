# Fourier-filtered central-field decay

This analysis uses the existing `milan_alapkod_data_20260930_105152.csv`; the simulation was not rerun and the source CSV was not changed. It fits the oscillation-amplitude envelope of the central field `phi_center = phi(t, r=0)`, not `<r^2>`.

## Exponential parameter

On the operational collapse/decay interval `4964.196 <= t <= 5145.515`, the amplitude-space least-squares fit is

```text
E(t) = 0.5666 * exp(-0.010562 * (t - 4964.196))
```

The effective decay rate is therefore `lambda = 0.01056` per simulation-time unit, corresponding to an e-folding time `1/lambda = 94.7` time units. The amplitude-space R-squared is 0.942 over 260 samples.

This is an effective rate over the full marked decay, not a constant rate that describes every part equally well. Fitting later portions gives a faster decay, so use the fit window whenever quoting lambda.

## Fourier and envelope method

The CSV is uniformly sampled at `dt = 0.70007`; its Nyquist frequency is 0.714 cycles per time unit. In the spectrum around the collapse, the strongest oscillation peak is at 0.217 cycles per time unit (period about 4.60). The central-field trace is Fourier low-pass filtered at 0.289 cycles per time unit, retaining that main oscillation while removing higher-frequency content. To turn the oscillations into a positive envelope, the filtered field is squared, its power is low-pass filtered at 0.05 cycles per time unit, and the square root is taken. Reflection padding is used around the finite trace before FFT filtering to limit edge wrap-around.

The fit interval is selected using the same operational markers as the existing decay plot: the 100-time-unit moving RMS first falls below 75% of its pre-decay baseline at `t=4964.196`, and later falls below 10% at `t=5145.515`. These mark an observed decay interval; the latter is not assumed to be a singularity time.

The filter choice is not driving the rate: with signal cutoffs from 0.25 to 0.40 cycles per time unit, lambda stays between 0.010562 and 0.010563. The fit interval matters much more:

| Fit interval | Lambda | E-folding time | Amplitude-space R-squared |
| --- | ---: | ---: | ---: |
| 4964.2–5145.5 (full marked decay) | 0.01056 | 94.7 | 0.942 |
| 4964.2–5100 | 0.00943 | 106.0 | 0.926 |
| 5000–5145.5 | 0.01421 | 70.4 | 0.984 |
| 5050–5145.5 | 0.01738 | 57.5 | 0.998 |

That trend shows the envelope decays faster later in the episode. If one wants a single number for the complete 75%-to-10% decay, use `lambda ~= 0.0106`; for the later collapse tail, the effective rate is closer to `0.014–0.017`.

A second common fitting convention regresses `log(E)` linearly against time. It gives `lambda=0.01353` (log-space R-squared 0.966; e-folding time 73.9). This differs because log-space fitting weights relative errors, while the reported primary fit minimizes errors in the envelope amplitude itself. Neither fitting convention removes the interval dependence.

## Outputs and reproduction

- [Filtered signal, spectrum, envelope, and fit](fourier_decay_fit.png)
- [Filtered envelope samples and fitted curve](fourier_decay_envelope.csv)
- [Numerical parameters and sensitivity checks](fourier_decay_fit_summary.json)
- [Analysis script](../fourier_decay_fit.py)

Run from the project root with:

```bash
MPLCONFIGDIR=/tmp/mplconfig_alapkod python3 milan_alapkod/fourier_decay_fit.py
```
