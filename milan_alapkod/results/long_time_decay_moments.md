# Baseline decay: central field and `<r^2>`

This uses the existing `milan_alapkod_data_20260930_105152.csv`. No simulation
rerun was needed, and the original CSV was left unchanged.

“Central field” means `phi_center = phi(t, r=0)`. The plot puts this signed,
oscillating field and its 100-time-unit moving RMS on the left axis, and the
CSV's energy-weighted `<r^2>(t)` on the right axis. The axes have separate
scales.

The stored moment is evaluated inside `R_core = 8`:

```text
<r^2>(t) = integral_0^Rcore r^2 E(r,t) 4*pi*r^2 dr
           / integral_0^Rcore E(r,t) 4*pi*r^2 dr
core_rms_radius = sqrt(<r^2>)
```

The decay markers use the median central-field RMS over `4000 <= t < 4700` as
the pre-decay reference. The marked start is when the smoothed RMS first falls
below 75% of that reference (`t = 4964.196`); the 10% marker is an operational
decay endpoint (`t = 5145.515`), not a finite-time blow-up time.

Around `t = 5000`, the existing data give `phi_center = 0.6384` and
`<r^2> = 17.0746`. The moment is `14.6789` at the marked decay start and
`35.2758` at the 10% endpoint; it later peaks near `49.2709` at `t = 5216.9`.
So the field at the origin quiets while the energy-weighted core spreads
outward.

No beta fit is made yet. For the previously discussed law,

```text
log(sqrt(<r^2>)) = beta*log(t_c - t) + constant
```

fitting `log(<r^2>)` directly would give slope `2*beta`, so that slope must be
halved when comparing with the beta in the equation above. The plotted
`t = 5145.515` marker is only an observable-based decay endpoint; choosing a
physically meaningful `t_c` and fit interval is a separate decision.

Files: [plot](long_time_decay_moments.png),
[plotted data window](long_time_decay_moments_window.csv), and
[numerical marker summary](long_time_decay_moments_summary.json). Recreate the
plot with `python milan_alapkod/plot_long_time_decay.py`.
