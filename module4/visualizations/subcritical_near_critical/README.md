# Near-critical subcritical animation

Run from the project root:

```bash
python module4/create_subcritical_animation.py
```

The default initial amplitude is the lower edge of the saved Module 4 critical
bracket, the closest tested amplitude that did **not** reach the
`|phi(0,t)| = 12` cutoff within `t = 100`. The animation is centered on the
largest central-field excursion and shows both the radial profile and the
central-field / signed-energy RMS-radius traces. The complete sampled scalar
time series and run summary are saved beside the GIF.

The RMS radius uses the notes' signed energy-weighted `<r^2>` inside the fixed
`R_core = 8`; it is not a positive-density or half-maximum width. The radial
profile is included so the core's shape can be inspected independently of
that moment. A subcritical run is not expected to blow up: it can approach a
focused state and then turn around or disperse.

For this run, `|phi(0,t)|` peaks at about `5.871` near `t=20.10`. In the
animated window (`t=14.20` to `32.00`) the RMS radius decreases from about
`2.93` to `2.05`, then grows to about `4.74`; the subcritical core contracts
transiently and later re-expands.

The simulation retains the `|phi(0,t)| = 12` safety cutoff. Continuing a
supercritical solution past that value should be treated as a separate
resolution-convergence experiment, not simply as waiting longer: a narrowing
core can become under-resolved on the fixed radial grid.
