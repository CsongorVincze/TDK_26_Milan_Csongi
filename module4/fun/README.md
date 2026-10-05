# Critical-threshold animation

This folder contains a 3×3 animation of four subcritical amplitudes, the saved
N=500 critical-amplitude estimate, and four supercritical amplitudes. The four
amplitudes increase from left to right, then top to bottom. The critical
estimate is centered, with the closest subcritical and supercritical cases on
its left and right.

The panels share one field and radius scale. Each supercritical run freezes
when the central field first reaches `|phi(0,t)| = 12`; the estimated critical
run is deliberately shown too. This is a finite numerical cutoff, not a claim
that the animation proves a mathematical singularity. The displayed `Acrit`
is the operational estimate from the N=500 threshold search.

For these chosen offsets, the four subcritical runs are followed through
`t = 100`; the estimated critical run reaches the cutoff at about `t = 20.69`,
and the four supercritical runs reach it earlier, from about `t = 15.15` down
to `t = 4.47` as the amplitude excess increases. None of the subcritical runs
reaches the cutoff by `t = 100`; their central fields are close to zero by then,
so these selected cases disperse rather than persist as long-lived oscillons.
This finite run cannot prove they will never collapse at still later times.

Open [the animation](results/critical_transition_3x3.gif). To regenerate it
from the PDE solver:

```bash
python3 module4/fun/critical_transition_animation.py
```

The exact amplitudes and cutoff times are recorded in
[case_summary.csv](results/case_summary.csv).
