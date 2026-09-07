# Parameter sweep

The first sweep varies only the two physical initial-condition parameters:

| parameter | values used | reason |
|---|---|---|
| `A` | `0.4, 0.8, 1.2, 1.6, 2.0` | Dimensionless amplitudes spanning the amplitude range shown for the 3D phi^4 model by Gleiser and Sicilia. |
| `R0` | `2.0, 2.5, 2.86, 3.5, 5.0` | The `2 <= r0 <= 5` interval is explicitly scanned by Honda and Choptuik; `2.86` is their reported long-lived example in the Gaussian description used by Gleiser and Sicilia. |

The sweep therefore contains 25 runs. `m`, `lambda`, `alpha`, `R_max`, `N`,
`CFL`, `dt`, `T`, and `output_every` remain fixed and are written into every
`run_info.txt`. The mass and coupling are kept at their dimensionless reference
values (`m = 1`, `lambda = 1`); changing them independently would change the
normalization rather than probe the same dimensionless model.

These are literature-informed starting values, not universal existence bounds.
The papers use slightly different vacuum conventions and potentials, so the
sweep is also testing how those initial-condition ranges behave in this code.
The current `T = 10` is an exploratory time window; it is much shorter than
the approximately `7100` lifetime reported for the long-lived example.

References:

1. E. P. Honda and M. W. Choptuik, *Fine Structure of Oscillons in the
   Spherically Symmetric phi^4 Klein-Gordon Model*, Phys. Rev. D 65, 084037
   (2002), [arXiv:hep-ph/0110065](https://arxiv.org/abs/hep-ph/0110065).
2. M. Gleiser and D. Sicilia, *A General Theory of Oscillon Dynamics*, Phys.
   Rev. D 80, 125037 (2009), [arXiv:0910.5922](https://arxiv.org/abs/0910.5922).

## Output layout

```text
output/
├── snapshots.dat             # ordinary single run
├── run_info.txt
└── sweep/
    ├── manifest.csv          # index of all sweep points
    ├── run_001_A0p40_R02p00/
    │   ├── snapshots.dat
    │   ├── run_info.txt
    │   └── video/             # created when this run is animated
    └── ...
```

To create the sweep:

```bash
make sweep
```

The current model is not identical to the literature models, so a high
amplitude case can collapse numerically. Such a case is kept in the manifest
with `status=nonfinite`; the other cases still finish and can be animated.
The first test sweep produced 9 `ok` runs and 16 `nonfinite` runs.

To animate one result later:

```bash
python3 animate.py --run output/sweep/run_001_A0p40_R02p00
```
