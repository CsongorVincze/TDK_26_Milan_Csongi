# Broad amplitude search: resolved windows

**6 finite-horizon non-crossing windows and 7 collapse bands were sampled. 6 window interiors pass the finer-grid and larger-domain checks; 7 collapse-band interiors agree across N=999,1997,3993.**

The search covers A=0.2 to 4.81044, including subcritical starts and raw positive offsets 1e-08 to 3. There are 3783 primary configurations, with adaptive sampling around state changes and two extra passes through non-crossing gaps.

A non-crossing window (G labels) is a sampled interval where field 48 is not reached by t=100, bounded by collapse samples on both sides. A collapse band (B labels) has a fixed final field sign and number of central-field zero crossings. Both are ordered by amplitude. A sign change or oscillation-count change is reported as a branch change, separately from a non-crossing window. This finite scan is not a proof of the exact number of continuum windows; smaller unresolved windows may exist.

The reference is Ac=1.2026097332023025 with estimated uncertainty 5e-12. All offset values are raw `(A-Ac)/Ac`, never amplitude percentages. The primary method is fourth order at N=999, CFL=0.4, Rmax=25, sponge starts at 20, R0=2, zero initial velocity. Collapse classification uses field 48; the plot retains t12 to match the original cutoff-time figure. Neither t12 nor t48 is the singularity time tc.

Among primary non-crossing runs, 0 reached field 12 without reaching 48. A zero count means the observed classification also agrees with the earlier field-12 criterion.

## Collapse bands

| Band | Sampled raw delta interval | Final sign | Zero crossings | Finer-grid interior confirmed |
| --- | --- | ---: | ---: | --- |
| B01 | [9.99999991e-09, 2.98887551e-06] | +1 | 6 | 1 |
| B02 | [3.5732098e-06, 3.47651335e-05] | -1 | 5 | 1 |
| B03 | [3.73386303e-05, 0.000286625791] | +1 | 4 | 1 |
| B04 | [0.000300727615, 0.00236158844] | -1 | 3 | 1 |
| B05 | [0.00244027042, 0.0266417] | +1 | 2 | 1 |
| B06 | [0.0268527291, 0.274792603] | -1 | 1 | 1 |
| B07 | [0.27486349, 3] | +1 | 0 | 1 |

## Non-crossing windows

| Window | Sampled raw delta interval | Validation confirmed |
| --- | --- | --- |
| G01 | [2.98888115e-06, 3.57320306e-06] | 1 |
| G02 | [3.47651611e-05, 3.73386007e-05] | 1 |
| G03 | [0.000286625962, 0.000300727436] | 1 |
| G04 | [0.00236158985, 0.00244026897] | 1 |
| G05 | [0.0266417158, 0.0268527131] | 1 |
| G06 | [0.274792767, 0.274863326] | 1 |

The endpoints above are the first/last sampled amplitudes of each group. `windows.csv` also lists the neighboring outside samples that bracket each boundary. They are N=999 sampled boundaries, not extrapolated continuum edge locations. Resolution checks validate representative interiors, not every boundary or the exact width of every window.


The extra near-critical extension contains 81 runs, with observed states ['+1/6']. State +1/6 means positive-field collapse after six central-field zero crossings. No additional window was detected in this extension if only that state is present.

## Validation and limitations

Each collapse-window representative is rerun at the same amplitude on N=1997 and N=3993. Each gap representative is checked on N=1997 and in a larger domain (radius at least 75, sponge starts at 60) through t=200. `resolution_checks.csv` records exact agreement or disagreement; a discrepancy is retained rather than hidden. Extra near-critical runs at N=1997 cover raw delta=1e-10..1e-8 separately; they do not inflate the main window count. A no-crossing result means only no field-48 event by the observation horizon.

The tested range is broader than the old endpoint 1.2 Ac: it extends to 4 Ac. Sampling is logarithmic near criticality, supplemented by a linear broad scan. Adaptive subdivision probes changes in sign/oscillation count and fills known gaps; a finite search cannot rule out arbitrarily narrow islands between equal-outcome samples.

## Files

- `window_map.png`: main cutoff-time/log(raw offset) figure and oscillation counts.
- `amplitude_overview.png`: the full linear amplitude range.
- `gap_details.png`: individual non-crossing intervals.
- `representative_trajectories.png`: finer-grid examples.
- `amplitude_samples.csv`: all primary amplitudes and outcomes.
- `windows.csv`: sampled window extents and exterior boundary brackets.
- `resolution_checks.csv`: finer-grid and larger-domain checks.
- `../data/runs.jsonl`: complete cached computations, including validation runs.
- `../data/representative_trajectories.npz`: compressed representative data; columns time, phi_center, pi_center, r50.
