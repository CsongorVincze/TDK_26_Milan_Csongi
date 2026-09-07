# Claim-source ledger (internal research record)

Retrieved and checked: 6 September 2026. This ledger records the consequential claims used in `report-source.md`. It is an internal audit artifact, not part of the delivered report.

| Claim or decision | Evidence used | Evidence boundary / uncertainty | Confidence |
|---|---|---|---|
| An oscillon is a finite-energy, localized, quasiperiodic real-field configuration that slowly radiates. | Fodor 2019 review, arXiv:1911.03340; Copeland, Gleiser, Mueller 1995, arXiv:hep-ph/9503217. | Terminology varies across communities; report explicitly distinguishes related objects. | High |
| A fundamental frequency below the mass threshold is localized, while harmonics above it can propagate. | Fodor et al. 2009, arXiv:0812.1919; Nagy and Takacs 2021, arXiv:2105.01089. | Which harmonics are present depends on symmetry and model. | High |
| Small-amplitude classical radiation is beyond all algebraic orders and exponentially weak. | Fodor et al. 2009, arXiv:0812.1919; Fodor 2019 review. | Controlled in the small-amplitude regime; not a universal lifetime formula. | High |
| Geometric decoupling and destructive interference can suppress radiation. | Cyncynates and Giurgica-Tiron 2021, arXiv:2104.02069, Sec. III. | Radiation minima are model- and frequency-dependent. | High |
| The U(1) symmetry and particle number are approximate in a real-scalar oscillon EFT. | Mukaida, Takimoto, Yamada 2017, arXiv:1612.07750; Levkov et al. 2022, arXiv:2208.04334. | Exact within finite-order averaged/gradient EFTs, broken by nonperturbative relativistic emission in the parent real theory. | High |
| Quantum emission or coupling to daughter fields can dominate a very small classical rate. | Hertzberg 2010, arXiv:1003.3459. | Perturbative/model-dependent; report does not turn it into a universal physical lifetime. | Medium-high |
| Post-inflation self-resonance can form oscillons in potentials shallower than quadratic away from a quadratic minimum. | Amin et al. 2012, arXiv:1106.3335. | Requires specified potential and sufficiently weak coupling to other fields. | High |
| Oscillons can form during subcritical bubble collapse. | Copeland, Gleiser, Mueller 1995, arXiv:hep-ph/9503217. | A demonstrated formation channel, not a universal abundance result. | High |
| Gauge-Higgs oscillons exist in a special numerical parameter regime. | Graham 2007, arXiv:hep-th/0610267. | Proof of principle; not the measured Standard Model mass ratio. | High |
| A mature exactly spherical oscillon does not itself emit gravitational waves. | Zhou et al. 2013, arXiv:1304.6094. | Formation, asymmetry, collisions, decay, or dynamical gravity can source waves. | High |
| A reported late high-frequency GW feature can be a numerical artifact. | Amin et al. 2018, arXiv:1803.08047. | Demonstrated in one hilltop setup; used as a numerical caution, not a universal null result. | High |
| Granular-media 'oscillons' are driven-dissipative objects distinct from scalar-field oscillons. | Umbanhowar, Melo, Swinney 1996, Nature 382, 793. | Shared name and phenomenological localization only. | High |
| In standard 3D quasibreather families, terminal 'collapse' often means rapid dissolution near an E(omega) minimum. | Saffin and Tranberg 2007, arXiv:hep-th/0610191; Nagy and Takacs 2021, arXiv:2105.01089. | Model-class and spherical-symmetry evidence; not a theorem for every potential. | High |
| Staccato bursts arise at successive n omega = m thresholds as omega increases. | Dorey, Romanczukiewicz, Shnir 2020, arXiv:1910.04128; Nagy and Takacs 2021. | Parity can remove even harmonics; strong staccato behavior is not universal across dimension/potential. | High |
| Cubic focusing NLS is mass-critical in d=2. | Standard NLS scaling; MIT 18.156 notes; Holmer and Roudenko 2007 for 3D blow-up context. | This is an envelope-PDE classification, not a direct oscillon taxonomy. | High |
| Small-amplitude localized quasibreather cores have Dcrit=4, while E ~ epsilon^(2-D) changes limiting behavior at D=2. | Fodor et al. 2008, arXiv:0802.3525. | Small-amplitude analysis; the paper explicitly does not prove absolute absence of oscillons at D>=4. | High |
| Levkov et al. obtain dN/domega = omega^-1 dE/domega < 0 as a linear-stability condition in their EFT. | Levkov et al. 2022, Appendix D, arXiv:2208.04334. | Long-wavelength EFT criterion; not a proof of singular collapse in full KG. | High |
| For cubic NLS scale balance, amplitude scales as 1/L, so L ~ sqrt(T-t) gives exponent -1/2 rather than -3/4. | NLS scaling; Akrivis et al. 2003, DOI 10.1137/S1064827597332041; Holmer and Roudenko 2007, arXiv:math/0703233. | The square-root rate describes the standard point-collapse branch, not every possible blow-up solution. | High |
| Parent-wave corrections can arrest an NLS collapse. | Fibich, Ilan, Tsynkov 2003, SIAM J. Appl. Math. 63, 1718. | Optical nonlinear-Helmholtz analogy, not direct evidence for KG oscillons; report labels it as such. | Medium |
| The local code uses V=phi^2/2-phi^4/4, with a barrier at |phi|=1 and no lower bound. | Direct inspection of `AI SLOP/cpp/potential.hpp` and `Csongi_gyakorlos_sim/solver.cpp`; elementary differentiation. | None for the algebra; physical interpretation of a nonfinite trajectory still requires resolution tests. | Very high |
| Fifteen of 25 coarse-sweep initial centers lie beyond the barrier; the manifest has 9 finite and 16 nonfinite runs at T=10. | `Csongi_gyakorlos_sim/parameter_sweep.md`; `Csongi_gyakorlos_sim/output/sweep/manifest.csv`; local run metadata. | 'Finite' is not an oscillon classification; T=10 is very short. | Very high |
| The sweep boundary and output design cannot establish a long-time oscillon lifetime. | Direct inspection of `Csongi_gyakorlos_sim/solver.cpp` and parameter notes. | Report states concrete limitations rather than inferring a corrected result. | High |
| The radial origin formula and method-of-lines architecture are useful foundations. | Direct code inspection; standard radial Laplacian limit. | The full implementation remains unverified until convergence, energy, and boundary tests are passed. | High |
| The exact local citation 'Salmi & Sutcliffe, J. Phys. A 45, 465201' is false. | Journal coordinates/DOI map to Akemann and Burda 2012, arXiv:1208.0187; likely intended oscillon paper is Salmi and Hindmarsh 2012, arXiv:1201.1934. | 'Likely intended' is an inference from topic and year; false-coordinate finding is exact. | Very high |
| The bounded double well V=(phi^2-1)^2/4 is an appropriate canonical first benchmark. | Fodor et al. 2006, arXiv:hep-th/0609023; Honda and Choptuik 2002, arXiv:hep-ph/0110065. | Other bounded benchmarks are valid; this one is chosen for comparability and clarity. | High |
| Reliable radial evolutions require energy/flux accounting, resolution convergence, and boundary-domain tests. | Honda and Choptuik 2002; Fodor et al. 2006; Langtangen and Linge 2017. | Exact discretization and tolerances remain a project decision. | High |

## Search and selection trail

- Started from the local PDFs' named concepts and citations: oscillon, quasibreather, staccato radiation, NLS criticality, Vakhitov-Kolokolov, dimensional collapse, radial KG numerics.
- Prioritized primary papers, author-supplied media, institutional course material, and Fodor's specialist review.
- Cross-checked publication metadata and the suspicious J. Phys. A coordinates independently rather than repairing the citation by title alone.
- Used popular/institutional material only for orientation and video discovery, not for consequential technical claims.
- Recorded direct links and publication dates in the report; claims that depend on a potential, dimension, symmetry, or approximation are labeled accordingly.

