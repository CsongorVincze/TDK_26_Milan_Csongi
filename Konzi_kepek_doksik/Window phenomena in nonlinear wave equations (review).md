# Window phenomena in nonlinear wave equations: mechanism and literature

Oct 23, 2026 · Gábor Takács

In many nonlinear wave equations the outcome of an evolution does not switch only once as a parameter of the initial data is varied. Narrow intervals of the opposite outcome, called windows, interrupt one regime: kinks that should be captured escape, pulses that should collapse disperse. This note explains the mechanism common to these phenomena, surveys the main examples in the literature and gives a reading list.

Section 1 defines a window and the one mechanism behind all the examples. Sections 2 and 3 cover kink collisions and solitons meeting defects, Section 4 the two ways in which windows can accumulate, including critical gravitational collapse, and Section 5 resonances in oscillon lifetimes. Section 6 collects what is known about the collapse threshold of the focusing Klein-Gordon equation, where a theorem rules windows out near the ground-state energy. Section 7 compares the phenomena in one table, and Section 8 has questions to think about. Appendix A covers collapse in anti-de Sitter space, which needs some general relativity; skip it on a first reading.

## 1. What is a window?

A window is a narrow interval of a control parameter p in which the outcome is the opposite of the outcome on both sides of it. Take a one-parameter family of initial data (an amplitude, a velocity, a radius) and a two-valued outcome (collapse or dispersal, escape or capture). In simple problems one critical value p\_c separates the two outcomes. Windows appear when the solution can come back to a decision point more than once.

All the examples in this note share one mechanism:

- The solution returns to a decision point several times: two kinks collide again, a core compresses again, a pulse comes back from a wall.
- At each return the outcome can be settled. Whether it is depends on the state of a second degree of freedom at that moment: an oscillating internal mode, or a growing unstable mode.
- A window is a parameter interval in which the decision falls the other way at the n-th return.

How the windows crowd towards their accumulation point is set by how the return time T(p) depends on p. Suppose window n requires T(p\_n) ≈ n T₀ + const. If T diverges as a power, as near a separatrix (a ball thrown just below escape velocity), then

$$
T(p) \propto (p_c - p)^{-a} \quad\Longrightarrow\quad p_c - p_n \propto n^{-1/a}.
$$

If T diverges logarithmically, as near an unstable state whose perturbations grow like e^{λt}, then

$$
T(p) \simeq \frac{1}{\lambda}\ln\frac{1}{|p-p_*|} \quad\Longrightarrow\quad |p_n-p_*| \propto e^{-\lambda T_0 n} = \rho^{-n}.
$$

The first case gives windows that pile up as a power law (kink collisions, Section 2). The second gives a geometric sequence with ratio ρ = e^{λT₀} (Choptuik's fine structure, Section 4). The second number that characterises a window sequence is its length: infinite if the returns can continue for ever, finite if the state that produces them does not last.

## 2. Bounce windows in kink–antikink collisions

In the 1D double-well (φ⁴) model, a kink and an antikink colliding below a critical speed v\_c ≈ 0.26 usually capture each other. In a sequence of escape windows they separate after two bounces instead, and the windows pile up at v\_c as a power law: v\_c² − v\_n² ∝ n⁻². This is the classic window phenomenon (Campbell, Schonfeld & Wingate 1983; see the review by [Kevrekidis & Goodman](https://arxiv.org/pdf/1909.03128)).

**The model.** With the potential V = (φ² − 1)²/4 in one dimension,

$$
\phi_{tt} - \phi_{xx} + \phi^3 - \phi = 0, \qquad \phi_K(x) = \tanh\frac{x}{\sqrt2} .
$$

Small oscillations about the kink have one bound mode, the shape mode, with frequency ω\_s = √(3/2) ≈ 1.22. It lies below the meson mass √2, where the continuum starts. Speeds are unchanged by rescaling x and t together, so v\_c is the same in every normalisation of the model.

**What happens.**

- Above v\_c ≈ 0.2598 the kinks bounce once and separate.
- Below v\_c they are mostly captured into a long-lived oscillating lump, a bion, which is the 1D cousin of an oscillon.
- In narrow windows they escape after two bounces. Three-bounce windows sit at the edges of the two-bounce windows, and windows with more bounces at the edges of those. The result is a nested, fractal structure (Anninos, Oliveira & Matzner 1991).

**The mechanism.**

1. The first collision moves part of the kinetic energy into the shape modes of the two kinks. They fly apart, but no longer have enough energy to escape, and return after a time T₂.
2. If the shape modes are in the right phase at the second collision, the energy flows back into motion and the pair escapes. This requires ω\_s T₂ ≈ 2πn + δ, with δ a fixed phase offset.
3. The time apart grows as the energy deficit shrinks, T₂ ∝ (v\_c² − v²)^{−½}, like a ball thrown just below escape speed. Together with step 2 this gives

$$
v_c^2 - v_n^2 \propto \frac{1}{(2\pi n + \delta)^2}.
$$

**Further reading.**

- [Goodman & Haberman 2005](https://epubs.siam.org/doi/10.1137/050632981) and [2007](https://arxiv.org/abs/nlin/0702048) reduce the collision to an iterated separatrix map from one bounce to the next. It predicts the window positions and the fractal nesting. See also [Goodman 2008](https://www.arxiv.org/abs/0710.3209) and, for nonlinear Schrödinger equations, [Zhu, Haberman & Yang](https://arxiv.org/pdf/0802.3595).
- In the φ⁶ model the kinks have no shape mode, yet windows appear ([Dorey, Mersh, Romanczukiewicz & Shnir 2011](https://arxiv.org/pdf/1101.5951)). The energy is parked in a bound state trapped between the kink and the antikink.

## 3. Windows when a soliton meets a defect

A single kink or soliton hitting a localised defect shows the same two-bounce windows. A mode bound to the defect plays the role of the shape mode. The set-up has one soliton and a fixed defect, so it is simpler than a collision and allows explicit formulas.

- **Sine-Gordon kink and a point defect** ([Goodman & Haberman 2004](https://arxiv.org/pdf/nlin/0310040)). Below a critical velocity the kink is trapped by the defect. The exceptions are windows in which it comes free after two bounces.
  - Melnikov theory gives the critical velocity. Matched asymptotics give the centre of the n-th window, which again obeys V\_c² − V\_n² ∝ 1/n².
  - For large n the window widths shrink as n⁻³.
  - Three-bounce windows hug the edges of the two-bounce windows.
- **NLS soliton and a defect** ([Goodman, Holmes & Weinstein 2004](https://www.arxiv.org/abs/nlin/0203057)). Depending on the parameters the soliton is reflected, transmitted or captured. The mechanism is resonant energy transfer to a nonlinear standing wave held by the defect. A collective-coordinate model reduces the problem to a Hamiltonian system with three degrees of freedom. The stable and unstable manifolds of its periodic orbits organise the outcomes.

## 4. Geometric versus power-law spacing

The way windows pile up tells you what the solution lingers near. Lingering near an unstable state gives geometric spacing. Lingering near a separatrix, as in a ball thrown just below escape velocity, gives power-law spacing.

**Near an unstable state.** Suppose the family of data crosses, at p = p\*, the stable manifold of a state with one unstable mode of growth rate λ. The deviation then grows like |p − p\*| e^{λt}, and the time to leave the state is

$$
t(p) \simeq \frac{1}{\lambda}\ln\frac{1}{|p-p_*|} + \text{const}.
$$

If the outcome can only be decided at discrete returns spaced T₀ apart, the logarithm turns into a geometric sequence of windows with ratio e^{λT₀}, as in Section 1.

**Near a separatrix.** Below v\_c the kink and antikink fly apart with a small energy deficit proportional to v\_c² − v². The time until they fall back grows like (v\_c² − v²)^{−½}, a power law. With the resonance condition this gives the n⁻² law of Section 2.

**Choptuik's critical collapse.** The best-known geometric structure at a threshold comes from gravitational collapse of a massless scalar field (Choptuik 1993, reviewed by [Gundlach & Martín-García 2007](https://arxiv.org/abs/0711.4620)). You need no general relativity to follow the logic.

- Fine-tuned data approach a critical solution that repeats itself on ever smaller scales, with period Δ = 3.44 in the logarithm of the scale. This critical solution has one unstable mode.
- Just above threshold the black-hole mass scales as M ∝ (p − p\*)^γ with γ ≈ 0.374.
- [Hod & Piran 1997](https://arxiv.org/abs/gr-qc/9606087) found a periodic wiggle on top of this power law. It is periodic in ln(p − p\*) with measured period ≈ 4.6.
- 4.6 = Δ/(2γ). The factor 2 arises because the scalar field changes sign after half an echo, while the geometry has period Δ/2 ([Gundlach's review](https://arxiv.org/pdf/gr-qc/0210101), Fig. 3). The same halving applies to any field whose critical solution changes sign after half a period.

**Type I collapse with oscillating critical solutions.** In gravitational collapse with a mass scale, the critical solution can be a star-like object that persists for a while and then either collapses or disperses. Such a critical solution may be periodic in time (a limit cycle), and the time spent near it follows the logarithmic law above ([Gundlach & Martín-García](https://arxiv.org/abs/0711.4620)).

- For a massive scalar field the critical solutions are oscillating soliton stars ([Brady, Chambers & Gonçalves 1997](https://arxiv.org/pdf/gr-qc/9709014)).
- For perturbed boson stars they resemble unstable boson stars. These persist for a finite time before either collapsing or dispersing ([Hawley & Choptuik 2000](https://arxiv.org/abs/gr-qc/0007039)).

## 5. Resonances in oscillon lifetimes

In the 3D double-well (φ⁴) model, the lifetime of an oscillon formed from a collapsing bubble depends on the initial radius r₀ through a forest of sharp resonances. Near each resonance the lifetime grows logarithmically, the single-unstable-mode law of Section 4.

- [Copeland, Gleiser & Müller 1995](https://arxiv.org/abs/hep-ph/9503217) found that bubbles above a minimum radius settle into oscillons. They also found that the lifetime depends on r₀ and explained why some oscillons live longer than others.
- [Honda & Choptuik 2002](https://ar5iv.labs.arxiv.org/html/hep-ph/0110065) resolved 125 resonances in r₀, each to one part in 10¹⁴. Near a resonance r₀\* the lifetime follows T ∼ γ ln|r₀ − r₀\*|. They interpret each resonance as a critical solution with one unstable mode, with γ tied to the inverse of its growth rate, by analogy with Type I critical collapse (Section 4).
- [Gleiser & Krackow 2019](https://arxiv.org/abs/1906.04070) studied lifetime spikes for particular Gaussian data in 3D. Their numerics are consistent with the conjecture that resonant oscillons may live for ever.

A lifetime resonance is a tuned point, not an interval, so it is the analogue of a threshold (a window edge) rather than of a window.

## 6. Collapse thresholds of the focusing Klein-Gordon equation

For the focusing cubic Klein-Gordon equation

$$
\phi_{tt} - \Delta\phi + \phi - \phi^3 = 0
$$

data either scatter (disperse) or blow up in finite time. The threshold between the two is organised by the ground state Q, the positive radial static solution of −ΔQ + Q − Q³ = 0, which has a single unstable mode. Windows of the collapse-or-dispersal kind need the solution to pass close to ±Q more than once.

**The one-pass theorem.** [Nakanishi & Schlag](https://www.arxiv.org/pdf/1005.4894) treat the radial equation in 3D with energy below E(Q) + ε², with ε small. A solution can enter and leave a small neighbourhood of ±Q at most once.

- Forwards and backwards in time each solution then scatters, blows up or stays trapped near ±Q, giving nine possible combinations.
- In this energy range the threshold is a single smooth surface (the centre-stable manifold of Q). Every solution makes at most one pass near ±Q, so the repeated-return mechanism of Section 1 cannot operate.
- [Krieger, Nakanishi & Schlag](https://arxiv.org/pdf/1010.3799) prove the analogous result for the massless energy-critical wave equation, radial, in 3 and 5 dimensions, slightly above the energy of the Aubin–Talenti solution W.

The theorem needs ε to be small, so it does not cover data well above E(Q), for example at twice the ground-state energy.

**Numerics at higher energy.** [Donninger & Schlag 2011](https://arxiv.org/abs/1011.2015) computed the 3D equation on two-parameter slices of initial data, for example (A f, B g) with Gaussian profiles f and g.

- Near E(Q) the boundary between blowup and global existence is a union of smooth curves.
- At higher energy, thin spikes and bands of the scattering region reach into the blowup region, and more appear as one zooms in. They state that it is completely unclear whether infinitely many emerge.
- Runs near the boundary oscillate about ±Q before deciding.

A one-parameter family of data is a line through such a slice. Where the line cuts one of the spikes, it shows a window.

**One dimension.** [Bizoń, Chmaj & Szpak](https://ar5iv.arxiv.org/html/1012.1033) studied the threshold of the 1D focusing Klein-Gordon equation with a power nonlinearity.

- In the cubic case (their α = 1), u\_tt − u\_xx + u − u³ = 0, the static solution is S = √2 sech x, with one unstable mode of growth rate √3.
- Data tuned to the threshold are trapped by S, then leave along its single unstable direction.
- They study single thresholds, with no repeated passes.
- They find slow convergence to S for α = 1, caused by a zero-energy resonance, and link very slow convergence for α < 1 to oscillatory modes.

**The energy-critical case.** In four dimensions the cubic nonlinearity is energy-critical, and the massive equation has no static solution of finite energy.

- The threshold theory for the critical Klein-Gordon equation is then set by the massless ground state W = 1/(1 + r²/8) ([Ibrahim, Masmoudi & Nakanishi](https://ar5iv.arxiv.org/html/1110.1709)).
- For the massless 4D cubic wave equation, smooth blowup by a concentrating W is known just above E(W) ([Hillairet & Raphaël](https://arxiv.org/pdf/1010.1768)).

## 7. Comparison at a glance

Five of the six phenomena rest on repeated returns to a decision point; oscillon lifetime resonances are tuned points instead. They differ in what brings the field back, what decides each return, and therefore how the windows are spaced.

| Phenomenon | Parameter | What brings the field back | What decides each return | Spacing | How many | Section |
| --- | --- | --- | --- | --- | --- | --- |
| Kink–antikink collisions | Speed v | Kinks fall back after the first bounce | Phase of the shape mode | Power law, v\_c² − v\_n² ∝ n⁻² | Infinite and nested | 2 |
| Kink or soliton on a defect | Speed | Trapping by the defect | Phase of the defect mode | Power law, n⁻²; widths n⁻³ | Infinite | 3 |
| Choptuik fine structure | Family parameter p | Echoes of the critical solution | Number of echoes before collapse | Periodic in ln(p − p\*), period 4.6 | A wiggle rather than windows | 4 |
| Oscillon lifetime resonances | Initial radius r₀ | Tuned points, not returns | One unstable mode per resonance | Isolated peaks | 125 resolved | 5 |
| Klein-Gordon scattering spikes | Amplitudes of the data | Oscillation about ±Q | Not analysed | Not measured | More appear on zooming; open whether infinitely many | 6 |
| AdS collapse staircase | Amplitude ε | Reflection from the boundary | Reaching Choptuik's solution after n reflections | Steps with ε\_n → 0 | As many as computed (10) | Appendix A |

## 8. Questions to think about

1. Suppose window n requires T(p\_n) = nT₀ + const and the return time diverges as T ∝ (p\_c − p)^{−a}. Show that p\_c − p\_n ∝ n^{−1/a}. What is a for kink–antikink collisions if you take p = v²?
2. Suppose the outcome is decided only at returns spaced T₀ apart, and near the edge the deviation from an unstable state grows like e^{λt}. Show that the window positions form a geometric sequence with ratio e^{λT₀}. What happens to the sequence if T₀ or λ drifts slowly from one return to the next?
3. Check that Hod and Piran's period 4.6 equals Δ/(2γ) for Δ = 3.44 and γ = 0.374. Why does a field that changes sign after half an echo give Δ/(2γ) rather than Δ/γ?
4. Why do three-bounce windows sit at the edges of the two-bounce windows in kink–antikink collisions?
5. The one-pass theorem forbids repeated passes near ±Q just above E(Q), yet a line through data space can cross a smooth threshold surface several times. How would you tell such crossings apart from windows produced by repeated returns, by looking at the solutions near the edges?
6. In the AdS staircase (Appendix A) every near miss collapses one round trip later. What would have to change for windows of dispersal to open between the steps?

## Appendix A. Collapse staircases in anti-de Sitter space (advanced)

In anti-de Sitter (AdS) space, a scalar pulse too weak to collapse at once comes back after reflecting off the boundary and collapses after n reflections. The outcome against amplitude is a staircase, and every step edge is a critical point.

**Background.** AdS is the maximally symmetric spacetime with a negative cosmological constant. Light rays reach its boundary in finite time and, with reflecting boundary conditions, come back, so AdS behaves like a box. In units where the AdS radius is 1, a light ray takes time π to go from the centre to the boundary and back.

**What Bizoń and Rostworowski found** ([Bizoń & Rostworowski 2011](https://arxiv.org/pdf/1104.3702)). They evolved a massless scalar field in 3+1 dimensional AdS from a Gaussian pulse of amplitude ε and fixed width.

- For large ε the pulse collapses on its first implosion, as in flat space.
- For smaller ε the pulse implodes, spreads, travels to the boundary, reflects and implodes again. Energy moves to shorter wavelengths with each pass, and a black hole forms after n reflections.
- The horizon radius as a function of ε is a sawtooth: continuous within a step, with a jump at each ε\_n.
- Each jump is a Choptuik critical point. Just above ε\_n, after n reflections, the solution approaches Choptuik's critical solution, and the horizon radius scales as (ε − ε\_n)^γ with γ ≈ 0.37.
- Successive steps are about π apart in collapse time, one round trip.
- The ε\_n decrease towards zero as n grows (computed up to n = 10). This suggests that arbitrarily weak pulses eventually collapse, so AdS would be unstable.

[Maliborski 2012](https://arxiv.org/abs/1208.2934) found the same staircase for flat space enclosed in a reflecting cavity. The reflecting wall, not the negative cosmological constant, produces the returns.

**How it compares.** The staircase has both ingredients of Section 1: returns (one light round trip, π) and a critical solution at every edge (Choptuik's self-similar solution, so the edges are Type II). It has no windows, though. The wall always sends the field back, so a near miss is never final: between neighbouring steps every pulse still collapses, one round trip later.

## References

**Klein-Gordon thresholds near the ground state**

- R. Donninger and W. Schlag, Numerical study of the blowup/global existence dichotomy for the focusing cubic nonlinear Klein-Gordon equation, Nonlinearity 24 (2011), [arXiv:1011.2015](https://arxiv.org/abs/1011.2015).
- P. Bizoń, T. Chmaj and N. Szpak, Dynamics near the threshold for blowup in the one-dimensional focusing nonlinear Klein-Gordon equation, J. Math. Phys. 52 (2011), [arXiv:1012.1033](https://arxiv.org/abs/1012.1033).
- K. Nakanishi and W. Schlag, Global dynamics above the ground state energy for the focusing nonlinear Klein-Gordon equation, [arXiv:1005.4894](https://www.arxiv.org/pdf/1005.4894).
- J. Krieger, K. Nakanishi and W. Schlag, Global dynamics away from the ground state for the energy-critical nonlinear wave equation, [arXiv:1010.3799](https://arxiv.org/pdf/1010.3799).

**Energy-critical case**

- S. Ibrahim, N. Masmoudi and K. Nakanishi, Threshold solutions in the case of mass-shift for the critical Klein-Gordon equation, Trans. Amer. Math. Soc. 366 (2014), [arXiv:1110.1709](https://ar5iv.arxiv.org/html/1110.1709).
- M. Hillairet and P. Raphaël, Smooth type II blow up solutions to the four dimensional energy critical wave equation, Anal. PDE 5 (2012), [arXiv:1010.1768](https://arxiv.org/pdf/1010.1768).

**Kink collisions**

- D. K. Campbell, J. F. Schonfeld and C. A. Wingate, Resonance structure in kink-antikink interactions in φ⁴ theory, Physica D 9, 1 (1983).
- P. Anninos, S. Oliveira and R. A. Matzner, Fractal structure in the scalar λ(φ²−1)² theory, Phys. Rev. D 44, 1147 (1991).
- R. H. Goodman and R. Haberman, Kink-antikink collisions in the φ⁴ equation: the n-bounce resonance and the separatrix map, SIAM J. Appl. Dyn. Syst. 4, 1195 (2005), [doi:10.1137/050632981](https://epubs.siam.org/doi/10.1137/050632981).
- R. H. Goodman and R. Haberman, Chaotic scattering and the n-bounce resonance in solitary-wave interactions, Phys. Rev. Lett. 98, 104103 (2007), [arXiv:nlin/0702048](https://arxiv.org/abs/nlin/0702048).
- R. H. Goodman, Chaotic scattering in solitary wave interactions: a singular iterated-map description, Chaos 18, 023113 (2008), [arXiv:0710.3209](https://www.arxiv.org/abs/0710.3209).
- Y. Zhu, R. Haberman and J. Yang, A universal map for fractal structures in weak solitary wave interactions, [arXiv:0802.3595](https://arxiv.org/pdf/0802.3595).
- P. Dorey, K. Mersh, T. Romanczukiewicz and Y. Shnir, Kink-antikink collisions in the φ⁶ model, Phys. Rev. Lett. 107, 091602 (2011), [arXiv:1101.5951](https://arxiv.org/pdf/1101.5951).
- P. G. Kevrekidis and R. H. Goodman, Four decades of kink interactions in nonlinear Klein-Gordon models: a crucial typo, recent developments and the challenges ahead (review), [arXiv:1909.03128](https://arxiv.org/pdf/1909.03128).

**Solitons and defects**

- R. H. Goodman and R. Haberman, Interaction of sine-Gordon kinks with defects: the two-bounce resonance, Physica D 195, 303 (2004), [arXiv:nlin/0310040](https://arxiv.org/pdf/nlin/0310040).
- R. H. Goodman, P. J. Holmes and M. I. Weinstein, Strong NLS soliton-defect interactions, Physica D 192, 215 (2004), [arXiv:nlin/0203057](https://www.arxiv.org/abs/nlin/0203057).

**Critical collapse**

- M. W. Choptuik, Universality and scaling in gravitational collapse of a massless scalar field, Phys. Rev. Lett. 70, 9 (1993).
- S. Hod and T. Piran, Fine structure of Choptuik's mass-scaling relation, Phys. Rev. D 55, R440 (1997), [arXiv:gr-qc/9606087](https://arxiv.org/abs/gr-qc/9606087).
- C. Gundlach, Critical phenomena in gravitational collapse (review), [arXiv:gr-qc/0210101](https://arxiv.org/pdf/gr-qc/0210101).
- C. Gundlach and J. M. Martín-García, Critical phenomena in gravitational collapse, Living Rev. Relativ. 10, 5 (2007), [arXiv:0711.4620](https://arxiv.org/abs/0711.4620).
- P. R. Brady, C. M. Chambers and S. M. C. V. Gonçalves, Phases of massive scalar field collapse, [arXiv:gr-qc/9709014](https://arxiv.org/pdf/gr-qc/9709014).
- S. H. Hawley and M. W. Choptuik, Boson stars driven to the brink of black hole formation, Phys. Rev. D 62, 104024 (2000), [arXiv:gr-qc/0007039](https://arxiv.org/abs/gr-qc/0007039).

**Oscillons**

- E. J. Copeland, M. Gleiser and H.-R. Müller, Oscillons: resonant configurations during bubble collapse, Phys. Rev. D 52, 1920 (1995), [arXiv:hep-ph/9503217](https://arxiv.org/abs/hep-ph/9503217).
- E. P. Honda and M. W. Choptuik, Fine structure of oscillons in the spherically symmetric φ⁴ Klein-Gordon model, Phys. Rev. D 65, 084037 (2002), [arXiv:hep-ph/0110065](https://arxiv.org/abs/hep-ph/0110065).
- M. Gleiser and M. Krackow, Resonant configurations in scalar field theories: can some oscillons live forever?, Phys. Rev. D 100, 116005 (2019), [arXiv:1906.04070](https://arxiv.org/abs/1906.04070).

**Anti-de Sitter space (Appendix A)**

- P. Bizoń and A. Rostworowski, Weakly turbulent instability of anti-de Sitter space, Phys. Rev. Lett. 107, 031102 (2011), [arXiv:1104.3702](https://arxiv.org/abs/1104.3702).
- M. Maliborski, Instability of flat space enclosed in a cavity, [arXiv:1208.2934](https://arxiv.org/abs/1208.2934).
