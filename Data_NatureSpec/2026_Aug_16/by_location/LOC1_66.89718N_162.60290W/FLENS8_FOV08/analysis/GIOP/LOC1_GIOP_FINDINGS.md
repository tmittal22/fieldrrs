# GIOP on LOC1 — what the concentrations are worth

Input: `FINAL_Rrs.csv`, the amplitude-normalised mean of 12 angle-matched scans,
Kotzebue 66.89718 N 162.60290 W, 2026-08-16 20:24–20:41 UTC. Physics and symbol
definitions in `THEORY_GIOP_NOTE.md`; how to reproduce in `FIELD_DAY_WORKFLOW.md` §6.

> **This file grew by correction.** §A below is the current answer, superseding
> everything under it. §B and §C are the two corrections that got us here, kept because
> each names a mistake worth not repeating. §1–§5 at the bottom are the ORIGINAL
> six-band analysis and are **superseded in almost every number** — kept for provenance
> only, not for quoting.
>
> **⚠ 2026-08-18: `FINAL_Rrs.csv` was regenerated** with `--glint nir_similarity`
> (00005/00007 corrected — `PAPER_READINESS.md` §2d0). §A1's table and the "1.93%"
> uncertainty figure throughout this file are updated to the corrected run, freshly
> re-verified from `giop_FINAL.csv`. Everything else in §A2–§A5 (S_dg routes, per-scan
> gain, the χ² map) was NOT individually re-run against the correction — a_dg/b_bp moved
> <0.3% and RMS misfit <0.2 pp, so those numbers are very unlikely to have moved
> meaningfully, but treat them as "probably still right", not "reconfirmed", until they
> are.

---

## A. THE ANSWER

### A1. Nothing in the GIOP family fits this water

Against the **measured** per-band uncertainty (median 1.80 % of R_rs), ν = 298:

| configuration | free parameters | χ²_ν | RMS misfit |
|---|---|---|---|
| constrained (GIOP-DC: S_dg = 0.018, η from QAA) | 3 amplitudes | **48.0** | 10.9 % |
| free (S_dg, η fitted — still OC4-seeded) | 5 | **14.0** | 8.5 % |
| **maximum freedom** (+ a*_φ family/seed released) | 5 + seed profiled | **13.9** | 8.6 % |
| *a good fit would be* | | *≈ 1* | *≈ 1.8 %* |

The best model in the family still misfits by **4.7× the measurement uncertainty**, and
the residual is one smooth curve (lag-1 ρ = 0.9964), identical in all 12 scans. Every
number below is a projection onto a basis that cannot represent this spectrum.

### A2. What to quote, and what not to

| quantity | value | why |
|---|---|---|
| **a_dg(443)** | **0.78 m⁻¹** | stable across both free arms (0.779 → 0.781, 0.2 %); per-scan sd 7 % |
| **b_bp(443)** | **0.042 m⁻¹** | stable across both free arms (0.0430 → 0.0416, 3 %); per-scan sd 12 % |
| **S_dg** | **0.0113–0.0130 nm⁻¹** | a retrieval, not an assumption — §A3 |
| u(λ), b_b/a | see §3.1 | assumption-light: only the AOP–IOP operator and the air–water transfer |
| R_rs shape | 1.7 % | `THEORY_SCALED_MEAN.md` |
| M_φ / chlorophyll | **do not quote** | 11.5 constrained, 2.3 free, 1.7 max-free — a factor 7 across arms that all fit better than the one giving the largest value |
| η | **do not quote** | rails at the −1 bound; the b_bp power law is the wrong functional form |

⚠ The earlier recommendation to quote **a_dg = 1.25 ± 0.10 and b_bp = 0.085 ± 0.009 is
withdrawn.** Those are conditional on S_dg = 0.018, which the data reject.

### A3. S_dg is the best-constrained parameter in the model, and 0.018 is wrong here

χ² mapped over (S_dg, η) has near-vertical contours: a sharp interior minimum in S_dg,
and almost no curvature in η. Three routes agree:

| route | what was free | S_dg |
|---|---|---|
| `giop5`, 16-node sweep | S_dg only, η from QAA | 0.0130 |
| `giop9` panel b, 34×34 map | S_dg and η on a grid | 0.0113 |
| the solver, `fit_shapes=True` | S_dg and η, profiled | 0.0118 |
| **GIOP's own `sdg='obpg'` option** | — (a published band-ratio formula) | **0.0121** |
| GIOP-DC default | — | 0.0180 |

**The `obpg` line is the important one.** A published parameterisation that shares no
machinery with our χ² minimisation lands within **7 %** of it, while the DC default is
**59 % high**. Two independent routes agree this water has a shallower CDOM slope than
GIOP-DC assumes — consistent with a sediment- and detritus-rich Arctic river mouth.

### A4. The free-vs-constrained gain is not an averaging artefact

Per spectrum, independently (`giop11_chi2_crossplot.png`):

| | |
|---|---|
| constrained χ²_ν, median over 12 scans | 72.0 |
| free χ²_ν, median | 20.6 |
| improvement factor | **median 4.25×, range 3.16–4.91×** |
| scans improved | **12 of 12** |
| correlation of gain with R_rs(555) | r = +0.20 |

Nearly the same factor everywhere, no dependence on brightness: **the constrained shapes
are wrong in the same structural way in every spectrum.** The plot also checks the
solver — free is nested inside constrained, so no point may lie above the 1:1 line, and
that region is empty.

### A5. Chlorophyll: the internal check passes, and it does not help

GIOP uses the OC4 seed's *shape* and discards its *amplitude*, then fits M_φ free — so
whether M_φ reproduces its own seed is a real consistency test, one GIOP never runs
because it does not iterate. Sweeping the seed gives two roots: an **unstable** one at
chl 1.22 (M_φ collapsed onto its zero bound) and a **stable** one at **11.18**, against
OC4's 9.84 — **agreement to 14 %**.

So the old "factor-2.4 internal contradiction" is dead twice over (once by going
hyperspectral, once by this). But that agreement belongs to the **constrained** model,
which fits at χ²_ν = 48. The arms that actually fit want **M_φ ≈ 2.2–2.3, about 5× less
phytoplankton than OC4 reports**. Two Case-1 relations agreeing with each other is
evidence about their shared calibration, not about this water.

### A6. Which assumption actually costs the fit

Releasing S_dg and η: χ²_ν **48.0 → 14.0**. Additionally releasing the a*_φ family and
seed: **14.0 → 13.9**. So the phytoplankton prescription — the assumption most obviously
wrong on Case-2 water — is **not** what the misfit is made of. The CDOM/detritus slope is.

> ⚠ `fit_shapes=True` is **not** assumption-free: it builds a*_φ once from the OC4 seed
> and holds it fixed inside the profile. Only the maximum-freedom arm releases it.

### A7. Where the misfit lives, and what would fix it

The residual is a **+9σ lobe at 560–600 nm and a −20σ notch at 690 nm**, the same in all
12 scans. GIOP cannot reach past 700 nm at all — Bricaud's a*_φ table stops there — so
the two clearest features of this water, the ~700 nm particulate peak and the ~810 nm
water-absorption shoulder, are outside the inversion entirely. The instrument recorded
them; the model cannot use them.

The fix is a turbid-water inversion that extends past 700 nm (QAA-turbid, or `titanspec`'s
registry with a mineral component), plus an independent constraint on backscatter, since
η running to its bound says the b_bp power law cannot make the red-rising backscatter the
data want. **Not** more replicates: the shape is already determined to 1.7 %.

### Files

`giop_FINAL.csv` is the headline table. `giop_assumption_arms.csv` is every arm ranked by
χ². `giop_per_scan.csv` is the 12 individual fits. Figures indexed in `GIOP/README.md`.

---

> ## B. ⚠ CORRECTION 2 — "does anything actually FIT?", and an optimiser defect it exposed
>
> Everything below §2 originally quoted RANGES over an assumption sweep without ever
> saying how well each arm fitted. Some fit far worse than others, and pooling them
> produced a "spread" that was really a mixture of one surviving model and twenty
> rejected ones. Figures `giop9_chi2_weighting.png`, and χ² now annotated on
> `giop5` and `giop8`.
>
> ### 2a. No arm fits. Not one.
>
> Against the **measured** per-band uncertainty (median 1.9 % of R_rs):
>
> | | χ²_ν (ν = 298) | RMS misfit |
> |---|---|---|
> | best arm of all 24 (shapes free) | **18.1** | 8.5 % |
> | GIOP default (S_dg = 0.018, OC4 seed) | **74.5** | 10.9 % |
> | worst arm | **130.0** | 16.6 % |
> | *a good fit* | *≈ 1* | *≈ 1.9 %* |
>
> So the honest answer to "do we ever get a really good fit" is **no — the best model in
> the family still misfits by 4.5× the measurement uncertainty**. Every retrieved number
> in this file is a projection of a spectrum onto a basis that cannot represent it.
>
> ### 2b. χ²-weighting collapses the spread rather than narrowing it
>
> Weighting by exp(−Δχ²/2ŝ) with errors inflated so the best arm has χ²_ν = 1
> (Avni 1976 — the standard treatment when misfit is model inadequacy, not noise) still
> puts **w = 1.000 on one arm and 0.000 on the other 23**. The gap to the next arm is
> Δχ²_ν = 1.9 — which sounds small, but over 298 dof that is Δχ² = 544, and the gap to
> the best *fixed-shape* arm is Δχ² = 15 468. The sweep is not an uncertainty band.
> Cutting instead at a stated χ²_ν ≤ 2 × best (2 arms of 24 survive):
>
> | | unweighted range, as first reported | **admissible only** |
> |---|---|---|
> | M_φ | 0 → 41.2 | 0 → 2.26 |
> | a_dg(443) | 0.779 → 1.481 | **0.779 → 0.822** |
> | b_bp(443) | 0.0391 → 0.0843 | **0.0391 → 0.0430** |
>
> ### 2c. χ² cannot be used as a likelihood here, and the figure shows why
>
> The residual of the best fit has **lag-1 autocorrelation 0.9964 across 301 bands**. It
> is one smooth curve — a +9σ lobe at 560–600 nm and a −20σ notch at 690 nm — not noise.
> Under an AR(1) reading that is n_eff ≈ 0.5 independent points, i.e. the 301 bands
> constrain the residual *shape* but contribute essentially one independent statement
> about its size. **χ² here ranks models; it does not assign probabilities**, and the
> ±2.2 % Monte-Carlo error bar in §Correction 1 is formal precision on a rejected model.
> `giop6_all_fits.png` corroborates it directly: all 12 scans show the *same* misfit —
> GIOP overshoots the 570 nm peak and inverts the 690 nm upturn in **every one**.
>
> ### 2d. S_dg IS determined by these data. §2.4 below is superseded.
>
> §2.4 says "S_dg is not measured, it is set to 0.018 by convention" and sweeps it to get
> chlorophyll 0 → 297. That sweep is real, but it was run **without χ²**, and the arms at
> its ends fit terribly. Mapping χ² over (S_dg, η) on a 34 × 34 grid:
>
> - **S_dg has a sharp interior minimum at 0.0113–0.0118 nm⁻¹**, χ²_ν 75 → 17. The
>   contours are near-vertical: S_dg is the best-constrained parameter in the whole model.
> - The GIOP default **0.018 is not the preferred value for this water**. 0.0115 is at the
>   detritus-rich end of the usual 0.010–0.021 range, which is what a sediment-laden
>   Arctic river mouth should look like.
> - **η rails at the −1.0 bound** with a nearly flat χ² in that direction. Negative η
>   means b_bp *rising* toward the red, which a Junge/Mie power law does not do. Read that
>   as **the b_bp power-law form is wrong**, not as a measurement of η.
>
> Three estimates of S_dg appear across the figures and they agree; they differ only in
> what else was free at the time, and the spread is the honest resolution on it:
>
> | | what is free | S_dg | χ²_ν |
> |---|---|---|---|
> | `giop5`, 16-node sweep | S_dg only, η from QAA per arm | 0.0130 | 24 |
> | `giop9` panel b, 34×34 map | S_dg and η on a grid | 0.0113 | 17.2 |
> | the solver, `fit_shapes=True` | S_dg and η, profiled | 0.0118 | 18.1 |
>
> The solver lands ~5 % off the grid minimum, so the profile's polish leaves a little on
> the table; that is a real limitation of the implementation, not of the data.
>
> ### 2e. What this costs the headline numbers
>
> Fitting S_dg instead of assuming it moves the two "robust" parameters by more than any
> assumption swept before: **a_dg(443) 1.254 → 0.779 (−38 %)** and
> **b_bp(443) 0.0836 → 0.0430 (−49 %)**. The earlier recommendation to quote
> a_dg = 1.25 ± 0.10 and b_bp = 0.085 ± 0.009 is therefore **withdrawn**: those are the
> values conditional on S_dg = 0.018, and the data reject that S_dg. Quote
> **a_dg(443) ≈ 0.78–0.82 m⁻¹** and **b_bp(443) ≈ 0.039–0.043 m⁻¹**, stating that they
> come from the best-fitting arm of a family in which nothing fits.
>
> ### 2f. The defect this uncovered in our own package
>
> The first version of this analysis reported that freeing S_dg and η "RAILS both and
> collapses every amplitude — the data cannot determine them", at χ²_ν = 2431. **That was
> an optimiser failure, not a property of the data.** Freeing parameters cannot make the
> optimum worse, because the fixed-shape solution is a point inside the free search box,
> and 2431 against 74.5 is 33× worse. Three faults in six lines of
> `giop.inversion._invert_fmin_shapes`: `n_starts` silently ignored, a single hardcoded
> start at the upstream oligotrophic `[0.01, 0.001, chl]` (~100× off in two of five
> coordinates on turbid water), and a `return 1e6` barrier ~10¹¹ above the cost at the
> solution, onto which the simplex collapsed. Rewritten as a **profile** — amplitudes
> solved exactly by the bounded trust-region solver at each trial (S_dg, η), only the two
> shapes searched, the configured shapes always among the starts, so the returned cost is
> ≤ the fixed-shape cost by construction. `tests/test_fit_shapes_guard.py` pins the
> nesting property and fails 3.3× against the old solver. Suite 140 passed.

> ## C. ⚠ CORRECTION 1, and it inverts part of what this file first said
>
> The first version of this analysis ran GIOP on **6 bands** and reported its
> conditioning as though 6 bands were GIOP's limit. **That was wrong.** GIOP solves on
> whatever grid you give it; the Bricaud a*_φ table is continuous **400–700 nm at 2 nm**.
> I used 6 bands because that is what upstream `run_giop.m` demos.
>
> Running it hyperspectral changes the answer and the uncertainty by more than an order
> of magnitude:
>
> | | 6 bands | hyperspectral (301) |
> |---|---|---|
> | M_φ | 23.74 ± **52 %** | **11.47 ± 2.2 %** |
> | a_dg(443) | 2.171 ± 25.9 % | **1.254 ± 1.7 %** |
> | b_bp(443) | 0.158 ± 39.1 % | **0.084 ± 1.1 %** |
> | S_dg swing (0.014→0.022) | chl 0 → 297 | chl 2.4 → 20.9 |
>
> Two claims below are therefore **withdrawn**: that the retrieval amplifies input error
> ~27× (true only at 6 bands), and the "factor-2.4 contradiction" with OC4 (M_φ = 11.5
> against OC4's 9.84 is agreement to 17 %). The S_dg concern **stands**, reduced from
> fatal to largest-single-error.
>
> Also note 10 nm and 20 nm sampling **fail outright** — GIOP-DC anchors on R_rs(443)
> and needs a band within 2.5 nm of it.

---

---

# ⚠ SUPERSEDED — the original six-band analysis, kept for provenance only

Everything from here down was written against a **six-band** GIOP run with S_dg fixed at
0.018. Corrections B and C overturn most of it, and §A is the answer. Do not quote a
number from below without checking it against §A first: §1's amplitudes are ~2x high
(six bands), §2.2's "factor-2.4 contradiction" is withdrawn (§A5 measures 14 %
agreement), §2.3's "27x amplification" is withdrawn (true only at six bands), and §2.4's
"S_dg is not measured" is exactly backwards (§A3 measures it).

## 1. What came out

| | GIOP-DC | bounded |
|---|---|---|
| M_φ (chlorophyll amplitude) | 46.3 | **23.7** |
| a_dg(443) | 3.03 | **2.17** m⁻¹ |
| b_bp(443) | 0.259 | **0.158** m⁻¹ |

Neither railed, which is better than the raw station means (3 of 5 of those returned
−999 or values above 10¹³). Retrieved totals: **a(443) = 3.46 m⁻¹**, **b_b(443) = 0.160 m⁻¹**.

OC4 gives **9.84 mg m⁻³** chlorophyll.

---

## 2. Four reasons not to believe the concentrations

### 2.1 The model does not fit

RMS relative misfit **7.6 %** across the six bands, on a 3-parameter fit to 6 points:

| nm | measured r_rs | modelled | error |
|---|---|---|---|
| 412 | 0.002693 | 0.003084 | **+14.5 %** |
| 443 | 0.004731 | 0.004363 | −7.8 % |
| 490 | 0.008510 | 0.008100 | −4.8 % |
| 510 | 0.010600 | 0.010348 | −2.4 % |
| 555 | 0.016179 | 0.017204 | +6.3 % |
| 670 | 0.009870 | 0.009695 | −1.8 % |

The residual is **structured, not random** — too high at both ends, too low in the middle.
That is the signature of a basis that cannot represent the spectrum, not of noise. With
three free amplitudes and six bands there is little left to hide a misfit in.

### 2.2 A factor-2.4 internal contradiction

OC4 says chlorophyll is **9.84**; GIOP's own M_φ says **23.7**. These are the same
algorithm family applied to the same six numbers. (On the Titan Ridge plume the same
comparison gave a factor 30; this is milder but the same failure.)

Both are almost certainly wrong in the same direction: **OC4 and the Bricaud a*_φ are
Case-1 relations**, calibrated where phytoplankton co-varies with everything else. This
water is manifestly Case-2 — the spectrum peaks at 570 nm with a 700 nm sediment peak.
A blue/green ratio in sediment-dominated water reads suspended mineral as chlorophyll.

### 2.3 The retrieval amplifies the input uncertainty ~27×

Propagating the **1.7 % shape uncertainty** (300 Monte-Carlo draws) through the bounded solver:

| | value | uncertainty |
|---|---|---|
| OC4 chl | 9.98 | ±7.2 % |
| **M_φ** | 25.9 | **±46.7 %** |
| a_dg(443) | 2.27 | ±24.1 % |
| b_bp(443) | 0.169 | ±35.3 % |

A 1.7 % input becomes a 47 % output. That is the ill-conditioning of splitting one smooth
spectrum into three smooth basis functions, and it is a property of the problem, not of
the solver.

The **11 % amplitude** uncertainty behaves quite differently, and informatively:

| scale | OC4 | M_φ | a_dg(443) | b_bp(443) |
|---|---|---|---|---|
| ×0.89 | 9.84 | 24.19 | 2.191 | 0.1421 |
| ×1.00 | 9.84 | 23.74 | 2.171 | 0.1580 |
| ×1.11 | 9.84 | 23.24 | 2.150 | 0.1733 |

OC4 is **exactly invariant** — it is a band ratio, and scaling cancels. The absorption
terms barely move (2–4 %). Nearly all of it lands on **b_bp (±22 %)**, which is right:
R_rs ≈ b_b/a, so an overall scaling is read as backscatter. This is the quantitative
version of the earlier point that band ratios inherit the 1.7 % and magnitudes inherit
the 11 %.

### 2.4 An assumed constant controls the answer

`S_dg` is **not measured**. It is set to 0.018 nm⁻¹ by convention. Over its ordinary range:

| S_dg | M_φ | a_dg(443) | b_bp(443) |
|---|---|---|---|
| 0.014 | **0.00** (railed at zero) | 1.16 | 0.055 |
| 0.018 | 23.74 | 2.17 | 0.158 |
| 0.022 | **297.1** | 13.77 | 1.376 |

**Chlorophyll goes from 0 to 297 as a constant nobody measured moves across its plausible
range.** This single line is the reason concentrations from this retrieval are not
reportable. Everything in §1 is conditional on S_dg = 0.018 being right, and nothing here
tests that.

---

## 3. What IS defensible

### 3.1 u(λ), which is what R_rs measures

Before any library assumption, R_rs determines only $u = b_b/(a+b_b)$:

| nm | u | b_b/a |
|---|---|---|
| 412 | 0.0277 | 0.0285 |
| 443 | 0.0479 | 0.0503 |
| 490 | 0.0838 | 0.0915 |
| 510 | 0.1028 | 0.1146 |
| 555 | 0.1513 | 0.1783 |
| 670 | 0.0963 | 0.1065 |

This carries only the AOP–IOP operator (Gordon/Lee) and the air–water transfer, both
validated independently, and no assumption about what is *in* the water. **Report this.**
u ≤ 0.15 everywhere, so the water is absorption-dominated at every GIOP band despite
being visibly turbid.

### 3.2 The shape

Determined to **1.7 %** over 450–700 nm (§`THEORY_SCALED_MEAN.md`). Green peak at
570 nm, sediment peak at 700 nm, a_w window shoulder at 810 nm. Two spectra can be
compared on this footing without any inversion.

### 3.3 Relative statements

Because the shape is 7× better determined than the magnitude, **differences between
stations in shape** are far more defensible than differences in retrieved concentration.

---

## 4. The structural problem: GIOP stops at 670 nm

R_rs(700)/R_rs(670) = **1.17** — the 700 nm peak is the single clearest sediment
signature in this spectrum, and **GIOP's reddest band is 670 nm**. The 700 nm peak and
the 810 nm shoulder are entirely outside its band set, and its Bricaud a*_φ is undefined
above 700 nm so it cannot simply be extended.

The instrument measured 350–2500 nm hyperspectral. GIOP used **6 numbers** and threw away
the two features that most clearly identify this water. That is the argument for a
turbid-water inversion (QAA-turbid, or titanspec's registry with a mineral component)
rather than for tuning GIOP.

---

## 5. Bottom line — SUPERSEDED by §A2

This table was written before the maximum-freedom arm and the self-consistency test
existed. **Use §A2.** Two entries in the version that stood here are now wrong: it cited
the "factor-2.4 self-contradiction" as the reason not to quote chlorophyll (§A5 measures
14 % agreement — the reason is the factor-7 spread across arms, not a contradiction),
and it gave a_dg/b_bp as ranges over "admissible arms" rather than as the single stable
value the free arms agree on.
