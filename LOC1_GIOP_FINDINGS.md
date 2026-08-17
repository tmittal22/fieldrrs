# GIOP on LOC1 — what the concentrations are worth

Input: `FINAL_Rrs.csv`, the amplitude-normalised mean of 12 angle-matched scans,
Kotzebue 66.89718 N 162.60290 W, 2026-08-16 20:24–20:41 UTC.

**Summary: do not report concentrations from this. Report `u(λ)`, and the shape.**

---

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

## 5. Bottom line

| quantity | report it? |
|---|---|
| u(λ), b_b/a | **yes** — measured, assumption-light |
| R_rs shape | **yes** — 1.7 % |
| R_rs magnitude | yes, with the 11 % environmental spread stated |
| a(443), b_b(443) totals | with caution — conditional on the operator, ±25–35 % |
| a_φ / a_dg **split** | **no** — 47 % on a 1.7 % input, and S_dg swings it 0→297 |
| chlorophyll | **no** — Case-1 algorithm on Case-2 water, factor-2.4 self-contradiction |
