# GIOP on LOC2b — the disturbed/turbid water at station start

Input: `FINAL_Rrs.csv`, the amplitude-normalised mean of **3** angle-matched scans
(00027/28/29 — see `../../LOC2_SPLIT.md` for why these are a separate, physically real
water state, not excluded data). Same pipeline as LOC1/LOC2a. Read `LOC1_GIOP_FINDINGS.md`
§A first.

> ⚠ **n = 3.** Every number below carries a weaker statistical basis than LOC1's (n=12)
> or LOC2a's (n=9), and §A6 below documents a real consequence of that: the retrieval is
> measurably less stable here than at either other station. This is not a reason to
> distrust the water-population finding (`analyse_water_scans.py`'s evidence for LOC2b
> being real, different water does not depend on GIOP at all) — it is a reason to hold
> the GIOP *composition* numbers from this station more loosely than LOC1's or LOC2a's.

---

## A. THE ANSWER

### A1. Nothing fits — and χ²_ν is enormous, for a reason that is NOT alarming

Against the measured per-band uncertainty (median **0.16 %** over 400–700 nm — the three
scans agree in shape to 0.53° spectral angle, far tighter than either other station,
because they were taken 90 seconds apart from the same settling plume):

| configuration | χ²_ν | RMS misfit |
|---|---|---|
| constrained (GIOP-DC) | 16 055 | 9.7 % |
| free (S_dg, η fitted) | 6 199 | 4.2 % |
| **maximum freedom** | **999** | **2.9 %** |
| *a good fit would be* | *≈ 1* | *≈ 0.3 %* |

**χ²_ν in the thousands does not mean the model fits catastrophically worse here than at
LOC1.** RMS misfit — the station-comparable quantity — is 9.7 % here against LOC1's
10.9 % and LOC2a's 8.2 %: **all three stations misfit by about the same relative
amount.** χ²_ν diverges because σ does: LOC2b's median σ over 400–700 nm (0.16 %) is
**~12× smaller** than LOC1's (1.93 %), and χ²_ν scales as 1/σ², so the same relative
misfit alone predicts roughly a **150×** larger χ²_ν here — the right order of magnitude
for the observed 16 055/74.5 ≈ 216× gap. The remaining factor is real (σ is not a single
flat number, it varies band to band, and LOC2b's RMS is slightly *better* than LOC1's,
which should if anything narrow the gap) but is not the point worth chasing further here:
**χ²_ν is a σ-dependent quantity and is not comparable across stations**, full stop —
that is the finding, not the exact multiplier.

### A2. What to quote, and the honest caveat on all of it

| quantity | value | caveat |
|---|---|---|
| **b_bp(443)** | **≈ 0.047–0.057 m⁻¹** | free 0.04678, max-free 0.05748 — a 23 % spread, LARGER than LOC1's (3 %) or LOC2a's (4 %) equivalent spread |
| a_dg(443) | **0.33–0.49 m⁻¹, NOT settled** | free 0.3297, max-free 0.4932 — a factor-1.5 disagreement. See A6. |
| S_dg | 0.0080–0.0084 nm⁻¹ | lowest of the three stations |
| M_φ / chlorophyll | **do not quote** | as everywhere else, and worst self-consistency of the three stations (A5) |

**a_dg is the parameter that broke here.** At LOC1 and LOC2a, a_dg agreed to within 1–4 %
between the free and maximum-freedom arms — the clearest "this is a robust retrieval"
signal in the whole dataset. At LOC2b it does not: 0.33 vs 0.49 m⁻¹, a 50 % swing. b_bp
also disagrees more than at the other two stations (23 % vs 3–4 %). **Do not report a
single a_dg value for LOC2b as if it carried LOC1's or LOC2a's confidence.**

### A3. S_dg is still sharply determined — this part of the result is robust

| route | S_dg | χ²_ν |
|---|---|---|
| GIOP-DC default | 0.0180 | 16 055 |
| grid minimum | 0.0100 | 4 612 |
| solver | 0.0084 | 6 199 |

Lower than either LOC1 (0.011–0.013) or LOC2a (0.011–0.013) — consistent with a
sediment-fresh, detritus-rich transient state having a shallower CDOM/detrital slope than
the settled water either main population shows. This is the one parameter that stayed
well-behaved across all three stations' analyses.

### A4. The free-vs-constrained gain reproduces in all 3 scans

| | |
|---|---|
| constrained χ²_ν, median | 15 986.7 |
| free χ²_ν, median | 6 170.4 |
| improvement factor | median 2.59×, range 2.52–2.66× |
| scans improved | **3 of 3** |

The narrow range (2.52–2.66×, tighter than LOC1's 3.16–4.91× or LOC2a's 2.05–2.66×) is
expected from n=3 near-duplicate scans (0.53° apart), not evidence of anything
qualitatively different about the fit.

### A5. Self-consistency is the worst of the three stations

Stable fixed point chl = 14.82, against OC4's 10.77 — **agreement to 38 %** (LOC1: 14 %,
LOC2a: 29 %). Consistent with a monotonic pattern: self-consistency degrades as n drops
and the misfit (in σ-normalised terms) grows. Not a new independent finding, a restatement
of A1/A2 at the level of one specific diagnostic.

### A6. Why a_dg specifically is unstable here, and it is NOT a bug

`giop_python/tests/test_fit_shapes_guard.py`'s nesting property (free ≤ constrained cost
by construction) is confirmed to hold here too — this is not the profile-solver defect
from earlier in the session recurring. The instability is a property of the **data**: at
n=3, tightly clustered scans (0.53°), the fit has far less independent information to
separate a_dg from the other two amplitudes than at n=9 or n=12, and A1 already
establishes the model does not fit well at any of the three stations — with 5–6 free
shape+amplitude parameters chasing a fixed, badly-fitting basis, the exact optimum a_dg
settles on is more sensitive to exactly which shapes are free at LOC2b than it is at the
better-replicated stations.

### A7. Cross-station comparison, and the honest read of it

`../../COMPARISON/giop_cross_station.png`, numbers in `LOC2a_GIOP_FINDINGS.md` §A8:

| | LOC2a (free / max-free) | LOC2b (free / max-free) |
|---|---|---|
| b_bp(443) | 0.041 / 0.043 | 0.047 / **0.057** |
| a_dg(443) | 0.49 / 0.49 | 0.33 / **0.49** |

**b_bp is higher at LOC2b than LOC2a in both configurations** — the direction more
suspended/disturbed sediment predicts, and consistent with the independent finding in
`../../LOC2_SPLIT.md` (elevated blue+NIR in the raw R_rs, timing evidence for disturbed
sediment). **This is suggestive corroboration from GIOP, not confirmation**: the
max-freedom b_bp gap (0.043 → 0.057, +33 %) is larger than the free-arm gap (0.041 → 0.047,
+14 %), and A2/A6 above establish that LOC2b's retrieval is the least stable of the three
— a real physical difference and a weaker retrieval both push in the same direction here,
and this analysis cannot fully separate them with n=3. a_dg gives no clean signal either
way (0.33 vs 0.49 free-arm values move in opposite directions depending on configuration).

**The R_rs-level finding stands on its own regardless.** The elevated blue/NIR signature,
the timing evidence, and the failure of every artefact test (§`LOC2_SPLIT.md`) do not
depend on GIOP at all. GIOP's b_bp result is a second, weaker, independently-suggestive
line of evidence — not the basis for the LOC2b finding, which was already established
before any inversion ran.

### Files

`giop_FINAL.csv`, `giop_assumption_arms.csv`, `giop_per_scan.csv` (n=3 rows).
