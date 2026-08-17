# GIOP on LOC3 — three sub-populations, one site, and a caveat that overrides all of them

Covers all three LOC3 GIOP runs: **FIBR15 main** (n=3, clean open water), **FIBR15
murky** (n=2, 00058/00059, the high-sediment pair), **FLENS8** (n=4, kept as one
population). Same pipeline, same rigor as `LOC1_GIOP_FINDINGS.md` /
`LOC2a_GIOP_FINDINGS.md` / `LOC2b_GIOP_FINDINGS.md` — read LOC1's §A first for the
physics and method; this file states only what's specific to LOC3.

> ## ⚠ READ THIS FIRST: `LOC3_BOTTOM_CAVEAT.md`
>
> All four FLENS8 water photographs show the substrate directly through the water
> column (sharpest in 00047). Nothing comparable appears at LOC1 or LOC2. **Every number
> below may be a mixture of water-column and bottom signal**, and GIOP — like every
> semi-analytical water-column model — has no way to separate them without a depth
> measurement, which this station does not have. Treat every retrieval on this page as
> conditional on optical depth, not as a confirmed water-column composition.

---

## A. THE HEADLINE TABLE

| sub-station | config | M_φ | a_dg(443) | b_bp(443) | χ²_ν (ν=298) | RMS |
|---|---|---|---|---|---|---|
| FIBR15 main | constrained | 31.5 | 1.709 | 0.1945 | 1745 | 16.7% |
| FIBR15 main | free | 1.81 | 1.100 | 0.0792 | 19.4 | 4.9% |
| FIBR15 main | max freedom | 1.22 | 1.086 | 0.0761 | 15.0 | 5.1% |
| **FIBR15 murky** | constrained | 71.8 | 1.131 | **0.3453** | 140 | 22.6% |
| **FIBR15 murky** | free | 1.32 | 1.356 | **0.1222** | 1.7 | 2.8% |
| **FIBR15 murky** | max freedom | 1.35 | 1.355 | **0.1223** | 1.7 | 2.8% |
| FLENS8 | constrained | 31.8 | 2.326 | 0.2178 | 78.8 | 14.5% |
| FLENS8 | free | 4.86 | 1.223 | 0.0908 | 8.7 | 5.2% |
| FLENS8 | max freedom | 5.61 | 1.156 | 0.0887 | 6.0 | 6.0% |

Remember χ²_ν is not comparable across sub-stations here either — the murky pair's
measured σ (4.79%) is nearly double FIBR15 main's (2.62%), so its already-tiny χ²_ν is
partly a σ effect, not purely a better fit — though its RMS (2.8%) genuinely is the best
of the three, and honestly, that is more consistent with n=2 near-duplicate scans giving
the model little to disagree with than with the fit being unusually good physics.

## A1. b_bp confirms the murky pair independently — again

**b_bp(443) = 0.34 m⁻¹ (constrained) / 0.12 m⁻¹ (free) at the murky pair, against
0.076–0.22 m⁻¹ at FIBR15 main and FLENS8.** This is the third independent line of
evidence for that population (after the spectral-angle clustering and the failed glint
collapse test in `analyse_water_scans.py`): elevated particulate backscatter is exactly
what a high-sediment target should show, in a completely different piece of code from
the one that first flagged it.

## A2. What to quote, with the bottom caveat attached to all of it

| quantity | FIBR15 main | FLENS8 | murky pair |
|---|---|---|---|
| a_dg(443) | ≈1.09–1.10 m⁻¹ (agree to 1%) | ≈1.16–1.22 m⁻¹ (agree to 5%) | ≈1.36 m⁻¹ (agree to <1%, n=2) |
| b_bp(443) | ≈0.076–0.079 m⁻¹ (agree to 4%) | ≈0.089–0.091 m⁻¹ (agree to 2%) | ≈0.122 m⁻¹ (agree to <1%, n=2) |
| M_φ | do not quote (1.2–1.8, unstable pattern as everywhere) | do not quote (4.9–5.6) | do not quote |

By the "does the free/max-freedom pair agree" test used at every other station, all
three LOC3 sub-populations pass — a_dg and b_bp are internally consistent retrievals.
**But "internally consistent" is not the same claim as "this is the water column",
given §the caveat above.** A bright, shallow bottom under clear water would also produce
a stable, well-fitting retrieval — stability tests the inversion's self-consistency, not
whether the light came from the place the model assumes.

## A3. Self-consistency: FIBR15 main and FLENS8 pass; the murky pair does not converge at all

- FIBR15 main: stable fixed point chl=24.9 vs OC4=8.3 — agreement only to a factor 3,
  the worst of any station in the dataset including LOC2's.
- FLENS8: stable fixed point chl=24.8 vs OC4=11.1 — similarly poor, factor ~2.2.
- **The murky pair has no stable root in the swept range at all** — only an unstable one
  at chl≈1.0. Read together with A1's very low χ²_ν, this is consistent with a spectrum
  that GIOP's amplitude solver fits precisely at n=2 without the shape-consistency check
  (self-consistency) ever landing anywhere — a further, independent reason not to trust
  M_φ from this population even by the standard already applied everywhere else.

## A3b. GIOP's "good" fit says nothing about the bottom question — it never sees the evidence

**GIOP fits only 400–700 nm** (Bricaud's a*_φ table stops there — `THEORY_GIOP_NOTE.md`
§5). The single most diagnostic piece of evidence for shallow water at this site — the
second R_rs peak near 805–810 nm — sits entirely **outside** that window. Visually
confirmed in `giop6_all_fits.png` at every LOC3 sub-station: the fit panels stop at
700 nm and never show the 800 nm feature at all. So the RMS numbers in §A (2.8–6.0%,
better than LOC1's 8.8%) describe how well GIOP's smooth basis reproduces the visible-
light shape — a fit can be "good" there while being completely blind to the part of the
spectrum that actually distinguishes bottom-influenced water from deep water. **A good
GIOP fit is not evidence against the bottom-reflectance caveat; it simply isn't a test
of it.**

What GIOP *does* fit less well, visible in every one of FLENS8's 4 individual panels: a
consistent systematic notch around 660–680 nm that even the free configuration
undershoots (compare the green FREE curve to the blue measured curve near the chlorophyll
red-absorption band). The same region misfits at LOC1 too (its −20σ residual notch,
`LOC1_GIOP_FINDINGS.md` §A7), so this may be a structural GIOP limitation general to this
whole dataset rather than something specific to LOC3's optical-depth question.

## A4. Why this matters for LOC3 as a location

Every diagnostic in this file — the self-consistency failures, the higher measured σ,
the elevated-but-uncertain b_bp — is at least as consistent with **shallow, optically
mixed water** as with a real water-column composition difference. This dataset cannot
resolve which. The single highest-value fix, repeated from `LOC3_BOTTOM_CAVEAT.md`
because it matters this much: **measure depth at this site before drawing any water-
composition conclusion from it.**

## Files

Each sub-station's `analysis/GIOP/giop_FINAL.csv`, `giop_assumption_arms.csv`,
`giop_per_scan.csv`, figures `giop0`–`giop11`, and an auto-generated `README.md` index
(`make_giop_figures.py::write_giop_index`, regenerated fresh every run — do not
hand-edit those).
