# Every validated station, one field day — cross-location synthesis

2026-08-16, Kotzebue, Alaska. Six validated R_rs products (LOC1, LOC2a, LOC2b, LOC3
FIBR15-main, LOC3 FIBR15-murky, LOC3 FLENS8) plus one non-water reflectance product
(LOC2c, algae). Figures: `by_location/COMPARISON/all_stations_overplot.png` (R_rs),
`by_location/COMPARISON/giop_cross_station_all.png` (GIOP composition).

## 1. R_rs: LOC1/LOC2 are optically deep; LOC3 is not

LOC1 and LOC2 (both sub-populations) collapse smoothly toward zero past ~700 nm, as
optically deep water must. **LOC3 does not** — every LOC3 sub-station shows a second,
distinct peak near 805–810 nm and retains substantial R_rs out to 900 nm. Confirmed in
the raw target radiance (not an R_rs-processing artefact): LOC3's L(810)/L(700) = 0.416
against LOC1's 0.226 — see `LOC3_66.89235N_162.59149W/LOC3_BOTTOM_CAVEAT.md`. **LOC3 is
a structurally different measurement regime from LOC1/LOC2 and its numbers should not be
pooled with theirs as equally certain water-column retrievals.**

## 2. GIOP composition: three real, distinct water states at LOC1/LOC2, plus LOC3's caveat

| | a_dg(443), free | b_bp(443), free | reading |
|---|---|---|---|
| LOC1 | 0.78 m⁻¹ | 0.043 m⁻¹ | its own site, more CDOM than either LOC2 population |
| LOC2a | 0.49 m⁻¹ | 0.041 m⁻¹ | main open water |
| LOC2b | 0.33 m⁻¹ | 0.047 m⁻¹ | disturbed/turbid — b_bp up, a_dg down (§`LOC2b_GIOP_FINDINGS.md`) |
| LOC3 (all 3) | 1.10–1.36 m⁻¹ | 0.079–0.122 m⁻¹ | elevated everywhere — consistent with a bright substrate adding apparent absorbing/scattering material GIOP has no other way to explain |

LOC3's a_dg and b_bp are the highest of any station **at every one of its three
sub-populations**, including the "clean" ones. Given §1, the more likely reading is not
"LOC3's water column carries more CDOM and particulates than LOC1/LOC2" but "GIOP is
partly absorbing a bright, shallow bottom into its water-column parameters, because that
is the only place in the model such a signal can go." This cannot be resolved without a
depth measurement (`LOC3_BOTTOM_CAVEAT.md`).

**Within LOC1/LOC2, which do not carry that caveat**, the composition differences are
real and traceable to real, independently-verified physical differences: LOC2b's water
was directly observed to be recently disturbed sediment (timing evidence,
`LOC2_SPLIT.md`), and its elevated b_bp is corroborated by GIOP for the third time
(spectral clustering, the glint-collapse test, and now composition, in three unrelated
pieces of code).

## 3. Fit quality: GIOP fits LOC2 much better than LOC1, and LOC3 worse still — for
   understood, different reasons

| | best RMS misfit (max freedom) |
|---|---|
| LOC1 | 8.8% |
| LOC2a / LOC2b | 2.9% / 2.9% |
| LOC3 FIBR15 / murky / FLENS8 | 5.1% / 2.8% / 6.0% |

LOC1's misfit is dominated by the blue band (12.75% RMS there alone, worst +38%),
traced to its high a_dg amplitude making a wrong CDOM-slope shape assumption costly (see
`LOC1_GIOP_FINDINGS.md`). LOC2's much better fit is real and reproduces per-scan (12/12
and 9/9 and 3/3 scans improve under the free configuration, `giop11_chi2_crossplot.png`
at each station). LOC3's numbers are not directly comparable to either, given §1 — a
"good fit" there may mean GIOP is fitting a bottom-contaminated spectrum well, not that
the water-column physics is better resolved.

**χ²_ν must never be compared across these six stations directly** — it scales with each
station's own measured shape uncertainty (as small as 0.16% at LOC2b, n=3 near-duplicate
scans, vs 1.93% at LOC1, n=12 genuinely independent replicates), so a χ²_ν of 16,000 at
LOC2b and 74 at LOC1 does NOT mean LOC2b fits 200× worse — its RMS misfit is actually
*better* (9.7% vs 10.9%). Full arithmetic in `THEORY_GIOP_NOTE.md` §4 point 4.

## 4. What generalizes beyond this specific dataset

- **A profile-solver bug (fixed this session) can look exactly like "the data can't
  determine the shapes."** Always check the nesting property (free ≥ as good as
  constrained) before trusting a "shapes are unconstrained" conclusion.
- **χ²_ν is not a station-comparable statistic** when σ is measured per station rather
  than assumed fixed. Use RMS misfit, or χ²_ν only within one station's own arm sweep.
- **A raw-radiance ratio (L(NIR-far)/L(NIR-near)) is a cheap, processing-independent
  screen for bottom contamination** in any above-water R_rs dataset — worth checking
  routinely, not only when photos happen to suggest it.
- **Photographic evidence and a quantitative population test can disagree**, and when
  they do (LOC3-FLENS8's real-but-photographically-ambiguous 2+2 substructure), the
  honest move is to document both rather than force one to win.

## 5. What is NOT yet resolved

- LOC3's water-column-vs-bottom split (needs a depth measurement — the single highest-
  value fix identified this session, cheaper than the wind-speed or footprint fixes
  below).
- The footprint/FOV question LOC3 was meant to answer (`LOC3_FOOTPRINT_COMPARISON.md`)
  — confounded by range, time, and the bottom-reflectance question all at once.
- Whether LOC2a's b_bp/a_dg would shift further under a full, correctly-weighted
  re-optimisation of `fit_shapes` at the grid optimum (§`LOC2a_GIOP_FINDINGS.md` A3: the
  solver lands ~5% off the true χ² grid minimum).
