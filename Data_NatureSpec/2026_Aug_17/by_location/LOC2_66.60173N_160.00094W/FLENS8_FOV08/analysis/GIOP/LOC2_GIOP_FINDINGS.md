# LOC2 (2026-08-17), Selawik village, 66.60173 N 160.00094 W: GIOP findings

Hand-written interpretation. Numbers come from `giop_FINAL.csv` and `../REPORT.txt` of
the run with `--glint nir_similarity`. Physics caveats: `THEORY_GIOP_NOTE.md`.

## Input

- 7 water scans (00013-18 at the dock, 00033 near the bank), 6 overcast sky scans, one
  panel reference. 00019/00020 are water in the photos but sit on the land/water
  threshold, so they are excluded (see `../../../../PROCESSING_NOTES.md`).
- Glint: 00013 and 00033 were glint-correctable in step 5; NIR similarity applied to all.
- R_rs shape, 450-700 nm: 3.9 % (Aug-16 LOC1: 1.9 %). Amplitude 14.1 %.
- R_rs(443) 0.00084, R_rs(555) 0.00264, R_rs(665) 0.00395 sr⁻¹. Maximum at 687 nm.
- Spectral features (fig14): a dip near 670 nm between maxima at ~650 and ~690 nm is
  chlorophyll-a red absorption. There is a secondary maximum near 810 nm at the
  water-absorption minimum (particles). 00033 runs ~25 % high but keeps the same shape.

## What the inversion says

| config | M_φ | a_dg(443) m⁻¹ | b_bp(443) m⁻¹ | S_dg nm⁻¹ | η | RMS |
|---|---|---|---|---|---|---|
| constrained (S_dg 0.018) | 3e-26 (pinned at 0) | 5.48 | 0.044 | 0.018 | 0.22 | 29.8 % |
| free | 0.46 | 2.66 | 0.046 | 0.0084 | −1 (bound) | 13.1 % |
| max freedom | 0.68 | 2.65 | 0.046 | 0.0084 | −1 (bound) | 13.3 % |

1. **The standard GIOP-DC configuration cannot fit this water.** It reaches χ²_ν = 47
   and puts all absorption in CDOM, driving M_φ to zero. Every Monte Carlo draw in
   giop4 lands on the same bound.
2. **The fit needs a much flatter CDOM slope.** S_dg is 0.0084 (mean fit) and
   0.0068 ± 45 % (per scan), against the 0.014-0.020 range usually seen. Possible
   reasons, not separated here: humic terrestrial DOM, detrital (NAP) absorption folded
   into a_dg (detritus slopes are flatter), or model error the slope absorbs.
3. **a_dg(443) ≈ 2.7 m⁻¹ is the quotable number.** Free and max-freedom agree to 0.4 %.
   Per-scan fits give 4.3 ± 36 %, so state ~2.7 with a factor-1.5 uncertainty. That is
   3.4-8x the Aug-16 Kotzebue stations (0.33-0.78 m⁻¹).
4. **b_bp(443) ≈ 0.046 m⁻¹ from the mean fit**, similar to Kotzebue LOC1/LOC2a.
   Per-scan fits scatter (0.19 ± 62 %) and η sits at its −1 bound, so the backscatter
   spectral shape is not constrained.
5. **Do not quote chlorophyll.** M_φ goes from 0 to 0.68 across configurations and the
   seed self-consistency test has no fixed point. The 670 nm dip shows phytoplankton is
   present but not how much.

## Bottom line

Humic, CDOM-dominated water: a_dg(443) ≈ 2.7 m⁻¹ with an unusually flat slope.
Backscatter is moderate and chlorophyll is unconstrained. Even the best fit misses by
13 % RMS, so treat the composition as model-dependent.
