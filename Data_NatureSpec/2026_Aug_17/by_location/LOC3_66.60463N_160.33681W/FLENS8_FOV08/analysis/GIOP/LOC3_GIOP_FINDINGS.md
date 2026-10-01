# LOC3 (2026-08-17), Selawik Lake west, 66.60463 N 160.33681 W: GIOP findings

Hand-written interpretation. Numbers come from `giop_FINAL.csv` and `../REPORT.txt` of
the run with `--glint nir_similarity`. Physics caveats: `THEORY_GIOP_NOTE.md`.

## Input

- 16 water scans (00035-41, 00049-57), 7 overcast sky scans, one panel reference;
  choppy, with whitecaps in 00052. The boat drifted 61 m (00057 alone at 60 m tolerance).
- Glint: 00038 and 00052 were glint-correctable in step 5; NIR similarity applied to all.
- R_rs shape, 450-700 nm: 6.8 %. Amplitude 24.3 %.
- R_rs(443) 0.00155, R_rs(555) 0.00325, R_rs(665) 0.00387 sr⁻¹. Broad maximum
  ~600-690 nm, peaking at 649 nm, with a 670 nm dip as at LOC2.
- **The blue end is the weak part.** At 400 nm per-scan R_rs runs from 0.0004 to
  0.0027 sr⁻¹ (6x), while 550-700 nm stays consistent. Blue is where L_w is smallest
  relative to ρ·L_sky, so sky-reflection error lands there. Scans 00049-57 sit higher
  than 00035-41.
- There is an unexplained step at ~525-530 nm in most scans and a weaker one at ~590 nm
  (also at LOC2). Not a known water feature, not interpreted.

## What the inversion says

| config | M_φ | a_dg(443) m⁻¹ | b_bp(443) m⁻¹ | S_dg nm⁻¹ | η | RMS |
|---|---|---|---|---|---|---|
| constrained (S_dg 0.018) | 3.5 | 3.68 | 0.048 | 0.018 | 0.44 | 35.2 % |
| free | 3.0 | 3.27 | 0.106 | 0.0078 | 0.48 | 3.9 % |
| max freedom (Ciotti sf=1) | 8.0 | 4.02 | 0.143 | 0.0076 | 0.89 | 2.2 % |

1. **The fixed CDOM slope fails here too**: χ²_ν = 63, RMS 35 %. Freeing S_dg gets RMS
   down to 3.9 %. The fitted slope, 0.0076-0.0078, matches LOC2. The flat slope shows
   up at two stations 15 km apart, so it is a property of this water or of the model in
   it, not one bad fit.
2. **a_dg(443) = 3.3-4.0 m⁻¹** (free vs max freedom); per-scan 5.3 ± 55 %. That is
   slightly above LOC2, and the scatter is larger because a_dg comes from the blue
   (see Input).
3. **b_bp(443) = 0.11-0.14 m⁻¹**, about 2.5x LOC2. LOC3 carries more particles, which
   fits the brighter 550-700 nm R_rs and the open, choppy lake.
4. **Do not quote chlorophyll.** M_φ runs 3.0-8.0 across configurations. The
   self-consistency fixed point is at chl 47.8 against OC4's 7.2, and the 8.2 fixed point
   is unstable.
5. **Depth was not measured.** If the lake is shallow enough for bottom reflectance to
   reach the sensor, part of b_bp and a_dg is bottom (compare Aug-16 `LOC3_BOTTOM_CAVEAT.md`).
   No photo shows the bottom, which argues against it but does not prove it.

## Bottom line

LOC3 has the same CDOM regime as LOC2 (a_dg(443) ~3-4 m⁻¹, flat slope) and roughly
2.5x more particulate backscatter. The fit is good (2-4 % RMS), but the composition
inherits the blue-end scatter, so a_dg carries about ±50 % per scan.
