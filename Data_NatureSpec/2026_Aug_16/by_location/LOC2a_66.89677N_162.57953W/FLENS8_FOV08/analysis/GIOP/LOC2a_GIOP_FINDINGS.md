# GIOP on LOC2a — the main open-water population

Input: `FINAL_Rrs.csv`, the amplitude-normalised mean of 9 angle-matched scans (main
open-water group; `--glint nir_similarity` applied — see `../water_scans/
SCAN_00035_GLINT.md`), Kotzebue 66.89677 N 162.57953 W, 2026-08-16. Same pipeline, same
rigor, same figure set as `LOC1_GIOP_FINDINGS.md`; read that file's §A first for the
physics and the method — this file states only what differs at this station.

---

## A. THE ANSWER

### A1. Nothing in the GIOP family fits well here either — and it fits noticeably better than LOC1

Against the **measured** per-band uncertainty (median 0.77 %, tighter than LOC1's
1.93 % — LOC2a's 9 scans agree in shape more closely than LOC1's 12 do):

| configuration | free parameters | χ²_ν | RMS misfit |
|---|---|---|---|
| constrained (GIOP-DC) | 3 | 221.8 | 8.2 % |
| free (S_dg, η fitted) | 5 | 104.3 | 4.7 % |
| **maximum freedom** (+ a*_φ family/seed) | 5 + seed | **23.6** | **2.9 %** |
| *a good fit would be* | | *≈ 1* | *≈ 0.8 %* |

> ⚠ **Do not compare χ²_ν across stations.** LOC1's constrained χ²_ν was 74.5; LOC2a's is
> 221.8 — three times worse-looking, despite LOC2a's RMS misfit (8.2 %) being *better*
> than LOC1's (10.9 %). χ²_ν = Σ(Δ/σ)², and σ is each station's own **measured** shape
> uncertainty. LOC2a's 9 scans happen to agree in shape more tightly than LOC1's 12 do
> (0.77 % vs 1.93 % median σ), so the *same* relative misfit produces a *larger* χ²_ν —
> roughly (1.93/0.77)² ≈ 6×, which is the right order of magnitude for the 3× gap
> actually observed. **RMS misfit is the number that is comparable across stations; χ²_ν
> is not**, and this is the clearest demonstration of that in the whole dataset.

**Maximum freedom fits much better here than at LOC1, relatively speaking.** At LOC1,
releasing the a*_φ family bought almost nothing (χ²_ν 18.1 → 17.2, RMS 8.5 % → 8.8 %,
actually *slightly worse* RMS). At LOC2a it buys a real gain (χ²_ν 104.3 → 23.6, RMS
4.7 % → 2.9 %). The phytoplankton prescription matters more to the misfit at LOC2a than
it did at LOC1 — see A5.

### A2. What to quote, and what not to

| quantity | value | why |
|---|---|---|
| **a_dg(443)** | **≈ 0.49 m⁻¹** | free 0.4906, max-free 0.4868 — agree to 0.8 %; per-scan sd 6 % |
| **b_bp(443)** | **≈ 0.041–0.043 m⁻¹** | free 0.04103, max-free 0.04286 — agree to 4 %; per-scan sd 5 % |
| **S_dg** | **0.0107–0.0124 nm⁻¹** | a retrieval, see A3 |
| M_φ / chlorophyll | **do not quote** | 3.6 (free) → 4.2 (max-free), and constrained gives 8.7 — same instability pattern as LOC1, worse self-consistency (A5) |
| η | **do not quote** | −1.00 (free) / −0.81 (max-free), both near the bound |

a_dg and b_bp are both **lower** at LOC2a than at LOC1 (0.49 vs 0.78 m⁻¹; 0.041 vs
0.043 m⁻¹ — b_bp is close, a_dg is markedly lower). LOC1 and LOC2a are different
locations 2.6 km apart, not two states of the same water, so this is a genuine
site-to-site difference, not something the split resolved.

### A3. S_dg — again sharply determined, and again far from the GIOP-DC default

| route | S_dg | χ²_ν |
|---|---|---|
| GIOP-DC default | 0.0180 | 221.8 |
| χ² grid minimum (34×34) | 0.0133 | 100.6 |
| the solver, `fit_shapes=True` | 0.0124 | 104.3 |

Contours are again near-vertical in the (S_dg, η) map: S_dg is well constrained, η is
nearly flat and runs to its bound (−1.00 free, −0.81 at maximum freedom — the only one of
the three stations where η does **not** rail exactly at −1, worth noting though not a
qualitatively different reading).

### A4. The free-vs-constrained gain reproduces in every scan

| | |
|---|---|
| constrained χ²_ν, median over 9 scans | 250.0 |
| free χ²_ν, median | 105.6 |
| improvement factor | median 2.37×, range 2.05–2.66× |
| scans improved | **9 of 9** |

Smaller factor than LOC1's (median 4.25×), consistent with LOC2a's constrained fit
already being relatively closer to its free fit (less headroom to improve) than LOC1's.

### A5. Self-consistency is WORSE here than at LOC1

Sweeping the chlorophyll seed against its own retrieved M_φ: **stable fixed point at
chl = 10.12, against OC4's 14.35 — agreement to 29 %** (LOC1: 14 %). Two mechanisms,
both real:

- **OC4 itself moved.** LOC2a's OC4 chlorophyll is 14.35 mg m⁻³ against LOC1's 9.84 — a
  band-ratio number, not a retrieval, but part of the seed chain (§2.1 in
  `../../../../THEORY_GIOP_NOTE.md`).
- **The phytoplankton shape matters more here.** A5 above already showed releasing a*_φ
  buys a real χ² improvement at LOC2a that it did not at LOC1 — consistent with the
  seed/self-consistency also being less settled.

### A6. Which assumption costs the fit — DIFFERENT from LOC1

Releasing S_dg/η: χ²_ν 221.8 → 104.3 (2.1× — a real but *smaller relative* gain than
LOC1's 4.1×). Additionally releasing a*_φ: 104.3 → 23.6 (**a further 4.4× gain** — far
larger than LOC1's negligible 1.05×). **At LOC1 the CDOM slope was the dominant
misspecification; at LOC2a the phytoplankton shape is at least as large a factor.** This
is a genuine cross-station difference in what the model gets wrong, not the same finding
restated.

### A7. A real, universal residual feature near 530–545 nm

The mean-spectrum residual (`giop9_chi2_weighting.png` panel c) shows a sharp −60σ
feature there — far more extreme than LOC1's smooth ±20σ lobes. Checked against
`giop6_all_fits.png` before writing this: **the underlying deviation is small and present
in all 9 individual scans** (a slight shoulder GIOP's fixed shapes don't reproduce, most
visible as the FREE curve tracking the data closely everywhere except a narrow dip around
530–545 nm). The extreme σ-normalised value is because the 9 scans agree almost exactly
in shape at that specific wavelength (locally even tighter than the 0.77 % median), so a
small absolute mismatch there is inflated by an unusually small local σ. Real, small,
universal, unexplained — flagged rather than investigated further this session.

### A8. Cross-station comparison — see `../../COMPARISON/giop_cross_station.png`

b_bp rises modestly from LOC2a to LOC2b (the disturbed-water station, its own findings in
`../../LOC2b_.../analysis/GIOP/LOC2b_GIOP_FINDINGS.md`) in the direction more suspended
sediment would predict — but a_dg's direction is inconsistent between the free and
max-freedom arms at LOC2b specifically, so this is **suggestive, not confirmed**. Full
discussion in that file's §A8, since LOC2b's retrieval is the less stable of the two.

### Files

`giop_FINAL.csv`, `giop_assumption_arms.csv`, `giop_per_scan.csv` as at LOC1. Figures
`giop0`–`giop11`, indexed in `README.md`.
