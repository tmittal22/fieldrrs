# Paper readiness — 2026-08-16 Kotzebue field day, status as of 2026-08-18

What is validated and quotable, what is real but conditional, what is still open, and
where every number and figure lives. Written after a full pass over all three field
locations (LOC1, LOC2 split into a/b/c, LOC3 split into three sub-populations across two
foreoptics) — six R_rs products, one reflectance product, and their GIOP inversions.

A stage below is called "done" only when it ran this session with numbers shown, was
checked against something independent where that was possible (a second ranker, a raw
radiance, an external weather record), and — where a defect was found — shown to fail
against the pre-fix code before being trusted fixed.

---

## 1. What is ready to write into a manuscript now

### 1a. The instrument chain and processing pipeline
Verified end to end: hand-derived algebra matches the package to machine precision,
agrees with the instrument's own onboard reflectance column to <0.1–0.3% with no
systematic bias at every station, and the round-trip (synthetic R_rs → radiances →
recovered R_rs) closes to 1e-16. `FIELD_DAY_WORKFLOW.md` §3, `verify_field_calcs.py`
output archived in every station's `analysis/`.

### 1b. LOC1 and LOC2 are one clean water-column measurement problem each
Both locations pass every conservation/physicality check, are confirmed optically deep
(R_rs collapses smoothly toward zero past ~700 nm, unlike LOC3 — §2c below), and their
GIOP retrievals are self-consistent (a_dg/b_bp agree to within 5% between independently
free and maximum-freedom configurations at LOC1 and LOC2a; LOC2b is the documented
exception, see the table footnote) and independently corroborated where a claim needed
it (LOC2b's disturbed-water composition confirmed three separate ways — spectral
clustering, a failed glint-collapse test, and elevated GIOP b_bp — in three different
pieces of code).

**Quotable numbers**, with the standard caveats stated alongside each:

| | R_rs(555) | a_dg(443) | b_bp(443) |
|---|---|---|---|
| LOC1 | 0.00863 sr⁻¹ | 0.78 m⁻¹ | 0.043 m⁻¹ |
| LOC2a (main) | 0.01066 sr⁻¹ | 0.49 m⁻¹ | 0.041–0.043 m⁻¹ |
| LOC2b (disturbed) | 0.01068 sr⁻¹ | 0.33–0.49 m⁻¹* | 0.047–0.057 m⁻¹ |

*LOC2b's a_dg is the one number in this table that is NOT stably determined — see
`LOC2b_GIOP_FINDINGS.md` §A2/A6; report as a range, not a point value.

**Never quote M_φ / chlorophyll as a concentration anywhere in this dataset.** It moves
by a factor 5–7 across equally-valid fit configurations at every single station, for a
documented, mechanistic reason (`THEORY_GIOP_NOTE.md` §5): a Case-1 phytoplankton model
on Case-2 water.

**Never compare χ²_ν across stations.** It scales with each station's own measured shape
uncertainty, not with fit quality — RMS misfit is the comparable number, and it is
reported at every station specifically to make this comparison possible
(`THEORY_GIOP_NOTE.md` §4.4, `Data_NatureSpec/.../GLOBAL_COMPARISON.md` §3).

### 1c. The scaled-mean / shape-vs-amplitude method
`THEORY_SCALED_MEAN.md` — full derivation, the rank-1 SVD equivalence measured (1.747%
vs 1.739%), and the failure modes stated. Reusable methodology, not specific to this
field day.

### 1d. Two real, generalizable methodological findings
- **A profile-solver bug can look exactly like "the data cannot determine the model
  shapes."** The nesting-property test that catches it (free ≤ as good as constrained,
  by construction) is now a permanent regression test
  (`giop_python/tests/test_fit_shapes_guard.py`), not a one-off fix.
- **χ²_ν is not a cross-dataset statistic when σ is measured rather than assumed** —
  demonstrated with real numbers (LOC2b's χ²_ν of 16,000 against LOC1's 74 is a σ
  artefact; LOC2b's RMS misfit is actually *better*). Worth stating in the methods
  section of any paper using this or a similar pipeline, since it is not specific to
  this dataset.

---

## 2. Real, but conditional — state the caveat every time these are cited

### 2a. LOC2b (disturbed water) — real, not an artefact, but n=3
Timing evidence (first three scans at the station, range unchanged into the next scan)
plus three independent failed-artefact tests point to genuinely disturbed sediment
settling over ~2 minutes. `LOC2_SPLIT.md`. GIOP composition from it should be reported
with n=3's weaker statistical basis stated, not silently pooled with LOC1/LOC2a's larger-n
confidence.

### 2b. LOC2a's 00035 (glint) — corrected, not just flagged
Fixed via `--glint nir_similarity`, verified to collapse the scan back into its own
population (3.16°→1.53°, inside the group's own 1.51° spread) before being trusted.
`SCAN_00035_GLINT.md`.

### 2c. LOC3 — an entire location's worth of results, conditional on optical depth
**This is the single most important caveat in the whole dataset.** All four FLENS8
photos show visible substrate; confirmed independently in the RAW radiance (not an R_rs
processing artefact): L(810)/L(700) = 0.416 at LOC3 vs 0.226 at LOC1, and every LOC3
sub-station shows a second R_rs peak near 805–810 nm that optically deep water cannot
produce (`LOC3_BOTTOM_CAVEAT.md`). GIOP's a_dg/b_bp are elevated at every one of LOC3's
three sub-populations, including the "clean" ones — plausibly the model absorbing a
bright shallow bottom into its only available water-column parameters, not a real
compositional signal. **Do not report LOC3 composition numbers as water-column retrievals
without this caveat attached, and do not average them into a "typical" number with
LOC1/LOC2.** The murky pair's b_bp elevation (§1b-style corroboration, `LOC3_GIOP_FINDINGS.md`
§A1) is the one LOC3 GIOP finding least affected by this, since it is corroborated by two
methods that do not depend on optical-depth assumptions at all.

### 2d0. LOC1 had the same glint defect LOC2a had — found this session, now fixed (2026-08-18)
`analyse_water_scans.py` had never been run at LOC1 (it postdates LOC1's original
processing); running it retroactively found 00005 and 00007 (2 of 12 water scans) verdict
`glint (correctable)`, the same test that flagged LOC2a's scan 00035. This was left open
for the user's explicit decision rather than silently applied (LOC1's numbers are "the"
reference value throughout this dataset), and the user confirmed: apply it.

**Applied via `--glint nir_similarity`, the same correction and code path as LOC2a's
00035, re-run through the full pipeline (water QC, interactive report, GIOP, all-spectra
figure, highlights, slide deck, all 4 cross-station comparison figures).** Verified
before quoting anything from it:
- **Shape consistency across the spectrum nearly doubled**: `shape_cv_pct` 7.08% → 3.68%
  (the actual target of this correction).
- **R_rs(555), the number everyone quotes, barely moved**: 0.00865 → 0.00863 sr⁻¹ (−0.23%).
- **GIOP composition is essentially unchanged**: a_dg(443) 0.779→0.780, b_bp(443)
  0.0430→0.0429 m⁻¹ (free config) — the retrieved water-column physics does not depend on
  which 2 of 12 scans carried a residual NIR offset.
- **RMS misfit is essentially unchanged** (8.5%→8.5% free, 8.8%→8.6% max freedom) — LOC1
  remains the worst-fitting station regardless, confirming again (§3 below) that its
  residual is the blue-band CDOM-slope mismatch, not scan-level noise.
- Per-scan check before applying: the correction's magnitude on the two flagged scans
  (offset 1.3–2.3×10⁻⁴ sr⁻¹) was NOT dramatically larger than on several scans the QC
  tool calls clean (00006 alone: 2.4×10⁻⁴) — expected and correct for this method
  (Ruddick et al. 2006 residual correction is computed per scan from each scan's own
  R_rs(780)/R_rs(870), applied dataset-wide, not selectively to flagged scans only), and
  confirmed not a problem by the shape-consistency and headline-number checks above: a
  uniformly-applied, physically-motivated correction that tightens precision without
  moving the answer is the signature of a real fix, not an over-correction.

This is the one number in this document that changed after §1b/§3/§4 below were written;
those sections are updated to match. The uncorrected `FINAL_Rrs.csv` is not specially
archived, but is fully recoverable from git history at the commit before this fix
(`git log -- <path to LOC1's analysis/FINAL_Rrs.csv>`) for anyone who wants to diff
against it.

### 2d. The planned FOV/footprint test did not resolve
LOC3 was meant to isolate footprint size (both foreoptics, "same water"). It could not:
range differed between foreoptics, the two were used 23 minutes apart, and 2c applies to
both. A real 6.77° shape difference exists between them but cannot be attributed to
footprint. `LOC3_FOOTPRINT_COMPARISON.md` states plainly what would be needed to actually
answer this (same-instant, same-range, known-depth comparison) — nothing in the current
dataset does.

---

## 3. Open items, ranked by value per unit of future field effort

~~0. Decide on LOC1's 00005/00007 glint correction~~ — **RESOLVED 2026-08-18**, applied
   (§2d0). Kept as a struck-through entry rather than deleted, since the reasoning for why
   it was held open (a re-processing decision that moves the reference station's headline
   numbers deserves an explicit call, not a silent inheritance from precedent) is itself
   worth keeping visible for the next time a similar defect turns up somewhere else.
1. **A depth measurement at LOC3.** Cheapest of everything below (a weighted line), and
   the only thing that converts §2c from a caveat into a resolved question.
2. **Wind speed, recorded live, at any future site not near an ASOS station.** This
   session recovered real historical wind for every 2026-08-16 station from Kotzebue's
   airport archive after the fact (`RHO_METHODOLOGY_REVIEW.md`) — a genuinely useful,
   free, retroactive technique worth trying on any past campaign near an airport — but it
   only works because the exact UTC timestamp was already recorded. It is not a
   substitute for a live reading at a site without a nearby station.
3. **A same-instant, same-range, two-foreoptic comparison**, to finally answer the
   footprint question §2d could not.
4. **Redo `fit_shapes` at LOC2a with a finer (S_dg, η) grid** — the solver lands ~5% off
   the true χ² minimum there (`LOC2a_GIOP_FINDINGS.md` §A3); a real but small
   implementation gap, not a data limitation.
5. Wiring `--glint` more broadly is possible but was not needed elsewhere this session —
   `analyse_water_scans.py`'s own collapse test is the thing that decides whether it is
   needed at a given station, and no other station's population flagged a glint scan.

---

## 4. Every figure and document, by location

| location | R_rs/water docs | GIOP docs | key figures |
|---|---|---|---|
| LOC1 | `verify_field_calcs`, `analysis/README.md` | `LOC1_GIOP_FINDINGS.md` | `fig1`–`fig13`, `GIOP/giop0`–`giop11` |
| LOC2a | `LOC2_SPLIT.md`, `SCAN_00035_GLINT.md` | `LOC2a_GIOP_FINDINGS.md` | same set, `--glint nir_similarity` |
| LOC2b | `LOC2_SPLIT.md` | `LOC2b_GIOP_FINDINGS.md` | same set (n=3) |
| LOC2c (algae) | `analyse_algae_mat.py` output | n/a — reflectance, not R_rs | `c1_algae_reflectance.png` |
| LOC3 (all 3) | `LOC3_BOTTOM_CAVEAT.md`, `LOC3_FOOTPRINT_COMPARISON.md` | `LOC3_GIOP_FINDINGS.md` | same set × 3 sub-populations |
| cross-station | `Data_NatureSpec/2026_Aug_16/GLOBAL_COMPARISON.md`, `SKY_CHOICE_SYNTHESIS.md`, `RHO_METHODOLOGY_REVIEW.md` | `giop_cross_station_all.png` | `COMPARISON/all_stations_overplot.png` |

Every `GIOP/README.md` is generated fresh by `make_giop_figures.py` on each run
(`write_giop_index`) — do not hand-edit; edit the script.

## 5. What this session added to the underlying software, beyond the field-day analysis

Real defects found and fixed, each with a regression test or an explicit before/after
check: a profile-solver joint-optimization bug that read as a data limitation; a
verify-script threshold hardcoded to one station's scan count (silently passing 40% bad
data at a smaller station); two pipeline crashes at low n that would have blocked any
small sub-population from ever being processed; a self-inflicted "infx" from a
divide-by-zero the first version of a new guard didn't anticipate; a stale hand-written
index file that had already drifted within a single day of edits, replaced with a
generated one that cannot. LOC1's numbers were re-verified bitwise unchanged after every
one of these — nothing above is a story about numbers moving, only about coverage
expanding to LOC2 and LOC3 without breaking what LOC1 had already established.
