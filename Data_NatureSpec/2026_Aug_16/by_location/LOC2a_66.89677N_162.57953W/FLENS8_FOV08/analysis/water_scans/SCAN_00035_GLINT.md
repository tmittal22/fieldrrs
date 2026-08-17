# Scan 00035 — what was found, what was done, and what was NOT done

One water scan at LOC2a, flagged by `analyse_water_scans.py`. This file exists because
the decision about it is more involved than "kept" or "discarded", and that decision
should be inspectable on its own rather than buried in a table.

## 1. What it is

`NaturaSpecPlus_SN25494G1_00035.sed`, LOC2a (main open-water group, 8 clean scans +
this one), tilt_y = 42.7° (view zenith 42.7°, the shallowest in the group — the group
spans 43.2–45.5°), paired to sky scan 00042 by angle matching, range 0.97 m.

## 2. Why it was flagged

`analyse_water_scans.py` computes the pairwise spectral **shape** angle (amplitude
removed) between every water scan at the station. 00035 sits **3.16°** from the main
group's mean, while the other 8 scans sit only **1.08°** apart from each other — roughly
3× the group's own internal spread. See `w1_contact_sheet.png` and `w2_per_scan.png`.

## 3. Why it was called GLINT, not "different water" (unlike 00027/28/29)

Two independent lines of evidence, not one:

**(a) The photograph.** `NaturaSpecPlus_SN25494G1_00035.jpg` shows visible surface
ripple/sheen across the upper part of the frame that is not present in the neighbouring
scans (00030–00038). This is the direct, non-spectral evidence.

**(b) The collapse test.** Sun glint is an *additive* surface term
(`rho * L_sky`-like, but from wave facets rather than the diffuse reflectance the
package already subtracts), so a correction that targets exactly that kind of term must
bring a glint-contaminated scan back toward the group it belongs to. A scan that is
different water for some other reason has no reason to collapse under a NIR correction.
Testing both packaged glint corrections (`fieldrrs.rrs.rrs_three_scan`,
Ruddick et al. 2006):

| method | main group's own spread | 00035's distance from the group |
|---|---|---|
| none (as measured) | 1.93° | 3.16° |
| `nir_zero` | 1.81° | 1.98° |
| **`nir_similarity`** | **1.51°** | **1.53°** |

Under `nir_similarity`, 00035 lands **inside** the main group's own spread (1.53° vs
1.51°). That is the signature of glint, and it is why 00027/28/29 were NOT called glint
at LOC2 — the same test left them at 2.84° against a group spread of 1.51°, i.e.
unchanged in kind, and their photographs show no comparable ripple.

## 4. What was actually done with it — READ THIS PART CAREFULLY

**00035 was kept in the LOC2a population, UNCORRECTED, at its raw measured value.**

This is deliberate, and it is *not* the same thing as saying the scan is fine:

- `analyse_location.py`, the script that forms `FINAL_Rrs.csv`, calls
  `rrs_three_scan(..., "none")` at every call site (`fig_rrs`, `fig_final`,
  `fig_scaled_method`, `fig_final_mean`, `match_by_angle`-based figures). It has **no
  glint-correction option wired into the final-product path.** Passing
  `glint="nir_similarity"` through the pipeline would require touching five call sites
  and re-validating LOC1 against it, which was out of scope for this pass.
- The **amplitude-normalised (scaled) mean** (`THEORY_SCALED_MEAN.md`) that produces
  `FINAL_Rrs.csv` corrects only the *amplitude* mismatch between scans, iteratively. It
  does **not** correct a *shape* mismatch — and glint changes shape (it adds a spectrally
  flat-ish term on top of a green-peaked water spectrum, so it is not a pure scaling).
  So 00035's glint contribution passes through into the final shape uncertainty
  unmodified.
- Excluding it outright (the way 00027/28/29 were made their own station, LOC2b) was
  considered and rejected: unlike those three, 00035 is not a different physical water
  state — the photograph and the collapse test agree it is the SAME water seen through a
  transient glint. Throwing it away would be discarding a real measurement for a
  correctable reason, at a station that already has one fewer scan than LOC1.

## 5. Measured cost of leaving it uncorrected

With 00035 included (9 scans), `fig11_scaled_mean_method.png` /
`fig12_FINAL_mean_Rrs.png` report a shape uncertainty of **0.9 %** over 450–700 nm —
*tighter* than LOC1's 1.7 % despite LOC1 having 12 clean scans against LOC2a's 9. So at
the level the amplitude-normalised mean actually reports, 00035's glint residual is not
visibly inflating the final uncertainty. That is evidence the cost is small, not proof it
is zero — a shape distortion concentrated in a few bands can hide inside a broadband
uncertainty number.

## 6. If this needs to be tightened later

The correct fix is threading a `glint` option through `analyse_location.py`'s final-path
calls (the four listed in §4), defaulting to `"none"` so LOC1's numbers are provably
unchanged (as `--exclude-water` was required to be, verified bitwise), and re-running
LOC2a with `glint="nir_similarity"`. Until that exists, treat any LOC2a number that is
sensitive to the exact shape near 400–450 nm or 850–900 nm (where the collapse test shows
00035 deviates most) with one extra grain of caution.

## Bottom line

| | |
|---|---|
| flagged for | shape angle 3.16° vs main group's 1.08° internal spread |
| photographic evidence | visible ripple in `00035.jpg`, absent in neighbours |
| collapse test | PASSES — `nir_similarity` brings it to 1.53°, inside the 1.51° group spread |
| verdict | glint, correctable, same water as the rest of LOC2a |
| what was done | **kept, uncorrected** — the final-product pipeline has no glint-correction path wired in |
| measured cost | shape uncertainty 0.9 % (LOC1: 1.7 %), so not visibly inflated, but not independently verified band-by-band |
| what to do next | thread `glint="nir_similarity"` through `analyse_location.py`'s 4 final-path call sites and re-run |
