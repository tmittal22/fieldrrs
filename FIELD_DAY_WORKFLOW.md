# Processing a field day, end to end

Every script, in order, with the exact command, what it writes, and what to look at.
Worked on the 2026-08-16 Kotzebue day (22 scans at LOC1, 23 at LOC2, 32 at LOC3), the
only field day this pipeline has processed so far — every path below is that day's, swap
in a new date's folder name and everything else is unchanged.

`fieldrrs/` itself is standard-library only. **The scripts in this directory are not** —
they need numpy, scipy, matplotlib, plotly, and (for GIOP) the sibling
[`giop-workbench`](https://github.com/tmittal22/giop-workbench) checkout.

```bash
source ~/miniforge3/etc/profile.d/conda.sh && conda activate claude-science-env
cd fieldrrs
```

---

## 0. Look before you process

```bash
python survey_field_data.py Data_NatureSpec/2026_Aug_16
```

Prints one line per `.sed`: time, lat/lon, tilt, range, foreoptic, integration, and the
**inferred role** (panel / sky / water / land). Roles come from spectral shape, not from
filenames or tilt — see §1. Read this before anything else; it is how the two land
targets (00014 algae-on-concrete, 00015 bare concrete) were caught.

> **Why not tilt?** The instrument records an unlabelled tilt magnitude and both sky and
> water scans were taken at 36–50°. Tilt does not separate them. Shape does: sky is blue,
> water peaks in the green, the panel is flat and bright.

## 1. Reorganise by location

```bash
python organize_by_location.py Data_NatureSpec/2026_Aug_16
```

Clusters scans at a **60 m** threshold on lat/lon and writes
`by_location/LOC<n>_<lat>N_<lon>W/<FOREOPTIC>_FOV<nn>/`, copying the `.sed` and its
`.jpg`. The key is `(location, foreoptic)`, not location alone, because a 15° and an 8°
foreoptic at the same water are **not** interchangeable and must not be pooled by
accident. LOC3 has both, which is what makes the FOV comparison possible there.

Writes `by_location/INDEX.md`.

## 2. Map of where the stations are

```bash
python make_location_map.py Data_NatureSpec/2026_Aug_16
```

Esri World Imagery basemap, Web-Mercator tiles mosaicked by `basemap.py`. Tiles cache to
`.tilecache/`, which is **gitignored deliberately** — the imagery is Esri's, not ours.

## 3. Verify the arithmetic BEFORE inverting anything

```bash
python verify_field_calcs.py Data_NatureSpec/2026_Aug_16/by_location/LOC1_*/FLENS8_FOV08
```

Six independent checks, all of which must pass:

| check | what it proves | LOC1 result |
|---|---|---|
| algebra | our R_rs equals the definition evaluated by hand | bitwise |
| vs DARWin | our reflectance matches the instrument's own `Reflect.` column | 0.219 % median, mean signed residual 0.030 % against 0.435 % scatter, i.e. **no bias** |
| transmittance | E_d implies a physical atmosphere | T = 0.745–0.919 |
| round trip | synthesise from a known R_rs, read back | 1.8e-16 |
| ρ sensitivity | how much the answer rides on the one assumed constant | 32 % |
| conservation | L_w ≤ L_t, no negative radiances | pass |

**Do not skip this.** The GIOP stage cannot tell a processing blunder from unusual water.
The pass/fail thresholds scale with the station's own scan count (`verify_field_calcs.py`
V2/V6), not a hardcoded number, so this is safe to run on a station with 2 scans or 30.

## 4. The per-location analysis — 13–14 figures

```bash
python analyse_location.py Data_NatureSpec/2026_Aug_16/by_location/LOC1_*/FLENS8_FOV08
```

Writes `analysis/` with `REPORT.txt`, `FINAL_Rrs.csv` and figures 1–13 (fig7,
land-targets, only appears if the station has scans classified `land` — LOC1 has 2,
LOC2/LOC3 have none). Each figure and what it settles is indexed in `analysis/README.md`,
written by the same run — never hand-edit that index, it goes stale within a day (found
this session: an early hand-written `GIOP/README.md` referenced a since-renamed figure
and quoted withdrawn numbers within the same session it was written; see §6).

Pairing is **angle-matched and location-locked**: each water scan is paired to the sky
scan whose view geometry mirrors it, drawn only from the same location *and* the same
panel-reference block. `match_by_angle(..., respect_blocks=True)`; `assert_same_dataset()`
raises rather than silently pairing across stations.

At n<4 water scans, `fig_sensitivity` (angle/range confound control) and at n<3,
`fig_variability` (PC1 shape/amplitude split) write an explicit "insufficient data"
placeholder instead of crashing or silently omitting the figure — this matters for small
split-off sub-populations like LOC3's murky pair (n=2).

## 5. Is this really one water population?

```bash
python analyse_water_scans.py Data_NatureSpec/2026_Aug_16/by_location/LOC1_*/FLENS8_FOV08
```

Writes `analysis/water_scans/`: per-scan shape-angle clustering (is every water scan
close to the group, or is one/some a different water body?), a NIR-similarity
glint-collapse test per deviant scan (Ruddick et al. 2006 — does a glint correction pull
it back into the group, or is it genuinely different water?), a contact sheet of the
photos, and a verdict table (`clean` / `glint (correctable)` / `deviant, NOT glint`).

**Run this at every station, not just ones you already suspect are mixed** — it is what
caught LOC2's disturbed-water sub-population (§below) and LOC1's own two
glint-affected scans (00005, 00007, both `glint (correctable)`, found only when this
tool was finally run there — see the LOC1 row of the status table). Needs ≥3 water
scans; below that (LOC3's murky pair, n=2) it writes an explanation to
`water_scans/REPORT.txt` instead of silently producing nothing, since clustering needs
something to cluster against.

**If it finds a real second population**, split manually: make a new
`by_location/LOC<n><letter>_.../<FOREOPTIC>_FOV<nn>/` folder and symlink (not copy) the
shared sky/panel scans and the deviant water scans into it, so there is exactly one copy
of every raw `.sed`/`.jpg` on disk (see LOC3-FIBR15's murky-pair split for the pattern).
Then re-run steps 3–8 on the new sub-station independently. `LOC2_SPLIT.md` and
`LOC3_BOTTOM_CAVEAT.md` are the two worked examples.

## 6. Interactive version

```bash
python make_interactive_report.py Data_NatureSpec/2026_Aug_16/by_location/LOC1_*/FLENS8_FOV08
```

`analysis/REPORT.html`, 8 plotly panels, zoomable, self-contained (no CDN), openable from
Dropbox on any machine.

## 7. GIOP inversion

```bash
python make_giop_figures.py Data_NatureSpec/2026_Aug_16/by_location/LOC1_*/FLENS8_FOV08
```

Writes `analysis/GIOP/` — 12 figures (`giop0`–`giop11`), `giop_per_scan.csv` (every
individual angle-matched fit) and `giop_assumption_arms.csv` (every arm of the assumption
sweep, ranked by χ², with its weight and whether it is admissible).
`analysis/GIOP/README.md` is **auto-generated fresh on every run**
(`make_giop_figures.write_giop_index`) — it pulls the live M_φ/a_dg/b_bp/χ²_ν/RMS numbers
straight from the run and searches upward from the output folder for a
`<LOC>_GIOP_FINDINGS.md`, so it can never point at a stale number the way a hand-written
index can. The physics and general caveats are in `THEORY_GIOP_NOTE.md`; the station-
specific verdict, written once per location by hand (not auto-generated, since it is
interpretation, not a number dump) is `LOC1_GIOP_FINDINGS.md` /
`LOC2a_GIOP_FINDINGS.md` / `LOC2b_GIOP_FINDINGS.md` / `LOC3_GIOP_FINDINGS.md`.

Runs on the **amplitude-normalised mean** from step 4 and on each angle-matched pair
separately, hyperspectral over 400–700 nm — GIOP's own Bricaud a*_φ table does not extend
past 700 nm, so it never sees anything redder than that (relevant at LOC3, §below).

> ⚠ Read the relevant `LOC*_GIOP_FINDINGS.md` before quoting any concentration.
> **Never quote M_φ / chlorophyll** — it moves 5–7× across equally-valid fit
> configurations at every station (`THEORY_GIOP_NOTE.md` §5, a Case-1 model on Case-2
> water). **Never compare χ²_ν across stations** — it scales with each station's own
> measured σ, not fit quality; RMS misfit is the comparable number.

## 8. Every individual spectrum, and a site-level summary

```bash
python make_all_spectra_figs.py Data_NatureSpec/2026_Aug_16/by_location/LOC1_*/FLENS8_FOV08   # one sub-case: fig14
python make_all_spectra_figs.py --site LOC2                                                    # a site: all sub-cases together
python make_all_spectra_figs.py --site LOC3
```

The first form writes `analysis/fig14_all_spectra.png` — every individual water scan on
disk for that sub-case, full 400–900 nm range, real (non-normalised) amplitude, with any
scan `FINAL_Rrs.csv` excluded shown dotted grey and labelled. This is deliberately
independent of `--exclude-water`: it shows everything, including what was cut and why
(read back from the station's own `REPORT.txt`), rather than only the clean survivors.

The `--site` form (LOC2 or LOC3 only — the two sites with more than one sub-case) writes
`by_location/COMPARISON/<SITE>_site_summary.png` (one panel per sub-case, every scan) and
`<SITE>_site_means_overlay.png` (just the sub-case means, one axis) — the site paths are
hardcoded in the `SITES` dict at the top of the script; add a new site there rather than
generalising the CLI, since each site's sub-case list is a one-off decision anyway.

## 9. Tests

```bash
python tests/test_fieldrrs.py                    # 111, standard library only, the package itself
python -m pytest tests/test_field_day.py -q       # 26, the field-day scripts
cd ../giop_python && python -m pytest -q          # 140 passed + 1 skipped, the inversion
```

The GIOP suite's one skip is `test_essd2023_raman` when `$ESSD_NC_DIR`'s HydroLight
`.nc` files aren't staged locally — not a field-day-specific test.

---

## The whole day in one block

```bash
D=Data_NatureSpec/2026_Aug_16
python survey_field_data.py $D
python organize_by_location.py $D
python make_location_map.py $D
for L in $D/by_location/LOC*/*_FOV*; do
    python verify_field_calcs.py      "$L"
    python analyse_location.py        "$L"
    python analyse_water_scans.py     "$L"
    python make_interactive_report.py "$L"
    python make_giop_figures.py       "$L"
    python make_all_spectra_figs.py   "$L"
done
python make_all_spectra_figs.py --site LOC2
python make_all_spectra_figs.py --site LOC3
```

This loop is written for a station laid out from the start — it does **not** perform a
split. If step 5 (or `--site`'s side-by-side figure) turns up a second population inside
one of the `LOC*_FOV*` folders the loop found, stop, do the manual split described in §5,
then re-run the loop (or just steps 3–9) on the new sub-station folders it creates.

## Starting a NEW field day

Nothing below this line assumes 2026-08-16 specifically — every command already takes the
day's folder as an argument. The 2026-08-16 walkthrough above is worked through with real
numbers so you can see what "correct" looks like; a new day follows the exact same 10
steps with a different `$D`.

1. **Copy the instrument's raw export** into `Data_NatureSpec/<new_day>/` — every `.sed`
   (+ its `.jpg`/`.RAW` siblings) flat in that one folder, plus `HeaderLog.csv`. This flat
   layout is what `survey_field_data.py`/`organize_by_location.py`/`make_location_map.py`
   (steps 0–2) read directly, non-recursively — they will not find scans nested any deeper.
2. **Run steps 0–9 exactly as in "the whole day in one block" above**, with
   `D=Data_NatureSpec/<new_day>`. If `--site` should show a cross-sub-case summary for this
   day too (only matters once a station has been split, §5), add an entry to the `SITES`
   dict at the top of `make_all_spectra_figs.py` — it is intentionally a short hardcoded
   dict, not auto-discovered, since deciding which sub-cases belong on one summary page is
   a one-off editorial call, not something to infer from folder names.
3. **Once `by_location/` is built and every station's `analysis/` has been checked**,
   the flat raw archive is fully redundant with it (`organize_by_location.py` copies every
   `.sed`/`.jpg` into `by_location/`, one real copy per scan) and can be deleted to save
   space — this IS what was done for 2026-08-16 (10 MB freed), verified first (`diff` every
   flat filename against `by_location/`, 0 missing) and confirmed after with a full
   pipeline re-run showing bitwise-identical output. **Two things this breaks, both
   recoverable but not for free:**
   - `survey_field_data.py`/`organize_by_location.py`/`make_location_map.py` (steps 0–2)
     can no longer be re-run for that day from within the repo, since their one input was
     the flat folder — regenerate anything from them (especially the site map,
     `fig1_map.png`) **before** deleting, not after. Git history still has the raw files if
     truly needed later (`git checkout <commit before the delete> -- Data_NatureSpec/<day>`).
   - Any test exercising real files for that day needs to read `by_location/` instead. See
     `tests/test_field_day.py`'s `REAL_BACKING_STORES`, written for 2026-08-16's specific
     layout (LOC1, LOC2's pre-split folder, LOC3's two foreoptics — the 4 real,
     non-symlinked copies, chosen to reconstitute all 60 original scans with none double-
     counted through a split's symlinks). A new day's equivalent test class would need its
     own list, built the same way: the real (non-symlink) `by_location/` folders for that
     day, not the split-off station folders that symlink into them.
   - If disk space isn't a concern, **keeping the flat archive alongside `by_location/` is
     also completely fine** — nothing breaks either way, and it keeps steps 0–2 re-runnable
     for that day indefinitely. Deleting it is an optional space/reproducibility trade, not
     a required step.

## Status, 2026-08-17

| location | scans | foreoptic | steps 0–6 | GIOP | notes |
|---|---|---|---|---|---|
| **LOC1** 66.89718 N 162.60290 W | 22 (12 water) | FLENS8 (8°) | **done** | **done** | one population; step 5 (run retroactively this session) flags 00005/00007 as glint-correctable — **not yet applied to `FINAL_Rrs.csv`, open item, see `Data_NatureSpec/2026_Aug_16/PAPER_READINESS.md`** |
| **LOC2a** (main) 66.89677 N 162.57953 W | 9 water | FLENS8 (8°) | **done** | **done** | 00035 glint-corrected via `--glint nir_similarity`, verified to collapse back into the group |
| **LOC2b** (disturbed) | 3 water | FLENS8 (8°) | **done** | **done** | real, not artefact, n=3 — report with the weaker-n caveat |
| **LOC2c** (algae mat) | 2 | FLENS8 (8°) | **n/a — reflectance, not R_rs** | n/a | one figure + REPORT.txt by design (`analyse_algae_mat.py`); not the 13-figure R_rs pipeline, see its own module docstring for why |
| **LOC3-FIBR15 (main)** 66.89235 N 162.59149 W | 3 water | FIBR15 (15°) | **done** | **done** | clean open water |
| **LOC3-FIBR15 (murky)** | 2 water | FIBR15 (15°) | **done** | **done** | 00058/00059, high-sediment pair, split out same pattern as LOC2b; step 5 structurally skips at n=2 (writes why) |
| **LOC3-FLENS8** | 4 water | FLENS8 (8°) | **done** | **done** | see `LOC3_BOTTOM_CAVEAT.md` before quoting anything — every LOC3 sub-station's composition numbers are conditional on optical depth, which is unmeasured |

`by_location/LOC2_66.89677N_162.57953W/` (the pre-split folder) is **NOT a duplicate to
delete** — every raw `.sed`/`.jpg` under LOC2a/2b/2c is a **symlink back into it**
(`organize_by_location.py` copied the originals there once; the split added symlinks
rather than a second copy), so it is the one real backing store for all three. Its own
`analysis/water_scans/REPORT.txt` is also cited by name in `LOC2_SPLIT.md` as the
original evidence that found the heterogeneity in the first place — keep it. Read
results from 2a/2b/2c, never from here, but do not remove this folder.

**LOC2 is not one water population** — `analyse_water_scans.py` found 3 of its 12 water
scans (00027–29) are physically different water (disturbed sediment, settling over
~2 min at station start, not an AC/geometry artefact), and 2 more scans classified
"land" are a floating algae mat, not open water or dry ground. Split into
`LOC2a_.../` (9 scans, main open water), `LOC2b_.../` (3 scans, the disturbed water,
its own valid station), `LOC2c_.../` (2 scans, algae reflectance) — full reasoning and
evidence in `Data_NatureSpec/2026_Aug_16/by_location/LOC2_SPLIT.md`, one water scan's
individual glint disposition in `LOC2a_.../analysis/water_scans/SCAN_00035_GLINT.md`,
and the LOC1/2a/2b/2c comparison in
`Data_NatureSpec/2026_Aug_16/by_location/COMPARISON/loc1_loc2abc_overplot.png`.

**LOC3-FIBR15 is not one water population either** — the same tool found 00058/00059
visibly murkier (elevated GIOP `b_bp`, a failed glint-collapse, real spectral-angle
separation from the other 3 FIBR15 water scans) and split them into their own
`FIBR15_FOV15_murky/` station, symlinked back to the shared raw scans rather than copied.

LOC3 is also the controlled FOV comparison — same water, two foreoptics — and was meant
to be the one place in the dataset where the footprint question could be answered
properly. It could not: range differed between foreoptics and they were used 23 minutes
apart, both confounds independent of FOV. At LOC1 range took just 3 values and correlates
with time at r = 0.92, so footprint and time are confounded there too, for a different
reason. See LOC3's own `LOC3_FOOTPRINT_COMPARISON.md` (next to `LOC3_BOTTOM_CAVEAT.md`,
in the LOC3 site folder — one file covering all three LOC3 sub-stations, not repeated per
sub-station) for what a real test would need.

## Where the cross-station and synthesis results live

Everything above is per-station. Once every station in a field day is done, these pull
it together, living next to that day's own data (not at the repo root, since they are
specific to this field day, not to the pipeline) — all cite live numbers back to a
specific run:

| doc | question it answers |
|---|---|
| `Data_NatureSpec/<day>/GLOBAL_COMPARISON.md` | R_rs and GIOP composition across every station, one table |
| `Data_NatureSpec/<day>/RHO_METHODOLOGY_REVIEW.md` | how to get a better ρ than the flat Mobley (1999) default — recovers real ASOS wind speed after the fact |
| `Data_NatureSpec/<day>/SKY_CHOICE_SYNTHESIS.md` | does the sky-scan choice matter (answer: no, at every station, for a stated geometric reason) |
| `Data_NatureSpec/<day>/by_location/LOC3_.../LOC3_FOOTPRINT_COMPARISON.md` | why the FOV/footprint question could not be answered cleanly this field day |
| `Data_NatureSpec/<day>/PAPER_READINESS.md` | what is quotable now, what is conditional, what is still open, ranked by value per unit of future field effort |
| `<station>/analysis/GIOP/LOC*_GIOP_FINDINGS.md` | the station-specific GIOP verdict — co-located with that station's own GIOP output, or (LOC3) at the shared site folder when one file covers several sub-stations |
