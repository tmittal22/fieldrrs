# Processing a field day, end to end

Every script, in order, with the exact command, what it writes, and what to look at.
Worked on the 2026-08-16 Kotzebue day (22 scans at LOC1, 23 at LOC2, 32 at LOC3).

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

## 4. The per-location analysis — 12 figures

```bash
python analyse_location.py Data_NatureSpec/2026_Aug_16/by_location/LOC1_*/FLENS8_FOV08
```

Writes `analysis/` with `REPORT.txt`, `FINAL_Rrs.csv` and figures 1–12. Each figure and
what it settles is indexed in `analysis/README.md`, written by the same run.

Pairing is **angle-matched and location-locked**: each water scan is paired to the sky
scan whose view geometry mirrors it, drawn only from the same location *and* the same
panel-reference block. `match_by_angle(..., respect_blocks=True)`; `assert_same_dataset()`
raises rather than silently pairing across stations.

## 5. Interactive version

```bash
python make_interactive_report.py Data_NatureSpec/2026_Aug_16/by_location/LOC1_*/FLENS8_FOV08
```

`analysis/REPORT.html`, 8 plotly panels, zoomable, self-contained (no CDN), openable from
Dropbox on any machine.

## 6. GIOP inversion

```bash
python make_giop_figures.py Data_NatureSpec/2026_Aug_16/by_location/LOC1_*/FLENS8_FOV08
```

Writes `analysis/GIOP/` — 10 figures, `giop_per_scan.csv` (the 12 individual fits) and
`giop_assumption_arms.csv` (every arm of the assumption sweep, ranked by χ², with its
weight and whether it is admissible). Indexed in `analysis/GIOP/README.md`; the physics
and the caveats are in `THEORY_GIOP_NOTE.md`, the verdict in `LOC1_GIOP_FINDINGS.md`.

Runs on the **amplitude-normalised mean** from step 4 and on each of the 12 angle-matched
pairs separately, hyperspectral over 400–700 nm.

> ⚠ Read `LOC1_GIOP_FINDINGS.md` Correction 2 before quoting any concentration. The best
> arm of 24 has χ²_ν = 18 against a measured 1.9 % band uncertainty. Nothing in the GIOP
> family fits this water well.

## 7. Tests

```bash
python tests/test_fieldrrs.py     # 58, standard library only, the package itself
python -m pytest tests/test_field_day.py -q    # the field-day scripts
cd ../giop_python && python -m pytest -q       # 140, the inversion
```

---

## The whole day in one block

```bash
D=Data_NatureSpec/2026_Aug_16
python survey_field_data.py $D
python organize_by_location.py $D
python make_location_map.py $D
for L in $D/by_location/LOC*/*_FOV*; do
    python verify_field_calcs.py    "$L"
    python analyse_location.py      "$L"
    python make_interactive_report.py "$L"
    python make_giop_figures.py     "$L"
done
```

## Status, 2026-08-17

| location | scans | foreoptic | steps 0–5 | GIOP |
|---|---|---|---|---|
| **LOC1** 66.89718 N 162.60290 W | 22 | FLENS8 (8°) | **done** | **done** |
| **LOC2** 66.89677 N 162.57953 W | 23 | FLENS8 (8°) | **done, split into 2a/2b/2c** | LOC2a/2b pending |
| LOC3 | 32 | FIBR15 (15°) **and** FLENS8 (8°) | not run | not run |

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

LOC3 is the controlled FOV comparison — same water, two foreoptics — and is the only
place in the dataset where the footprint question can be answered properly. At LOC1 range
took just 3 values and correlates with time at r = 0.92, so footprint and time are
confounded there and neither can be separated.
