# No pooled analysis here — see `../../LOC2_SPLIT.md`

This folder holds `water_scans/` only: the QC run (`analyse_water_scans.py`) that
established the 12 water scans here are **not one population** — 3 are physically
different water (disturbed sediment, settling over ~2 minutes at station start) and 2 of
the "land"-classified scans are a floating algae mat, not open water.

A pooled `FINAL_Rrs.csv` averaging all 12 as one spectrum was generated during that
investigation and has been **removed** — averaging across populations is not a valid
water product, and leaving it in this folder would risk it being quoted by mistake.

The valid products are in the split stations, symlinked from this folder's own raw
`.sed`/`.jpg` files (nothing here was duplicated):

- `../../LOC2a_66.89677N_162.57953W/FLENS8_FOV08/analysis/` — main open water, n=9
- `../../LOC2b_66.89677N_162.57953W/FLENS8_FOV08/analysis/` — disturbed water, n=3
- `../../LOC2c_66.89677N_162.57953W/FLENS8_FOV08/analysis/` — algae mat reflectance, n=2

Full reasoning: `../../LOC2_SPLIT.md`.
