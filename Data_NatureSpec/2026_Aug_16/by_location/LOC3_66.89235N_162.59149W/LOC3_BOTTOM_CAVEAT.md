# LOC3 — the bottom is visible, and this is a different regime from LOC1/LOC2

**Found while doing the per-scan photo QC pass, and it changes how every number from
this location must be read.** In every one of FLENS8's four water photographs
(00045–00048), the substrate — cobbles/rock, mottled dark-and-light texture — is
directly visible through the water column, sharpest in 00047, softer but still present
in the other three. FIBR15's water photos (00055–00059) do not show the same sharp
texture, but they are also close-up/blurred at the panel and cannot rule it out.
Nothing comparable appears in LOC1's or LOC2's water photographs at any station.

**Consequence: LOC3 may not be optically deep.** GIOP, QAA, and every semi-analytical
water-column model applied elsewhere in this project assume the water column alone
determines R_rs. If the sensor is also seeing light that reflected off the bottom and
came back up through the water column, R_rs at LOC3 is a **mixture of water-column and
bottom signal**, and no water-column-only inversion (GIOP included) can be trusted to
separate them. This is a structurally different problem from anything else in the
dataset — turbidity, glint, and disturbed sediment (LOC2) are all still water-column
effects; bottom reflectance is not.

## Update, 2026-08-18: confirmed independently, in the raw radiance, not just the photos

Overplotting every validated station's R_rs (`COMPARISON/all_stations_overplot.png`)
made the difference impossible to miss: **LOC1 and LOC2 collapse smoothly toward zero
past ~700 nm, the way deep, optically clean water must** (pure water absorption is
strong and rises steeply through the NIR). **LOC3's spectra do not — they show a
second, distinct peak near 805–810 nm** and retain substantial signal out to 900 nm.

Checked directly against the RAW target radiance (`rad_target`, before any R_rs
processing) so this cannot be an artefact of E_d, ρ, or the panel chain:

| | L(700) | L(810) | **L(810)/L(700)** |
|---|---|---|---|
| LOC1 | 0.00383 | 0.00087 | **0.226** |
| LOC3-FLENS8 | 0.00897 | 0.00373 | **0.416** |

LOC3 reflects nearly **twice as much at 810 nm relative to 700 nm** as LOC1 does, in the
sensor's own raw measurement. A sand/rock/algae substrate is far more NIR-reflective
than water; deep water alone cannot produce this. This is now the strongest single piece
of evidence for the bottom-reflectance hypothesis in the dataset — independent of, and
more decisive than, the photographic evidence above — and it applies to LOC3 broadly,
not only to the scans flagged for other reasons.

## What this does and does not explain

- **It is consistent with, and a more complete explanation than, the "different water
  population" findings below.** The two sub-clusters found in FLENS8's 4 scans
  (`{00045,00046}` vs `{00047,00048}`, 5–7° apart in shape, and a 700/675 nm ratio of
  ~1.19 vs ~1.55) line up with 00047 being the photo with the sharpest, most visible
  bottom — plausibly a real difference in local water depth or clarity across the ~1 m
  the sensor moved between scans, not a water-column composition change at all.
- **It does NOT explain FIBR15's 00058/00059** as cleanly. Those two scans were
  identified as a separate population by the same spectral-angle/glint-collapse test
  used at LOC2 (`analyse_water_scans.py`; not glint — does not collapse under NIR
  correction), and their photos don't show more or less bottom texture than the "clean"
  three. Kept as the documented exclusion, on the original evidence, not the bottom
  hypothesis.
- **It is not confirmed.** No depth measurement exists for this station (the field
  protocol does not record water depth), and "texture visible in a phone photo" is not a
  quantitative bathymetry. This is flagged as the most likely explanation for a real
  spectral pattern, not asserted as fact.

## What was done about it

**Nothing was silently corrected**, because there is no water-column-only correction for
bottom reflectance without a depth/albedo model this dataset cannot support (no depth
measurement, and a single R_rs spectrum cannot by itself separate a shallow bright bottom
from a genuinely different water column). Separating the two requires exactly the kind of
model this project does not have the data to run: Lee, Carder, Chen & Peacock (2001,
*Properties of the water column and bottom derived from AVIRIS data*, JGR-Oceans 106(C6),
11639–11651, doi:10.1029/2000JC000554) jointly retrieve depth, bottom albedo, and water
IOPs from a semi-analytical model — but that requires either a depth prior or enough
independent spectral information to break the depth/bottom-albedo/water-IOP degeneracy,
neither of which a single R_rs spectrum with no depth measurement provides. Instead:

- FLENS8's water scans are processed as one population (analysed as usual, all 4
  included), with this caveat attached rather than an unjustified split.
- Every LOC3 GIOP finding in this pass carries an explicit "shallow water, results
  conditional on optical depth" caveat, and none of the LOC3 GIOP composition numbers
  are used in the cross-station R_rs/GIOP comparison as if they were on equal footing
  with LOC1/LOC2's confirmed-deep-water numbers.
- **The single highest-value fix for a future visit to this site is a depth
  measurement** (a weighted line, or even a phone-camera stereo estimate) at each
  station — cheaper than almost anything else in `NEXT_CAMPAIGN.md` and the only thing
  that would let this question be closed rather than flagged.
