# How to get a better ρ for LOC1/LOC2 (and what LOC3 tells us about doing this at all)

ρ carries 20–100% of the R_rs error budget (`THEORY.pdf` p3) and is the largest single
term in this whole pipeline. What's already implemented, what was actually missing, and
what was recoverable after the fact.

## 1. What's already done, and its stated limit

`fieldrrs.rrs.rho_at_angle` scales Mobley's (1999) ρ = 0.028 (40°, 135° azimuth, wind
<5 m/s, clear sky) by the ratio of Fresnel reflectance at the achieved angle to Fresnel
reflectance at 40°. Measured to work: on LOC1, correcting per-scan angle removed a
p=0.006 tilt-vs-R_rs trend and cut scatter 11.1%→8.6% (`view_zenith_from_tilt`
docstring). This captures the **angular** trend only. The function's own docstring
already states the limit precisely: *"it captures the angular trend only, not the wind
or sky-state dependence... Beyond roughly ±15°, or above ~5 m/s wind, the real Mobley
(2015) tables are needed and those are not redistributable."*

`rho_advice(wind_ms, sky)` exists to give that guidance operationally, but it is
advisory only — never wired into the processing pipeline — and with no wind speed
recorded at any station this season, every run defaults to *"No wind speed recorded...
record the wind, it is the largest error in this measurement."*

**So the honest starting point is: nothing was wrong, but the single most-flagged gap in
this whole project (`NEXT_CAMPAIGN.md` item 2) was never actually closed.**

## 2. What was recoverable: real wind speed, after the fact, from public records

Kotzebue has an ASOS station at the airport (PAOT/OTZ), whose 5-minute observations are
archived and were retrieved for 2026-08-16 via the Iowa Environmental Mesonet's public
ASOS request API. This is the first time in this project that a genuinely independent
external measurement — not derived from the NaturaSpec data at all — has been brought in
to check a field assumption.

| station | UTC window | wind speed | Cox & Munk slope variance σ² |
|---|---|---|---|
| **LOC1** | 20:24–20:41 | mean 4.73 m/s, range 3.55–5.32 | 0.021–0.030 |
| **LOC2** | 21:05–21:18 | mean 4.54 m/s, range 4.14–4.73 | 0.024–0.027 |
| **LOC3-FLENS8** | 21:44–21:47 | 5.32 m/s | 0.030 |
| **LOC3-FIBR15** | 22:07–22:10 | mean 5.62 m/s, range 5.32–5.92 | 0.030–0.033 |

(Cox & Munk 1954: σ² = 0.003 + 0.00512·W, the mean-square wave-facet slope a wind-driven
sea surface presents; 0.003 is the flat-sea/zero-wind floor.)

**Every station in this dataset sits at or above the 5 m/s edge of ρ = 0.028's stated
validity, not comfortably below it as had been assumed by default.** LOC2 is the
closest to safely inside (4.14–4.73 m/s). LOC1 straddles the line — its own range spans
3.55 to 5.32 m/s across the 17-minute station, meaning the wind was *rising through* the
threshold during the very station that fixed the angle correction. LOC3-FIBR15 is
clearly above it at every sample.

## 3. What this means for the sign and (roughly) the size of the bias

Wind above the Mobley reference condition means MORE wave slope variance, which means
the **true** ρ is larger than 0.028, not smaller — a calmer sea reflects sky more like a
mirror at the nominal angle, a rougher sea samples a wider range of sky radiances and
(under a clear, non-uniform sky) that generally raises the effective reflectance
(`rho_advice`'s own physical reasoning, §the docstring in `rrs.py`).

Using ρ = 0.028 when the true value is higher means **under-subtracting** ρ·L_sky, so
`L_w = L_t − ρ·L_sky` runs systematically too high, and — because sky radiance is
strongly blue-weighted (Rayleigh scattering) — **the residual lands mostly in the blue
end of R_rs.** This is not a new claim invented for this document: it is exactly the
mechanism already isolated empirically this session. In the LOC2b (00027/28/29)
investigation, a candidate "sky reflection" term (proportional to L_sky/E_d, i.e. exactly
what an under-corrected ρ leaves behind) explained **51%** of an anomalous blue-elevated
residual, more than any other single candidate tried (`LOC2_SPLIT.md`).

**We cannot recover the exact numeric correction.** That needs the actual Mobley (2015)
wind- and geometry-dependent ρ table, which the package's own documentation states is not
redistributable, and this project has no local copy. What can honestly be said:

- The **direction** of the residual bias is known: R_rs is likely biased slightly HIGH,
  concentrated in the blue, at every station, worst at LOC3-FIBR15.
- The **order of magnitude** is bounded by the Cox & Munk slope-variance ratio: σ² at
  the measured winds is 7–11× the flat-sea floor (0.003), so the facet-averaging effect
  Mobley's table captures is not a small correction at these wind speeds — but slope
  variance is not itself ρ, and translating one into the other needs the actual
  radiative-transfer table.

## 4. What to actually do

**Retroactively, for LOC1/LOC2/LOC3 as recorded**: nothing further can be done inside
this package without the Mobley table. The angle correction already applied is real and
already measured to help; it is not wrong, it is simply not the whole ρ story. State the
wind-speed finding as an explicit caveat alongside every R_rs number from this field day,
the same way the shape/amplitude split and the shallow-water caveat are already stated.

**For the next campaign, in order of value per unit effort:**

1. **Record wind speed in the field.** Already `NEXT_CAMPAIGN.md`'s #2 item; this
   session's finding makes it more urgent, not less — every station this season needed
   it and none had it.
2. **Even without an anemometer, record the exact UTC timestamp precisely enough to
   repeat this trick.** This session shows historical ASOS data can be recovered *after
   the fact*, for free, for any station within reach of an airport METAR/ASOS site —
   which is most coastal Alaska field sites. Recording only "the wind was calm" in a
   field notebook throws away the ability to do this; recording the precise time (already
   done, since the instrument timestamps every scan) is what makes it possible. **This is
   a free, retroactive partial fix for past data too** — if a past campaign's site is near
   any ASOS/METAR station, the same recovery is worth trying before assuming the
   information is lost.
3. **If reliable field internet or a pre-downloaded Mobley (2015) table become
   available**, wire a real wind-dependent ρ into `fieldrrs.rrs`, gated behind a wind
   speed argument the way `rho_advice` already anticipates — the interface exists, only
   the table is missing.
4. **Photograph the water surface deliberately for sea state**, not only the target,
   framed to catch the horizon or a reference object — `NEXT_CAMPAIGN.md` item 8 already
   recommended this for a different reason (corroborating wind), and it is now the
   cheapest available substitute when ASOS recovery isn't possible (e.g. a genuinely
   remote site).

## References

- Mobley, C.D. (1999), *Estimation of the remote-sensing reflectance from above-surface
  measurements*, Applied Optics 38(36), 7442–7455, doi:10.1364/AO.38.007442
- Cox, C., Munk, W. (1954), *Statistics of the sea surface derived from sun glitter*,
  Journal of Marine Research 13(2), 198–227 (no DOI assigned; pre-DOI publication)
- Wind data: Iowa Environmental Mesonet ASOS archive, station PAOT (Kotzebue Ralph Wien
  Memorial Airport), retrieved for 2026-08-16 19:00–22:20 UTC,
  `mesonet.agron.iastate.edu/cgi-bin/request/asos.py`
