# What to change next time

Derived from the LOC1 analysis, ordered by how much each would improve the result per
unit of field effort. Every number cited is measured in this repository, not general
advice.

---

## Tier 1 — changes that unlock something currently impossible

### 1. Take one water sample per station. This is the single biggest gain.

Everything that failed in `LOC1_GIOP_FINDINGS.md` failed for the same reason: **nothing
independent constrains the retrieval.** With three free amplitudes fitted to six smooth
bands, a 1.7 % input uncertainty becomes 47 % on chlorophyll, and an *assumed* S_dg moves
it from 0 to 297.

One filtered sample per station fixes this permanently:

| measurement | what it pins | what it kills |
|---|---|---|
| **SPM** (dry mass on a filter, gravimetric) | the amplitude term directly | the entire b_bp scale ambiguity |
| **Chlorophyll** (filter, extraction) | M_φ | the Case-1/Case-2 contradiction, factor 2.4 |
| **CDOM absorption** (0.2 µm filtrate, spectrophotometer) | **S_dg and a_g(443)** | the 0→297 swing |

The CDOM one is the highest value per unit effort: it converts the single most damaging
assumed constant into a measured quantity, and needs only a syringe filter, a vial, and a
bench spectrophotometer later.

*Cost:* about 5 minutes per station.

> **Updated 2026-08-17, and the reason has changed.** The original argument was that
> S_dg is *unconstrained* by R_rs, so a sample is the only way to pin it. That is now
> known to be wrong: with the shape solver fixed, **S_dg is the best-constrained
> parameter in the whole model** — a sharp interior χ² minimum at 0.011–0.013 nm⁻¹,
> χ²_ν 74 → 17 against the assumed 0.018 (`LOC1_GIOP_FINDINGS.md` §2d). The sample is
> still the top item, but now as a **falsification test**: a measured a_g(λ) on filtrate
> gives S_dg directly, and if it disagrees with the retrieved 0.0115 then the retrieval
> is absorbing something else into the CDOM slope — which, given that η simultaneously
> rails at an unphysical −1, is a live possibility. One vial converts a
> self-consistency argument into an external check.

### 2. Record the wind speed

ρ carries 20–100 % of the R_rs error budget (THEORY.pdf p3) and is the largest single
term by a wide margin — larger than everything else combined. ρ = 0.028 is only defensible
below ~5 m s⁻¹, and **nothing in the current files records whether that held.** Without it,
the data cannot be reprocessed with a better ρ later, ever.

A handheld anemometer, or Beaufort off the water using the field card's scale. Record
sky state too (clear / uniform overcast / broken).

*Cost:* 30 seconds per station.

**If you forget anyway, it may not be lost.** This session recovered real historical wind
for every 2026-08-16 station after the fact from Kotzebue airport's public ASOS archive
(`RHO_METHODOLOGY_REVIEW.md`) — free, no equipment, works for any site within reach of an
airport METAR/ASOS station (most coastal sites). It only works because the exact UTC
timestamp was already recorded (the instrument does this automatically). This is a
fallback for a missed reading, not a substitute for one — it can't tell you gusts or
local shielding an anemometer would catch.

---

## Tier 2 — changes that would have made this dataset cleaner

### 3. Hold the view angle steady, or at least record which datum you used

Measured at LOC1: over the 8° spread of achieved tilt, R_rs(443) moves **+27 %**
(r = +0.74, p = 0.006), and the trend survives controlling for time and range. That is
larger than the entire choice of sky scan.

Correcting per-scan with ρ(θ_v) removes it (p → 0.244) and cuts scatter 11.1 % → 8.6 %.
But it depended on determining the tilt datum *from the data*, because the instrument
reports an unlabelled magnitude and the two readings differ by 90°. **A note saying
"tilt is from nadir" would have removed that inference.**

Practical: aim for 40° and accept ±3° rather than ±8°, and write the datum down once.

### 4. Spread the SKY scans over a wider range of angles, or stop taking so many

At LOC1 the 8 sky scans span only ~5°. Consequence: matching a sky to each water by
mirrored geometry is **not better than picking one at random** (72 % of random draws are
as tight). The sky choice contributes 9 % against the water's 33 %.

So either
- take **fewer** sky scans (3 bracketing the sequence is enough) and spend the time on
  water replicates or a sample, **or**
- deliberately vary the sky angle to span the water angles, making the mirror-matching
  meaningful.

Currently 8 sky scans buy very little.

### 5. Keep the range fixed, or vary it deliberately

Range took only 3 values (0.98/1.17/1.23 m), and it correlates with time at r = 0.92, so
footprint and time are confounded and neither can be tested properly. Either hold the
range constant, or vary it over a real span (say 1–4 m) *within* a fixed time window, so
that footprint becomes a testable factor rather than a nuisance.

---

## Tier 3 — worth doing, lower payoff

### 6. Bracket each station with a panel scan

The panel was re-referenced once mid-station at LOC1 and the two agreed to 1.6 %, which
is reassuring but is a sample of two. A panel before *and* after each station bounds the
illumination drift over the actual measurement window rather than inferring it.

### 7. Repeat one station later in the day

Nothing in this dataset separates **spatial** patchiness from **temporal** change, because
each location was visited once. The LOC1 variability shows no time trend within 17
minutes (r = +0.03), but a return visit 2–3 hours later would establish whether the
factor-1.45 amplitude spread is a fixed property of the site or drifts with tide and light.

### 8. Photograph the water surface, not just the target

The photos are already valuable — they are what identified 00014/00015 as algae and
concrete. A shot showing surface state (glassy / rippled / whitecapping) would
independently corroborate the wind record.

---

## What NOT to change

- **The three-scan protocol.** It works. Closure is exact and the arithmetic agrees with
  the instrument's own reflectance column to 0.219 % median with no bias.
- **The panel.** Panel replicates hold to 0.6 %, which is 18× below the water variability
  and is why the variability can be attributed to water at all.
- **The solar timing.** Every station sat inside the preferred 30–60° window. At 66.9 °N
  in August the sun peaks at 36°, so that window is available for only part of the day and
  it was used correctly.
- **Number of water replicates.** 12 was ample: the shape converged to 1.7 %, and the
  residual spread is real water, not noise, so more replicates would measure the patchiness
  better but not the spectrum better.

---

## Tier 1b — added 2026-08-17: measure past 700 nm, or accept that the misfit stays

The best model in the GIOP family still misfits these spectra by **8.5 % RMS against a
1.9 % measurement uncertainty**, and the residual is one smooth curve with a +9σ lobe at
560–600 nm and a −20σ notch at 690 nm — *identical in all 12 scans*, so it is model
structure, not noise or scan variability. The instrument already records to 2500 nm; GIOP
stops at 700 nm because Bricaud's a*_φ table does.

So the fix is not more replicates. It is either

- a **turbid-water inversion** that extends past 700 nm (QAA-turbid, or titanspec's
  registry with a mineral component), which costs no field effort at all, or
- **an independent constraint on backscatter**, since η running to the −1 bound says the
  b_bp power law cannot make the red-rising backscatter the data want. A transmissometer
  or a backscatter sensor would settle it; failing that, an SPM sample plus a stated
  particle-size assumption bounds it.

This is now the largest single item, ahead of anything in Tier 2, because it is the term
that no amount of better field technique will reduce.

---

## The one-line version

**Bring a syringe filter, a vial and an anemometer.** The optics are already good enough
that the limiting factor is no longer the measurement — it is that nothing independent
constrains the inversion.
