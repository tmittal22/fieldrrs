# Does footprint matter? LOC3's FIBR15-vs-FLENS8 comparison, and why it isn't clean

LOC3 is the only station in this dataset with both foreoptics (FIBR15, 15°; FLENS8, 8°),
which is why it was the planned test for whether footprint size systematically affects
R_rs — at LOC1, range correlated with time at r=0.92, so footprint and time were
confounded there and the question could not be asked at all.

## The result

| | FIBR15 (15°, n=3 clean) | FLENS8 (8°, n=4) |
|---|---|---|
| footprint diameter | **0.231 m** | **0.152 m** (1.52× smaller) |
| R_rs(443) | 0.00335 | 0.00308 (+8.5%) |
| R_rs(555) | 0.00970 | 0.00994 (−2.4%) |
| R_rs(700) | 0.00965 | 0.01122 (−14.0%) |
| shape angle between them (amplitude removed) | **6.77°** | |
| a_dg(443), GIOP free | 1.10 m⁻¹ | 1.22 m⁻¹ (+11%) |
| b_bp(443), GIOP free | 0.079 m⁻¹ | 0.091 m⁻¹ (+15%) |

**There is a real difference — 6.77° of spectral shape angle is large.** For comparison,
LOC1's own 12-scan internal replicate spread was 1.93°, LOC2a's 9-scan spread 1.08°.
This is 3.5–6× the scale of ordinary replicate noise, similar in magnitude to the
difference *between* LOC2a and LOC2b's genuinely different water populations.

## Why it cannot be attributed to footprint size

**This is not the controlled comparison it was meant to be**, for three compounding
reasons, none of them anticipated when LOC3 was planned:

1. **Range was not held equal.** FIBR15 was operated at 0.877 m, FLENS8 at 1.085 m — a
   deliberate or incidental choice, but it means the footprint difference (1.52×) is not
   FOV alone; it is FOV *and* range acting together, same confound as LOC1's, just
   smaller.
2. **The two foreoptics were not used simultaneously.** FLENS8 was scanned 21:44–21:47
   UTC; FIBR15 was scanned 22:07–22:10 UTC — **23 minutes later**. Anything that changed
   at the site in that window (tide, a passing patch of turbidity, sun angle moving from
   36.6° to 35.2°) is now inseparable from a footprint effect.
3. **`LOC3_BOTTOM_CAVEAT.md` applies to both.** If R_rs at this site is a variable
   mixture of water-column and bottom signal, a shape difference between two scans could
   reflect different depths or different bottom visibility at the two spots/times the
   sensor was aimed, which has nothing to do with footprint at all.

Any one of these three would be enough to block a clean footprint conclusion; having all
three at once means the honest answer is **not resolvable from this dataset**, not
"footprint doesn't matter" and not "footprint matters" — genuinely unknown.

## What would actually answer the question

A same-spot, same-instant (or immediately sequential, within a few seconds) pair of
scans with **both foreoptics at the identical range**, ideally over water known (by a
depth measurement — see `LOC3_BOTTOM_CAVEAT.md`) to be optically deep. That isolates FOV
as the only variable. Nothing in this dataset does that; it is the natural addition to
`NEXT_CAMPAIGN.md` alongside the wind-speed and depth recommendations, and cheaper than
either — it costs one extra scan per station, swapping foreoptics with the tripod
otherwise undisturbed.

## What this does NOT undermine

Every within-foreoptic result elsewhere in this project (LOC1, LOC2a/b/c, and LOC3's own
FIBR15-vs-FLENS8-internally-consistent numbers) is unaffected — this note is about
whether footprint *size* is a confound across a *foreoptic change*, not about whether any
single foreoptic's own measurements are trustworthy.
