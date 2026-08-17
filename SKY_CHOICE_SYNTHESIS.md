# Does the sky choice matter? Answered once at LOC1, now confirmed at every station

`NEXT_CAMPAIGN.md` item 4 (written after LOC1 alone) said: *"the 8 sky scans span only
~5°... matching a sky to each water by mirrored geometry is not better than picking one
at random."* Six more stations later, across two more locations and two foreoptics, that
finding holds everywhere, for the same underlying reason.

## The numbers

| station | sky angular span | median sky-choice spread | verdict |
|---|---|---|---|
| LOC1 | 5.3° | 1.47% | does NOT matter |
| LOC2 (a+b combined) | 6.0° | 0.87–0.93% | does NOT matter (LOC2b: worst scans exceed its own tiny shape floor, see below) |
| LOC3-FIBR15 | 5.6° | 0.44% | does NOT matter |
| LOC3-FLENS8 | 2.9° | 0.36% | does NOT matter |

**Every single station shows the sky scans aimed within a 3–6° band.** Not a
coincidence of one operator on one day at one site — it happened identically at three
locations, two different foreoptics, and (per `RHO_METHODOLOGY_REVIEW.md`) across a wind
range that actually varied station to station. The aiming discipline was consistent even
when conditions were not.

## Why this is causal, not just correlated

Angle-matching can only outperform a random sky choice if the sky scans span enough
angular range that different picks would actually predict different ρ-corrected results.
With every station's sky spread under 6°, `rho_at_angle`'s own Fresnel-ratio correction
moves by only a percent or so across that whole range — smaller than the shape
uncertainty at every station. **There was never enough geometric diversity in the sky
reference for the matching procedure to have anything to bite on.** This is not a
statement about whether angle-matching is a sound idea (the physics is right, and it
demonstrably fixed a real angle-vs-R_rs trend on the WATER side at LOC1 — a 27% effect
over an achieved water-tilt spread that *did* vary by ~8°). It is a statement about this
particular field day's sky-aiming practice specifically.

## The one apparent exception, and why it isn't one

LOC2b is the only station where the sky choice "matters" by the printed verdict (worst
scans exceed the shape uncertainty). Its absolute sky-choice spread (0.93%) is
unremarkable — smaller than LOC1's (1.47%) and comparable to LOC2a's (0.87%) and LOC3's
(0.36–0.44%). It only crosses its own threshold because LOC2b's shape uncertainty is
exceptionally tight (0.2%, from 3 near-duplicate scans 0.53° apart — see
`LOC2b_GIOP_FINDINGS.md`). **This is the same "denominator effect" already documented for
χ²_ν comparisons across stations** (`THEORY_GIOP_NOTE.md` §4 point 4): a fixed numerator
compared against a station-specific, sample-size-dependent floor. Not a real difference in
how much the sky choice cost at LOC2b — the same denominator caution applies here too.

## What to actually change

1. **Stop taking 8–9 sky scans per station.** At every station this season they spanned
   a few degrees and contributed no information angle-matching could use. 3 sky scans
   bracketing the water-scan sequence, as `NEXT_CAMPAIGN.md` already recommended, is
   sufficient — the freed time is better spent on more water replicates or (per
   `RHO_METHODOLOGY_REVIEW.md` and `LOC3_BOTTOM_CAVEAT.md`) a wind reading or a depth
   measurement.
2. **If the value of angle-matching itself is ever to be tested properly**, the sky scan
   needs to be *deliberately* aimed across a wide range (say 20–40°) at least once per
   campaign — a dedicated calibration station, not mixed into ordinary water stations.
   Nothing in this dataset can currently distinguish "angle-matching doesn't help" from
   "this dataset never gave it the chance to."
3. **The water-side angle correction is a different, already-validated story and should
   not be cut.** It fixed a real, large effect (LOC1, 27% over the achieved water-tilt
   spread) precisely because the water tilt *did* vary enough, station to station, for
   the correction to matter. Only the sky side shows this "no information" pattern.
