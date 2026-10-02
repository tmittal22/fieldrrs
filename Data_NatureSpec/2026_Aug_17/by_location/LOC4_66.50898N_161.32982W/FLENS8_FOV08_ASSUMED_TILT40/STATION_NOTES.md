- **Tilt ASSUMED.** The tilt sensor logged `n/a` for all 7 scans (hand-held). The scans
  here are derived copies of `../FLENS8_FOV08/` with `Tilt (Y): +40.0°` (the protocol
  view angle) and a `Comment:` line saying so; nothing else in the files differs. ρ and the
  sky pairing therefore rest on that assumption, and with every tilt equal, pairing
  reduces to "same panel block".
- No range recorded, so no footprint.
- The camera on this hand-held setup is offset from the foreoptic: the photos (arm,
  boat, cooler) do not show the 8° field of view.
- Glint correction removed ~0.014 sr⁻¹, about 5x LOC5. Most likely ρ·L_sky is
  under-subtracted (an arm is in the sky photos and the angle is assumed). Use the shape;
  treat the absolute level as tentative. The uncorrected run is in `analysis_glint_none/`.
- E_d transmittance 0.21: heavy evening overcast.
