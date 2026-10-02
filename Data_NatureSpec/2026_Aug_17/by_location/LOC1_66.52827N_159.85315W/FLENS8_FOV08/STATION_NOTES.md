- **Drifting occupation.** The boat drifted 308 m in 7 min (60 m clustering split it
  into three fragments). Processed as one station with `--max-span-m 350`.
- **00012 excluded** (`--exclude-water 00012`). Its stored panel reference is a dark
  target (L_ref(450-650) = 0.005 against 0.16-0.19 for the station's other references), so E_d
  is 32-38x too small and R_rs(555) comes out at 0.21. The radiance itself is fine; the
  reference is not.
- **Three panel-reference blocks**: pairing never crosses a block.
- **00004** (tilt 61°, ρ = 0.072) has a sky term larger than the measured radiance at
  443 nm, so it is negative before glint correction, and it dips below zero in the
  720/760/820 nm atmospheric bands after: sky/panel mismatch. It is kept (step 5 rated it
  glint-correctable), but it is the noisiest scan.
- Overcast sky scans 00000-02, 08, 10-11 were classified by the overcast rule in
  `classify()`, confirmed against the photos.
- GIOP not triaged.
