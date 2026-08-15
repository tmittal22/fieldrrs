# FIELD PROTOCOL — above-water R_rs with the NaturaSpec Plus

Print this. The software cannot fix the wrong three scans.

---

## Before you leave

- [ ] Python installed on the tablet, **"Add python.exe to PATH" ticked**
- [ ] `run_gui.bat` double-clicked once, GUI opens, closes cleanly
- [ ] DARWin radiometric calibration set to **RADIANCE**, units `W/m^2/sr/nm`
- [ ] Spectralon panel clean, dry, in its case. Note its reflectance (0.99 typical;
      use the calibration certificate value if you have one)
- [ ] Instrument battery charged, spare charged
- [ ] Compass or phone compass (you need sun azimuth)
- [ ] Anemometer, or a way to record wind speed. **This matters more than you think**
- [ ] Notebook: station, time, wind, cloud, sea state, water appearance

---

## The three scans

You cannot measure R_rs directly. Pointed at the water, the instrument sees
water-leaving radiance **plus reflected skylight**. Three scans separate them.

```
                          SUN
                           \
                            \   135 deg in azimuth
                             \      between sun and your view
        SKY SCAN               \
        40 deg from ZENITH      \
              \                  \
               \  *               \
                \   *              \
                 \    *             \
   ~~~~~~~~~~~~~~~\~~~~~*~~~~~~~~~~~~\~~~~~~~~~~  water surface
                        *
                         *  WATER SCAN
                          *  40 deg from NADIR
```

| scan | point at | what it gives |
|---|---|---|
| **PANEL** | white reference panel, level, in full sun, no shadow on it | E_d = pi x L_panel / R_panel |
| **SKY** | sky at **40 deg from zenith**, **same compass bearing** as the water view | the skylight that reflects into your water view |
| **WATER** | water at **40 deg from nadir**, **135 deg in azimuth from the sun** | L_t = water + reflected sky |

### About the panel

**"Spectralon"** is a brand name (Labsphere) for a sintered PTFE plaque: a hard white
disc or square that reflects about 99 % of the light hitting it, and reflects it almost
perfectly diffusely, across 350-2500 nm. In practice it means **the white reference panel
that came with the instrument**. Generic equivalents are Zenith Polymer or a
barium-sulfate plate. Any of them works as long as you know its reflectance.

- Enter its reflectance in the GUI. **0.99** is the standard value; if your panel came
  with a calibration certificate, use that number. It is a straight multiplier on every
  R_rs you produce, so a 5 % error in it is a 5 % error in everything.
- Keep it clean and dry. Do not touch the white face. A fingerprint or a water drop is a
  local reflectance error you cannot correct afterwards.

**"Horizontal" means level, face up, perpendicular to gravity.** Not tilted toward you,
not tilted toward the sun, not held at an angle in your hand.

E_d is *by definition* the downwelling irradiance onto a horizontal plane, and the
equation `E_d = pi x L_panel / R_panel` is only true if the panel is in that plane. Tilt
it and it intercepts a different amount of sunlight, by the cosine of the tilt. With a
low sun a 5 deg tilt is a several-percent error that propagates into every wavelength.

- Put it on a flat surface, or use a bubble level (many panel cases have one; a phone
  level app is fine).
- View it from **near-nadir**, within about 10-15 deg of straight down.
- Stand so that neither your shadow nor your reflection lands on it. Dark clothing helps.
  A bright jacket leaning over the panel adds its own light to your E_d.

Then: **L_w = L_t − rho x L_sky**, and **R_rs = L_w / E_d**

The sky scan is the mirror of the water scan through the horizontal. Water 40 deg down
on bearing X, sky 40 deg up on bearing X. Same bearing, not the opposite one.

**Why 40 / 135.** Mobley (1999). 135 deg from the sun keeps you out of the specular
glint lobe without looking straight downsun; 40 deg from nadir is where rho is smallest
and most stable and where you avoid your own shadow. If you change the geometry, rho
changes and the number in the software is no longer right.

---

## DARWin workflow at each station

DARWin has a **reference** scan and **target** scans. Use that: scan the panel as the
REFERENCE, then take sky and water as TARGETS. Every saved file then carries the panel
in its `Rad. (Ref.)` column, so the panel is automatically contemporaneous with the
target and you do not need a separate panel file.

1. **Panel** — level, unshaded, you standing so your shadow misses it. Take the
   REFERENCE scan.
2. **Sky** — 40 deg from zenith, on the water bearing. TARGET scan. Save.
3. **Water** — 40 deg from nadir, 135 deg from the sun. TARGET scan. Save.
4. **Repeat 2 and 3** four more times. Five replicates minimum.
5. **Re-scan the panel** as reference every ~10 minutes, and immediately whenever the
   light changes (cloud crosses the sun, you move station).

**Name your files so the GUI can pre-sort them.** Put `water`, `sky` or `panel` in the
DARWin Comment field or the filename. The GUI guesses the role from those words and you
confirm with one tap. It never assumes silently.

---

## Record for every station

| field | why it matters |
|---|---|
| **wind speed (m/s)** | sets rho, which is the largest error in the whole measurement |
| cloud cover, sun obscured y/n | rho assumes clear sky |
| time (local + UTC) | solar zenith |
| sun azimuth, view azimuth | confirms the 135 deg |
| water appearance: clear / green / brown / sediment | decides which glint correction is legal |
| whitecaps present y/n | whitecaps break the method |
| panel reflectance used | direct multiplier on every number |

---

## What ruins a measurement

- **Looking into the sun glint.** If you can see the sun's reflection in your view, move.
- **Ship, platform or your own shadow** on the water patch or on the panel.
- **Whitecaps or foam** in the field of view.
- **Sky scan on the wrong bearing.** Most common error. Same bearing as the water, not
  opposite, not toward the sun.
- **The light changing between panel and target.** Re-scan the panel.
- **Panel not level**, or tilted toward you.
- **Wind above ~5 m/s.** The rho = 0.028 default stops being valid. Record the wind and
  flag the station; the software will tell you it cannot supply the right rho.

---

## Sanity checks you can do in the field, on the tablet

Run the GUI at the station, not at home. Compute after the first station and look at:

1. **Is R_rs positive across 400–700 nm?** Negative visible values mean over-subtraction:
   rho too high, wrong sky bearing, or the glint correction is wrong for this water.
   The GUI prints a warning; do not ignore it.
2. **Does the spectrum peak where the water looks?** Blue-peaking for clear water,
   green-peaking (~550–570) for productive or turbid water, red shoulder for sediment.
   If clear blue water gives a green peak, check the sky scan.
3. **Do the five replicates lie on top of each other?** If they scatter, the surface is
   too variable; take more, or move.
4. **Magnitude.** Open ocean R_rs(443) is order 1e-3 to 1e-2 sr^-1. Coastal or turbid
   is higher. If you get 0.1 or 1e-6, something is wrong.

---

## Which residual glint correction to use

Set this in the GUI. It is a property of **your water**, not a software preference.

| water | use |
|---|---|
| clear, oceanic, blue | `nir_zero` is defensible |
| turbid, sediment, river-influenced, brown | **do not use `nir_zero`** — it deletes real signal. Use `nir_similarity` or `none` |
| unsure | `none`, and decide later at the desk. You can always reprocess |

`none` is the default deliberately. Recording the raw scans means you can change your
mind; a correction applied in the field cannot be undone if you only keep the CSV.

---

## Keep the raw files

Copy the whole `.sed` folder off the tablet at the end of the day, before anything else.
The CSVs are derived; the `.sed` files are the measurement.
