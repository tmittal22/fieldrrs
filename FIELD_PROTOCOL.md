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

### 135 deg or 90 deg? Read this before you set up

The IOCCG Protocol Series v3.0 (Ch. 5) endorses Mobley's 40 deg / 135 deg as the geometry
that best minimises sun glint, and then immediately qualifies it:

> "the use of phi = 135 deg may easily become the source of perturbations in L_T
> measurements because the radiometer necessarily looks at the sea close to the
> deployment structure or at its shadow. This limitation, which becomes more severe with
> large sun zenith angles, would suggest that phi = 90 deg is a better solution"

Their Fig. 5.1 shows **40 deg / 90 deg** as the geometry "commonly applied".

- **From a boat, pier or any structure**: use **phi = 90 deg**. Looking back at 135 deg
  puts the hull, its shadow or its wake in the field of view, and that ruins the scan
  more thoroughly than the extra glint at 90 deg.
- **From a small craft, a rock or wading, with nothing behind you**: **135 deg** is better.

Set it in the GUI. Whichever you choose, use the SAME geometry all day: the protocol
notes that a single consistent geometry gives more consistent products than switching
geometry to chase conditions.

### Using your phone as the compass

The solar azimuth the GUI computes is a **true** bearing. Your phone compass reads
**magnetic**, and the two differ by the local magnetic declination (about -11 deg in
Pennsylvania, tens of degrees elsewhere). Pointing 11 deg off is a real error in the
relative azimuth.

You do not need to look declination up. Calibrate against the sun itself:

1. Press **WHERE IS THE SUN?** in the GUI. It gives the true solar azimuth.
2. Point the phone at the sun's bearing (at its horizontal direction; **do not look at
   the sun**, and do not sight it through the instrument optics) and read the compass.
3. declination = true solar azimuth − compass reading.

That one sighting absorbs both the local declination and any constant offset in the
phone's magnetometer, which a published declination value does not. Enter it and the GUI
prints the magnetic bearings your phone should actually read. Redo it if you move a long
way, or after standing near anything ferrous: engine block, railing, winch.

Your phone's inclinometer or level app also sets the 50 deg tilt below horizontal.

### Wind speed: what "velocity instrument" means, and what you actually need

Wind speed is standard required metadata in every above-water protocol, including IOCCG
v3.0, which lists "sea state; wind speed and direction; air and water temperature" plus
photographs of the conditions. It is not bureaucratic: rho is a function of
(theta, phi, theta_sun, **W**), so a measurement without a wind speed cannot be
reprocessed with a better rho later.

**You do not need a dedicated instrument.** In order of preference:

1. **Handheld anemometer** — a basic one is about USD 25, a Kestrel is the research
   standard. Best if you will do this repeatedly.
2. **The boat's instruments**, if there are any.
3. **Beaufort estimate from the sea surface.** Standard, accepted, and free:

| Beaufort | wind | what the water looks like | rho status |
|---|---|---|---|
| 0 | < 0.5 m/s | mirror | 0.028 fine |
| 1 | 0.5-1.5 m/s | ripples, no foam crests | 0.028 fine |
| 2 | 2-3 m/s | small wavelets, glassy crests, no breaking | 0.028 fine |
| 3 | 3.5-5 m/s | large wavelets, crests begin to break, **scattered whitecaps** | at the limit |
| 4 | 5.5-8 m/s | small waves, **frequent whitecaps** | **0.028 invalid** |
| 5+ | > 8 m/s | moderate waves, many whitecaps, spray | do not measure |

**The practical rule: when whitecaps start appearing, you have reached the edge of where
rho = 0.028 is valid.** Beyond that, record the sea state, take the measurement if you
must, and flag it. The GUI will tell you the same thing if you enter the wind speed.

4. **A weather app or nearby station** — adequate on an open coast, poor in sheltered or
   fetch-limited water where local wind differs from the regional value.

Also worth recording per IOCCG: **wind direction**, and a **photograph of the sky and the
water** at each station. A photo settles arguments later about cloud and sea state that
no number can.

### Fully overcast: not a lost day

A uniformly overcast sky is a **usable measurement condition**, and in one respect it is
better than clear sky. What changes:

**Better under full cloud**

- **No sun glint.** The specular sun beam is the dominant error in above-water
  radiometry, and under thick cloud there is no beam. The single hardest thing to avoid
  goes away.
- **rho is more stable and barely wind-dependent.** The wind sensitivity of rho comes
  from wave facets sampling *different parts of a non-uniform sky*: tilt a facet under a
  clear sky and it swings between bright horizon and dark zenith. Under a uniform sky
  there is little to sample between, so facet orientation stops mattering. rho = 0.028 is
  a defensible choice under overcast even in wind that would rule it out on a clear day.

**Worse under full cloud**

- **E_d is much lower**, so signal-to-noise drops. Increase integration time and take
  more replicates.
- **No satellite match-up is possible.** The satellite cannot see the water through the
  cloud either.
- **Do NOT apply the BRDF normalisation.** The Morel f/Q tables describe a clear-sky
  light field with a direct beam. Under full overcast the field is entirely diffuse and
  those tables do not describe it.

**Changes to the protocol**

- The **135° relative azimuth stops meaning anything** — there is no sun direction to
  point away from. Keep the 40° view angle, and pick the bearing purely to avoid the
  platform, its shadow and its wake.
- Keep everything else: three scans, level panel, sky at the mirror angle.
- Record **"overcast"** in your notes. It changes how rho and the BRDF should be treated
  downstream, and nothing in the file records the sky state for you.

### ⚠ Broken cloud is the worst case, worse than either extreme

Not thick cloud: **patchy, moving cloud.** The panel method assumes E_d is the same at
the moment of the panel scan and the moment of the target scan. Under broken cloud it is
not, and the error is unbounded and invisible in the output.

If the sky is patchy, either wait for it to become uniformly overcast or uniformly clear,
or use the irradiance channel below.

### Your separate irradiance instrument solves exactly this

The irradiance sensor is a **second physical instrument** with a wide (hemispherical)
field of view, logging E_d at the **same moment** as the radiance scans. It is not a
diffuser foreoptic swapped onto the NaturaSpec. In the GUI it is the fourth slot,
**LOAD ED**, and loading it switches which physics runs:

| E_d slot | what runs |
|---|---|
| empty | `E_d = π L_panel / R_panel` (the panel route) |
| filled | E_d measured; **the panel drops out of R_rs entirely** |
| both filled | the above, **plus** an inter-instrument cross-calibration |

```python
res = rrs_from_separate_ed(wavelength, l_target, l_sky, ed_measured, ed_wavelength)
```

The two instruments need not share a wavelength grid; E_d is interpolated onto the
radiance grid, and radiance bands outside the irradiance sensor's range are reported
rather than silently extrapolated.

**What it buys:**

- **No panel-to-target time lag.** E_d is measured with the target, so changing light
  stops being an error. This is the fix for broken cloud.
- **The panel reflectance drops out entirely.** No 0.99 assumption, no certificate.
- **Panel levelness stops mattering.** The collector defines the horizontal plane itself.

**⚠ What it costs, and this is new.** R_rs now divides a radiance measured by instrument A
by an irradiance measured by instrument B, so **any offset between their absolute
calibrations is a direct multiplicative bias on every R_rs**. It does not average out and
nothing in the spectrum reveals it. Measured: a 6 % gain offset moves R_rs(443) by
**−5.7 %**, uniformly and silently.

**So keep taking the panel scan.** It stops being the E_d source and becomes the
**transfer standard** that ties the two instruments together:

```
C(λ) = [π · L_panel(λ) / R_panel(λ)] / E_d_measured(λ)
```

Load the panel and the E_d file together and the GUI computes C and prints a verdict.
A **flat** disagreement behaves like a gain offset. A **spectrally structured** one is not
a gain at all and points at a stale calibration, a tilted collector or a partly shaded
panel. C is **reported, not applied automatically**: multiplying it back in re-introduces
the panel reflectance the separate sensor was there to remove, so that is your call.
Run it **once per deployment**, and again after anything gets knocked.

Conditions: the irradiance channel must be radiometrically calibrated in W m⁻² nm⁻¹, and
the collector must be **level** for the same reason the panel had to be. Its cosine
response error is better behaved under diffuse light than under a low direct sun, so
overcast is where it performs best. A file with no irradiance column is **refused**, not
silently used: dividing by a radiance would be wrong by a factor of π and would look
entirely plausible in the output.

Verified by test: the measured-E_d path recovers a known R_rs across mismatched
wavelength grids, contains no panel term at all, and the cross-calibration recovers an
injected 6 % gain offset as C = 1/1.06 with zero spectral spread.

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
