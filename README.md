# fieldrrs

Spectral Evolution `.sed` scans in, remote-sensing reflectance out. Runs on a Windows
field tablet with **nothing installed but Python**. No pip, no internet, no numpy.

## What this does, in one table

Which instruments you carry decides which path runs. The exe supports **both**, and picks
by whether you loaded an E_d file. Full derivation in **[`THEORY.pdf`](THEORY.pdf)**.

| | **A · NaturaSpec Plus alone** | **B · + separate irradiance sensor** |
|---|---|---|
| scans | panel, sky, water (one at a time) | sky, water + **simultaneous E_d** |
| E_d comes from | the **panel**, via `E_d = π L_p / R_p` | **measured directly** |
| panel needed? | **yes — it IS the irradiance** | no — but load it anyway, see below |
| you get | R_rs | R_rs, **PAR**, atmospheric transmittance, E_d-stability check |
| you carry | panel reflectance, panel levelness, **panel-to-target time lag**, ρ | **inter-instrument calibration**, ρ |
| broken cloud | disqualifying (light moves between scans) | **workable** — E_d is simultaneous |
| in the GUI | load 3 slots | load 4 slots; the E_d slot switches the physics |

**The trade is not free.** Setup B removes the panel and the time lag, but R_rs now divides
a radiance from instrument A by an irradiance from instrument B, so **any offset between
their absolute calibrations is a direct multiplicative bias on every R_rs**. Measured: a
6 % gain offset moves R_rs(443) by **−5.7 %**, uniformly and silently.

**So load the panel in setup B too.** It stops being the E_d source and becomes the
**transfer standard** that ties the two instruments together:

```
C(λ) = [π·L_panel(λ)/R_panel(λ)] / E_d_measured(λ)          (F4)
```

Load panel *and* E_d and the GUI prints the verdict automatically:

```
CROSS-CALIBRATION  mean C = 0.943   spread 0.0 %   (301 bands)
The two instruments differ by -5.7 % (mean C = 0.943) but the disagreement is
spectrally flat, so it behaves like a gain offset. Decide which absolute scale
you trust before correcting; applying C re-introduces the panel reflectance you
were trying to avoid.
```

A **flat** disagreement is a gain offset. A **spectrally structured** one is not a gain at
all and points at a stale calibration, a tilted collector or a shaded panel. C is
**reported, never applied automatically** — multiplying it back in re-introduces the panel
reflectance setup B existed to remove. Run it **once per deployment**, and again after
anything gets knocked.

## When can you measure? Sun 30–60° above the horizon

**PREFERRED: solar elevation 30–60°.** That number is chosen so it is also solar *zenith*
30–60° — the window is symmetric about 45°, so it is right whichever convention you mean.
This matters because the instrument logs "Solar Angle" as an **elevation** with no label,
and elevation/zenith get confused constantly. **USABLE 20–70°**; nothing fails outside
30–60, ρ = 0.028 is just less well determined.

Both ends were **computed**, not adopted — `python analysis_solar_window.py`, plotted in
`figures_solar_window.png` and THEORY.pdf p6:

- **Top end is GLINT, and the direction is the surprise.** A wave facet puts the sun in
  the field of view when its normal bisects the sun and view directions. That required
  tilt is *minimised* at ~66° elevation, at only 15.2°, where Cox & Munk (1954) slope
  statistics put such facets just **1.6σ** from the mean at 5 m/s — common. Going from 60°
  to 30° elevation raises the required tilt to 25.8° and cuts the glint weight **15×**.
  **High sun is the glint problem, not low sun. Noon is the worst time, not the best.**
- **Bottom end is SIGNAL.** E_d on the horizontal goes as sin(elevation): 50 % of overhead
  at 30°, 34 % at 20°, and the airmass roughly doubles over that span so the light is both
  weaker and more diffuse.
- **Not a constraint, though everyone assumes it is: your own shadow.** Reach is the wrong
  test — the shadow runs along the anti-solar bearing while the target sits 45° off it, so
  the miss distance is a fixed **1.74 m lateral**, independent of solar elevation. At 10°
  elevation a 2 m operator throws an 11.3 m shadow and still misses by 1.74 m. A *boat or
  pier* is a bigger object and can shadow the spot — that is the IOCCG argument for 90°
  azimuth, where the clearance is 2.47 m.

So: **mid-morning and mid-afternoon.** `WHERE IS THE SUN?` in the GUI prints the
elevation, the tier, and which limit you are near.

---

**Print [`FIELD_CARD.pdf`](FIELD_CARD.pdf)** — A3 landscape, **two pages, print
double-sided**, and take it with you.

*Page 1*: side view with the angles, plan view with the compass bearings, the
sky-condition decision strip, and the three-step sequence.
*Page 2*: working under cloud, the E_d sensor path, and the products the absolute
calibration supports (PAR, nLw).

**Read [`THEORY.pdf`](THEORY.pdf)** (6 pages, A3) for the physics behind all of it:
the three governing equations and why ρ is not a Fresnel coefficient; the A-vs-B
capability chart above with its error budget; every cross-check and what each one can
and cannot catch; the products an absolute irradiance channel supports; and where the sun has to be.
Read [`FIELD_PROTOCOL.md`](FIELD_PROTOCOL.md) before collecting data.

Regenerate with `python make_field_card.py` and `python make_theory_pdf.py`. Those two
scripts are the only things in the repository that need matplotlib; the field package
itself stays pure standard library.

---

## Getting it onto the tablet: two routes

### Route A — standalone `fieldrrs.exe`, no Python on the tablet

**PyInstaller cannot cross-compile**, so a Windows executable has to be built on
Windows. It cannot be produced from Linux or macOS. Build it once, anywhere, then copy
the single file.

On any Windows machine with Python and internet, double-click **`build_exe.bat`**. It
runs the test suite first, refuses to build from a broken tree, installs PyInstaller,
and produces **`dist\fieldrrs.exe`**. Copy that one file to the tablet; it needs no
Python there.

No Windows machine to hand? Push this folder to GitHub and the included
`.github/workflows/build-exe.yml` builds the exe on a Windows runner and attaches it as
a downloadable artifact under the Actions tab.

The spec builds with a **console window on purpose** (`CONSOLE = True` in
`fieldrrs.spec`). A windowed build looks tidier, but if the exe fails to start at a
station it shows you nothing at all, and a traceback you can read is worth more than a
clean taskbar. Flip the flag if you disagree.

### Route B — run from source (already tested, works today)

1. Copy the whole `fieldrrs` folder to the tablet, e.g. `C:\fieldrrs`.
2. Install Python from <https://www.python.org/downloads/windows/>.
   **Tick "Add python.exe to PATH"** on the first installer screen.
3. Double-click `run_gui.bat`.

That is the entire install. `tkinter` ships with the python.org installer, and everything
else is standard library. If Python is missing, the `.bat` says so and tells you what to do.

**Test it tonight** with the demo data:

```
python make_demo_data.py
```

then open the `demo_scans` folder in the GUI and press COMPUTE. You should get a
blue-peaking spectrum with R_rs(443) near 0.0029 sr^-1 and no warnings.

---

## Using it

**Single station.** Four buttons at the top, one per scan:

- **LOAD WATER** — the 40°-from-nadir, 135°-from-sun scan
- **LOAD SKY** — the 40°-from-zenith scan on the same bearing
- **LOAD PANEL** — *optional*. Leave it empty and the panel radiance is read from the
  water file's own `Rad. (Ref.)` column, which is the DARWin reference-scan workflow.
- **LOAD ED** — *optional*, and it is the **setup A / setup B switch**. Leave it empty
  and E_d comes from the panel (setup A). Load it and E_d is the measured irradiance
  (setup B), the panel drops out of R_rs entirely, and — if a panel is also loaded — the
  two instruments are cross-calibrated and the verdict printed. The file needs a real
  irradiance column (`Irr. (Target)` / `Irr. (Ref.)`); a radiance file is **refused**
  rather than silently used, because that error is a clean factor of π and looks
  entirely plausible in the output. The irradiance sensor may be on its own wavelength
  grid; it is interpolated, and radiance bands outside its range are reported.

Each slot shows the loaded filename, band count and wavelength range, and turns green.
If the filename or DARWin Comment disagrees with the slot you are loading it into (a file
called `..._sky` loaded as WATER), it asks before accepting. Both paths, with and without
an explicit panel file, are tested to give the same R_rs.

**Many stations.** Use **OPEN WHOLE FOLDER (batch)** instead: the table lists every
`.sed`, pre-guesses roles from filename and Comment, and you correct them by selecting
rows and tapping a role button. Multiple water scans become replicates of one sky scan
and the mean is computed and plotted with the individuals.

Then, either way:

3. **Settings** — panel reflectance, wind speed, rho, residual glint correction.
3b. **Geometry** — view zenith (default 40° from nadir) and relative azimuth (default
   135° from the sun), plus latitude, longitude, date, local time and UTC offset.
   *fill from the loaded water scan* takes lat/lon/date/time straight out of the `.sed`
   header. **WHERE IS THE SUN?** then computes the solar zenith and azimuth and prints
   the two compass bearings you can point at, the tilt below horizontal for the water
   scan, and the tilt above horizontal for the sky scan. It warns if the solar zenith is
   outside the usable 20–60° window.
   Tap *check rho vs wind*: above ~5 m s⁻¹ it will tell you the 0.028 default is no
   longer valid rather than quietly using it.
4. **COMPUTE Rrs** — plots every water scan plus the replicate mean, and prints
   R_rs at 443/490/555/670 and every warning.
5. **Save** — one CSV per spectrum, plus a wide `rrs_all_stations.csv` with one row per
   station, which is the file to hand to an inversion later.

Multiple water scans are processed as replicates against the same sky scan, and the mean
is computed and plotted with the individuals.

## The physics

```
E_d(l)  = pi * L_panel(l) / R_panel(l)
L_w(l)  = L_target(l) - rho * L_sky(l)
R_rs(l) = L_w(l) / E_d(l) - delta
```

- `rho` is the effective sea-surface radiance reflectance. **Not** the Fresnel
  coefficient. Default 0.028 (Mobley 1999) for 40°/135°, wind under ~5 m s⁻¹, clear sky.
  Under **uniform overcast** it is more stable and barely wind-dependent; call
  `rho_advice(wind, sky="overcast")`.
- `delta` is the residual glint offset: `none` (default), `nir_zero` (clear water only),
  or `nir_similarity` (Ruddick et al. 2006, usable in turbid water).

`source="irradiance"` uses **measured E_d** from a cosine collector instead of inferring
it from the panel. Preferred whenever the light is changing: it removes the
panel-to-target time lag, the panel reflectance and the panel-levelness error at once.
See `overcast_notes()` and FIELD_PROTOCOL.md for when this matters.

**No wind-dependent rho table is bundled.** Mobley (2015) published one and it is not
redistributable. Above ~5 m s⁻¹ the software refuses to invent a value and tells you what
your options are. This is the largest uncertainty in the measurement and it is not
hidden.

## Two file workflows

- **Two files per station (recommended).** Scan the panel as the DARWin REFERENCE, then
  sky and water as TARGETS. Both files carry the panel in `Rad. (Ref.)`, so no separate
  panel file is needed and the panel is guaranteed contemporaneous. Leave `panel`
  unassigned in the GUI.
- **Three files.** A separate panel scan, assigned as `panel`.

Both are tested to give the same answer.

## From the command line

```python
from fieldrrs import read_sed, rrs_from_sed, write_rrs_csv

water = read_sed("station1_water.sed")
sky   = read_sed("station1_sky.sed")
res   = rrs_from_sed(water, sky, panel_reflectance=0.99, rho=0.028)

for note in res.notes:
    print("WARNING:", note)
write_rrs_csv("station1_rrs.csv", res)
```

## Tests

```
python tests/test_fieldrrs.py
```

58 tests, standard library only. They build `.sed` files from a **known** R_rs, read them
back through the full chain, and check the physics inverts to 9 decimal places, with
controls confirming that a wrong rho, a wrong panel reflectance, and the turbid-water
`nir_zero` failure mode all change the answer in the direction they should. One test
asserts the package imports nothing outside the standard library, because that is the
property the field deployment depends on.

## What this is not, and what to use next

It produces R_rs. It does not invert for IOPs.

The CSVs it writes are read directly by
**[giop-workbench](https://github.com/tmittal22/giop-workbench)**, which runs the GIOP
semi-analytical inversion. It parses the `#` header too, not just the two data columns, so
ρ, the glint method, the geometry, the wind speed and the footprint travel with the
spectrum and the inversion can warn about them. It also offers BRDF normalisation using
the geometry recorded here, which is what you need before comparing against a satellite
product.

```python
from giop.io import read_fieldrrs_csv
fs = read_fieldrrs_csv("station1.csv")
print(fs.review())        # what should change your interpretation
```

## References

- Mobley, C.D. (1999), *Estimation of the remote-sensing reflectance from above-surface
  measurements*, Applied Optics 38(36), 7442–7455, doi:10.1364/AO.38.007442
- Ruddick, K.G., De Cauwer, V., Park, Y., Moore, G. (2006), *Seaborne measurements of
  near infrared water-leaving reflectance*, L&O 51(2), 1167–1179,
  doi:10.4319/lo.2006.51.2.1167
