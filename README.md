# fieldrrs

Spectral Evolution `.sed` scans in, remote-sensing reflectance out. Runs on a Windows
field tablet with **nothing installed but Python**. No pip, no internet, no numpy.

**Print [`FIELD_CARD.pdf`](FIELD_CARD.pdf)** (A3 landscape) and take it with you: side
view with the angles, plan view with the compass bearings, and the three-step sequence.
Read [`FIELD_PROTOCOL.md`](FIELD_PROTOCOL.md) before collecting data.

Regenerate the card with `python make_field_card.py`. That script is the only thing in
the repository that needs matplotlib; the field package itself stays pure standard
library.

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

**Single station.** Three buttons at the top, one per scan:

- **LOAD WATER** — the 40°-from-nadir, 135°-from-sun scan
- **LOAD SKY** — the 40°-from-zenith scan on the same bearing
- **LOAD PANEL** — *optional*. Leave it empty and the panel radiance is read from the
  water file's own `Rad. (Ref.)` column, which is the DARWin reference-scan workflow.

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
- `delta` is the residual glint offset: `none` (default), `nir_zero` (clear water only),
  or `nir_similarity` (Ruddick et al. 2006, usable in turbid water).

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

52 tests, standard library only. They build `.sed` files from a **known** R_rs, read them
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
