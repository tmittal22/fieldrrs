"""Write synthetic .sed scans in the EXACT structure of a real NaturaSpec Plus export.

    python make_demo_data.py

Creates ./demo_scans with two stations, each with panel / sky / water files.

WHY THE STRUCTURE MATTERS. These files are modelled on a real
NaturaSpecPlus_SN25494G1 export from 2025-08-07, not on a reconstruction from the
DARWin string table. That means: 2151 bands at 1 nm over 350-2500 nm, CRLF line
endings, UTF-8 with degree signs, the `<Metadata>` USER_FIELD block carrying Range /
Tilt (X) / Tilt (Y) / Solar Angle, a `Columns [4]:` line, an EMPTY Comment, GPS Time
alongside a separate instrument clock, and the four columns
`Wvl`, `Rad. (Ref.)`, `Rad. (Target)`, `Reflect. %`.

Each file is built from a KNOWN R_rs, so the closure test has a right answer. After
pressing COMPUTE you should get:

    station1   R_rs(443) = 0.00294 sr^-1   blue peak near 490 nm
    station2   R_rs(555) = 0.00600 sr^-1   green peak near 560 nm
"""

import math
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "demo_scans")

# Real instrument grid: 350-2500 nm at 1 nm, 2151 channels.
WL = [350.0 + i for i in range(2151)]
RHO = 0.028
PANEL_R = 0.99
LAT, LON = 40.89674, -78.20499          # from the real file
DEG = "\u00b0"


def darwin_float(v, digits=6):
    """DARWin's 3-digit-exponent style, e.g. 1.119097e-002."""
    s = "%.*e" % (digits, v)
    mant, _, exp = s.partition("e")
    return "%se%s%03d" % (mant, exp[0], int(exp[1:]))


def write_sed(path, rad_ref, rad_target, tag, solar_elev, tilt_y):
    lines = [
        "Comment:",
        "Version: 2.4 [2.1.9284]",
        "File Name: C:\\Users\\W0168\\Documents\\SpectralEvolution\\"
        "NaturaSpecPlus_25494G1\\2025_Aug_07\\%s" % os.path.basename(path),
        "<Metadata>",
        "USER_FIELD1: Range: 3.837m",
        "USER_FIELD2: Tilt (X): + 4.8%s" % DEG,
        "USER_FIELD3: Tilt (Y): +%.1f%s" % (tilt_y, DEG),
        "USER_FIELD4: Solar Angle: %.2f%s" % (solar_elev, DEG),
        "</Metadata>",
        "Instrument: NaturaSpecPlus_SN25494G1 [3]",
        "Detectors: 1024,512,512",
        "Measurement: REFLECTANCE",
        "Date: 08/07/2025,08/07/2025",
        "Time: 05:56:02.07,06:07:42.07",
        "Temperature (C): 32.90,20.23,-15.67,34.84,20.23,-15.67",
        "Battery Voltage: 15.29,15.16",
        "Averages: 20,20",
        "Integration: 60,100,30,100,100,30",
        "Dark Mode: AUTO,AUTO",
        "Foreoptic: FLENS8 {RADIANCE},FLENS8 {RADIANCE}",
        "Radiometric Calibration: RADIANCE",
        "Units: W/m^2/sr/nm",
        "Wavelength Range: 350,2500",
        "Latitude: %.5f" % LAT,
        "Longitude: %.5f" % LON,
        "Altitude: 440.84",
        "GPS Time: 13:08:18",
        "Satellites: 15/12",
        "Calibrated Reference Correction File: Built-In [built-in]",
        "Channels: %d" % len(WL),
        "Columns [4]:",
        "Data:",
        "\t".join(["Wvl", "Rad. (Ref.)", "Rad. (Target)", "Reflect. %"]),
    ]
    for w, r, t in zip(WL, rad_ref, rad_target):
        refl = 100.0 * t / r if r > 0 else 0.0
        # Wavelength is a plain decimal in a real file ("350.0"), NOT scientific
        # notation. Writing it as %.1e collapses 350.0 and 351.0 onto the same value
        # and silently shifts every retrieved R_rs; the closure test caught it at 10.7 %.
        lines.append("\t".join(["%.1f" % w, darwin_float(r),
                                darwin_float(t), "%.4f" % refl]))
    # CRLF and UTF-8, as the instrument writes.
    with open(path, "w", encoding="utf-8", newline="\r\n") as fh:
        fh.write("\n".join(lines))


def station(tag, rrs_fn, tilt_y=28.9, solar_elev=32.27):
    # Smooth solar-ish Ed with a realistic magnitude; the deep SWIR is deliberately
    # driven toward zero, as the real file does beyond ~1800 nm.
    ed, l_panel, l_sky, l_water = [], [], [], []
    for w in WL:
        base = 1.0 * math.exp(-((w - 550.0) ** 2) / (2 * 300.0 ** 2)) + 0.2
        swir = math.exp(-max(0.0, w - 1750.0) / 120.0)      # dies off in the SWIR
        e = 0.11 * base * swir
        ed.append(e)
        l_panel.append(e * PANEL_R / math.pi)
        l_sky.append(0.02 * e / math.pi)
    for i, w in enumerate(WL):
        l_water.append(rrs_fn(w) * ed[i] + RHO * l_sky[i])

    for role, tgt in (("panel", l_panel), ("sky", l_sky), ("water", l_water)):
        write_sed(os.path.join(OUT, "%s_%s.sed" % (tag, role)),
                  l_panel, tgt, tag, solar_elev, tilt_y)


def main():
    if not os.path.isdir(OUT):
        os.makedirs(OUT)
    station("station1", lambda w: 4e-3 * math.exp(-((w - 490.0) ** 2) / (2 * 60.0 ** 2)))
    station("station2", lambda w: 6e-3 * math.exp(-((w - 560.0) ** 2) / (2 * 70.0 ** 2))
            + 1.2e-3 * math.exp(-((w - 650.0) ** 2) / (2 * 40.0 ** 2)))
    print("Wrote demo scans to: %s" % OUT)
    print("Structure copied from a real NaturaSpecPlus_SN25494G1 export "
          "(2151 bands, CRLF, UTF-8, <Metadata> block).")
    print("Start the GUI, LOAD WATER + LOAD SKY, press COMPUTE.")
    print("Expect station1 R_rs(443) ~ 0.00294 and station2 R_rs(555) ~ 0.00600 sr^-1.")


if __name__ == "__main__":
    main()
