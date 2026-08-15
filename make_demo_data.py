"""Write a folder of synthetic .sed scans so you can test the GUI before the field.

    python make_demo_data.py

Creates ./demo_scans with two stations, each with panel / sky / water files in DARWin's
own format. Station 1 is clear blue water, station 2 is greener and more turbid.

The files are built from a KNOWN R_rs, so after pressing COMPUTE you should get:

    station1   R_rs(443) = 0.00294 sr^-1   blue peak near 490 nm
    station2   R_rs(555) = 0.00600 sr^-1   green peak near 560 nm

If you get those, the whole chain works on this machine. If you get warnings about
negative R_rs or over-subtraction, something is wrong with the install.
"""

import math
import os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "demo_scans")
WL = [350.0 + i for i in range(651)]          # 350-1000 nm at 1 nm
RHO = 0.028
PANEL_R = 0.99


def darwin_float(v, digits=10):
    s = "%.*e" % (digits, v)
    mant, _, exp = s.partition("e")
    return "%se%s%03d" % (mant, exp[0], int(exp[1:]))


def write_sed(path, rad_ref, rad_target, comment):
    lines = [
        "Version: 2.1",
        "File Name: %s" % path,
        "Instrument: NaturaSpecPlus_DEMO",
        "Detectors: 512,256,256",
        "Measurement: RADIANCE",
        "Date: 08/15/2026,08/15/2026",
        "Time: 12:00:00,12:00:30",
        "Temperature (C): 25.00,-5.00,-5.00",
        "Averages: 10,10",
        "Integration: 20,1,1",
        "Radiometric Calibration: RADIANCE",
        "Units: W/m^2/sr/nm",
        "Latitude: 40.7934",
        "Longitude: -77.8600",
        "Comment: %s" % comment,
        "Channels: %d" % len(WL),
        "Data:",
        "\t".join(["Wvl", "Rad. (Ref.)", "Rad. (Target)", "Reflect. %"]),
    ]
    for w, r, t in zip(WL, rad_ref, rad_target):
        lines.append("\t".join([darwin_float(w), darwin_float(r), darwin_float(t),
                                darwin_float(100.0 * t / r)]))
    with open(path, "w", encoding="latin-1") as fh:
        fh.write("\n".join(lines))


def station(tag, rrs_fn):
    ed = [1.0 * math.exp(-((w - 550.0) ** 2) / (2 * 300.0 ** 2)) + 0.2 for w in WL]
    l_panel = [e * PANEL_R / math.pi for e in ed]
    l_sky = [0.02 * e / math.pi for e in ed]
    l_water = [rrs_fn(WL[i]) * ed[i] + RHO * l_sky[i] for i in range(len(WL))]
    write_sed(os.path.join(OUT, "%s_panel.sed" % tag), l_panel, l_panel, "%s panel" % tag)
    write_sed(os.path.join(OUT, "%s_sky.sed" % tag), l_panel, l_sky, "%s sky" % tag)
    write_sed(os.path.join(OUT, "%s_water.sed" % tag), l_panel, l_water, "%s water" % tag)


def main():
    if not os.path.isdir(OUT):
        os.makedirs(OUT)
    station("station1", lambda w: 4e-3 * math.exp(-((w - 490.0) ** 2) / (2 * 60.0 ** 2)))
    station("station2", lambda w: 6e-3 * math.exp(-((w - 560.0) ** 2) / (2 * 70.0 ** 2))
            + 1.2e-3 * math.exp(-((w - 650.0) ** 2) / (2 * 40.0 ** 2)))
    print("Wrote demo scans to: %s" % OUT)
    print("Now start the GUI (run_gui.bat), open that folder, and press COMPUTE.")
    print("Expect station1 R_rs(443) ~ 0.00294 and station2 R_rs(555) ~ 0.00600 sr^-1.")


if __name__ == "__main__":
    main()
