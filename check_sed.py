"""Diagnose a real .sed file against what fieldrrs expects. Zero dependencies.

    python check_sed.py  path\\to\\scan.sed
    python check_sed.py  path\\to\\folder
    python check_sed.py  path\\to\\folder  --report sed_report.txt

WHY THIS EXISTS. The reader was built from the column vocabulary extracted from
DARWin2.exe, not from a real export, so every assumption in it is an inference until a
genuine file has been through it. This prints exactly what your instrument writes, marks
which columns fieldrrs recognises and which it would silently ignore, and says whether
the three-scan calculation can run on it.

Run it on your own data, then send the report. It contains no spectra, only structure.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fieldrrs.sed import COLUMN_ALIASES, KNOWN_UNUSED, guess_role, read_sed  # noqa: E402

KNOWN = set()
for _aliases in COLUMN_ALIASES.values():
    KNOWN.update(_aliases)

# Header keys DARWin is known to write, from the binary's string table.
EXPECTED_HEADER = [
    "Version", "File Name", "Instrument", "Detectors", "Measurement", "Date", "Time",
    "Temperature (C)", "Battery Voltage", "Averages", "Integration", "Dark Mode",
    "Foreoptic", "Radiometric Calibration", "Units", "Wavelength Range", "Latitude",
    "Longitude", "Altitude", "GPS Time", "Satellites", "Channels", "Comment",
    "Calibrated Reference Correction File",
]

OUT = []


def say(msg=""):
    OUT.append(msg)
    print(msg)


def diagnose(path):
    """Returns (ok, can_do_radiance, can_do_reflectance)."""
    say("=" * 78)
    say("FILE: %s" % os.path.basename(path))
    say("=" * 78)

    try:
        s = read_sed(path)
    except Exception as exc:
        say("  ** FAILED TO PARSE **")
        say("  %s: %s" % (type(exc).__name__, exc))
        say("")
        say("  Send this file (or its first 30 lines) so the reader can be fixed.")
        say("")
        return False, False, False

    say("  PARSED OK")
    say("")

    # ---- header
    say("  HEADER  (%d keys)" % len(s.header))
    for k in sorted(s.header):
        flag = "" if k in EXPECTED_HEADER else "   <-- not in the expected set"
        say("    %-34s = %s%s" % (k, s.header[k][:44], flag))
    missing = [k for k in ("Instrument", "Measurement", "Units", "Comment")
               if k not in s.header]
    if missing:
        say("    NOTE: no %s key. Not fatal, but the GUI uses it for labelling."
            % ", ".join(missing))
    say("")

    # ---- columns
    say("  COLUMNS  (%d)" % len(s.raw_columns))
    for lab in s.raw_columns:
        if lab in KNOWN:
            canon = [c for c, al in COLUMN_ALIASES.items() if lab in al][0]
            say("    %-26s RECOGNISED -> %s" % (repr(lab), canon))
        elif lab in KNOWN_UNUSED:
            say("    %-26s known, not needed (ignored deliberately)" % repr(lab))
        else:
            say("    %-26s ** UNRECOGNISED - ignored. If this column matters, say so **"
                % repr(lab))
    say("")

    # ---- wavelength axis
    wl = s.wavelength
    n = len(wl)
    steps = sorted({round(wl[i + 1] - wl[i], 4) for i in range(min(n - 1, 400))})
    say("  WAVELENGTH AXIS")
    say("    bands            : %d" % n)
    say("    range            : %.2f to %.2f nm" % (wl[0], wl[-1]))
    say("    spacing (first 400 gaps): %s nm"
        % (", ".join(str(x) for x in steps[:6]) + ("..." if len(steps) > 6 else "")))
    say("    monotonic        : %s" % all(wl[i + 1] > wl[i] for i in range(n - 1)))
    covers = wl[0] <= 400 and wl[-1] >= 700
    say("    covers 400-700   : %s%s" % (covers, "" if covers else "   <-- GIOP needs this"))
    say("    reaches 750-800  : %s   (needed for the nir_zero glint correction)"
        % (wl[-1] >= 800))
    say("    reaches 870      : %s   (needed for nir_similarity)" % (wl[-1] >= 870))
    say("")

    # ---- what can be computed
    has_rad = s.has("rad_target") and (s.has("rad_ref") or True)
    has_refl = s.has("reflectance")
    say("  WHAT fieldrrs CAN DO WITH THIS FILE")
    say("    radiance path (preferred): %s" % ("YES" if s.has("rad_target") else "NO"))
    if s.has("rad_target") and not s.has("rad_ref"):
        say("      but no Rad. (Ref.) column, so a SEPARATE panel scan is required")
    elif s.has("rad_ref"):
        say("      Rad. (Ref.) present, so the panel can come from this file itself")
    say("    reflectance path         : %s" % ("YES" if has_refl else "NO"))
    if has_refl:
        say("      scale detected: %s"
            % ("percent, divided by 100" if s.reflectance_scale == 100.0
               else "already a fraction"))
    say("")

    # ---- value sanity
    say("  VALUE SANITY")
    for canon in ("rad_ref", "rad_target", "reflectance"):
        if not s.has(canon):
            continue
        v = s.columns[canon]
        finite = [x for x in v if x == x]
        neg = sum(1 for x in finite if x < 0)
        say("    %-12s min %-12.6g max %-12.6g  negatives %d/%d"
            % (canon, min(finite), max(finite), neg, len(v)))
        if canon == "reflectance":
            mx = max(finite)
            if s.reflectance_scale == 100.0 and mx <= 1.5:
                say("      ** label says percent but values look like a fraction "
                    "(max %.3f). CHECK THIS. **" % mx)
            if s.reflectance_scale == 1.0 and mx > 1.5:
                say("      ** label says fraction but values look like percent "
                    "(max %.1f). CHECK THIS. **" % mx)
    say("")

    # ---- metadata the protocol needs
    say("  FIELD METADATA")
    say("    lat / lon        : %s / %s" % (s.latitude, s.longitude))
    say("    date / time      : %s" % (s.when or "(absent)"))
    say("    comment          : %s" % (s.comment or "(empty)"))
    say("    role guessed as  : %s" % guess_role(s))
    if guess_role(s) == "unassigned":
        say("      put 'water', 'sky' or 'panel' in the filename or DARWin Comment and")
        say("      the GUI will pre-sort your scans automatically")
    say("")

    return True, s.has("rad_target"), has_refl


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    report = None
    if "--report" in sys.argv:
        i = sys.argv.index("--report")
        if i + 1 < len(sys.argv):
            report = sys.argv[i + 1]

    if not args:
        print(__doc__)
        return 2

    target = args[0]
    if os.path.isdir(target):
        paths = [os.path.join(target, f) for f in sorted(os.listdir(target))
                 if f.lower().endswith(".sed")]
        if not paths:
            print("No .sed files in %s" % target)
            return 2
    else:
        paths = [target]

    say("fieldrrs .sed format diagnostic")
    say("%d file(s) to check" % len(paths))
    say("")

    ok = rad = refl = 0
    for p in paths:
        a, b, c = diagnose(p)
        ok += a
        rad += b
        refl += c

    say("=" * 78)
    say("SUMMARY")
    say("  parsed successfully      : %d / %d" % (ok, len(paths)))
    say("  usable via radiance      : %d" % rad)
    say("  usable via reflectance   : %d" % refl)
    if ok == len(paths) and rad == len(paths):
        say("  VERDICT: these files work with fieldrrs as written.")
    elif ok == len(paths):
        say("  VERDICT: parseable, but check the radiance columns above.")
    else:
        say("  VERDICT: at least one file did not parse. Send it so the reader")
        say("           can be corrected before it matters.")
    say("=" * 78)

    if report:
        with open(report, "w") as fh:
            fh.write("\n".join(OUT) + "\n")
        print("\nWrote %s  (structure only, no spectra) -- safe to send." % report)

    return 0 if ok == len(paths) else 1


if __name__ == "__main__":
    sys.exit(main())
