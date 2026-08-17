"""Reorganise a day's scans into folders by LOCATION and FOREOPTIC.

    python organize_by_location.py Data_NatureSpec/2026_Aug_16 [--apply]

Prints the proposed layout by default and writes nothing. `--apply` COPIES files into
the new tree; the originals are never moved or deleted.

Why foreoptic is part of the key, not just a note:

    This day mixes a FLENS8 (8 deg) lens with a FIBR15 (15 deg) fibre. A station mean
    that averages an 8 deg sky scan with a 15 deg water scan is not a measurement of
    anything -- the two see different solid angles, so they average over different
    footprints and different wave-facet populations, and rho is not the same for both.
    Grouping on the panel reference alone silently mixed them here.

Location clustering uses a 60 m tolerance, which is far larger than the GPS scatter
within one occupation (metres) and far smaller than the separation between occupations
(hundreds of metres), so the answer is insensitive to the exact threshold.
"""

import argparse
import csv
import math
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fieldrrs.sed import read_sed
from process_field_day import classify

TOL_M = 60.0


def foreoptic(spec):
    """'FLENS8 {RADIANCE},FLENS8 {RADIANCE}' -> 'FLENS8'. Both channels always agree."""
    raw = spec.header.get("Foreoptic", "?")
    return raw.split("{")[0].strip().rstrip(",").strip() or "?"


def fov_deg(name):
    """Full-angle field of view. FLENS8 is 8 deg; FIBR15 is a 15 deg fibre."""
    digits = "".join(c for c in name if c.isdigit())
    return float(digits) if digits else float("nan")


def survey(folder):
    out = []
    for f in sorted(os.listdir(folder)):
        if not f.lower().endswith(".sed"):
            continue
        spec = read_sed(os.path.join(folder, f))
        role, diag = classify(spec)
        out.append({"file": f, "stem": os.path.splitext(f)[0],
                    "n": f.split("_")[-1][:5], "spec": spec, "role": role,
                    "diag": diag, "fo": foreoptic(spec),
                    "lat": spec.latitude, "lon": spec.longitude,
                    "gps": spec.gps_time, "range": spec.range_m,
                    "sun": spec.solar_elevation_deg})
    return out


def cluster(scans, tol_m=TOL_M):
    lat0 = sum(s["lat"] for s in scans) / len(scans)
    kx = 111320.0 * math.cos(math.radians(lat0))
    locs = []
    for s in scans:
        for c in locs:
            if math.hypot((s["lon"] - c["lon"]) * kx,
                          (s["lat"] - c["lat"]) * 111320.0) < tol_m:
                c["scans"].append(s)
                c["lat"] = sum(x["lat"] for x in c["scans"]) / len(c["scans"])
                c["lon"] = sum(x["lon"] for x in c["scans"]) / len(c["scans"])
                break
        else:
            locs.append({"lat": s["lat"], "lon": s["lon"], "scans": [s]})
    locs.sort(key=lambda c: min(x["gps"] for x in c["scans"] if x["gps"] is not None))
    return locs


def label(lat, lon, i):
    return "LOC%d_%07.5fN_%08.5fW" % (i, lat, -lon)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--apply", action="store_true",
                    help="copy files into the new tree (originals untouched)")
    ap.add_argument("--with-photos", action="store_true")
    a = ap.parse_args()

    scans = survey(a.folder)
    locs = cluster(scans)
    print("%d scans -> %d locations (%.0f m tolerance)\n" % (len(scans), len(locs),
                                                             TOL_M))
    manifest = []
    for i, c in enumerate(locs, 1):
        name = label(c["lat"], c["lon"], i)
        gps = [x["gps"] for x in c["scans"] if x["gps"] is not None]
        print("%s   %d scans   %s" % (name, len(c["scans"]),
                                      "%.2f-%.2f UTC" % (min(gps), max(gps))
                                      if gps else "no GPS time"))
        by_fo = {}
        for s in c["scans"]:
            by_fo.setdefault(s["fo"], []).append(s)
        for fo in sorted(by_fo):
            grp = by_fo[fo]
            roles = {}
            for s in grp:
                roles[s["role"]] = roles.get(s["role"], 0) + 1
            usable = roles.get("sky", 0) > 0 and roles.get("water", 0) > 0
            print("    %-8s  FOV %2.0f deg  n=%2d  %-42s  %s"
                  % (fo, fov_deg(fo), len(grp), roles,
                     "OK" if usable else "** NOT SELF-CONTAINED **"))
            for s in grp:
                manifest.append({"file": s["file"], "location": name,
                                 "lat": "%.6f" % s["lat"], "lon": "%.6f" % s["lon"],
                                 "foreoptic": fo, "fov_deg": "%.0f" % fov_deg(fo),
                                 "role": s["role"], "gps_hours": s["gps"],
                                 "range_m": s["range"], "sun_elev_deg": s["sun"],
                                 "blue_green": "%.4f" % s["diag"]["blue_green"],
                                 "nir_vis": "%.4f" % s["diag"]["nir_vis"]})
        print()

    mixed = [c for c in locs if len({s["fo"] for s in c["scans"]}) > 1]
    if mixed:
        print("WARNING  %d location(s) mix foreoptics. Sky and water scans taken with "
              "DIFFERENT\n         fields of view must not be combined: they average "
              "over different\n         footprints and different facet populations, so "
              "rho is not shared.\n" % len(mixed))

    out = os.path.join(a.folder, "by_location")
    if not a.apply:
        print("DRY RUN. Nothing written. Re-run with --apply to create %s/" % out)
        return

    for i, c in enumerate(locs, 1):
        name = label(c["lat"], c["lon"], i)
        for s in c["scans"]:
            dest = os.path.join(out, name, "%s_FOV%02.0f" % (s["fo"], fov_deg(s["fo"])))
            os.makedirs(dest, exist_ok=True)
            shutil.copy2(os.path.join(a.folder, s["file"]),
                         os.path.join(dest, s["file"]))
            if a.with_photos:
                jpg = s["stem"] + ".jpg"
                src = os.path.join(a.folder, jpg)
                if os.path.exists(src):
                    shutil.copy2(src, os.path.join(dest, jpg))
    path = os.path.join(out, "MANIFEST.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0]))
        w.writeheader()
        w.writerows(manifest)
    print("copied %d .sed files into %s/" % (len(manifest), out))
    print("wrote %s" % path)
    print("originals untouched in %s" % a.folder)


if __name__ == "__main__":
    main()
