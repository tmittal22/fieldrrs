"""Curate a small, high-level figure set per station and for the whole field day.

    python make_highlights.py Data_NatureSpec/2026_Aug_16

Every station's `analysis/` folder holds 13-14 diagnostic figures plus a GIOP/ subfolder
of 12 more -- the right depth for verifying a result, the wrong depth for finding "the
plot that shows the answer" at a glance or for pulling into a slide deck. This does NOT
prune or delete any of that: every diagnostic figure stays exactly where it is and is
still the evidence `*_GIOP_FINDINGS.md`/`PAPER_READINESS.md` cite by path. It only ADDS
`analysis/highlights/` (headline R_rs + GIOP result, symlinked not copied, so there is
never a second copy to go stale) per station, and a day-level `highlights/` folder
(the site map, the cross-station comparison figures, and the field photos that establish
what kind of target each station's water actually was).
"""

import argparse
import os

# (headline filename, where it lives relative to a station's own analysis/ folder)
STATION_HIGHLIGHTS = [
    "fig1_pooled_measurements.png",     # what was measured, before any correction
    "fig10_final_product.png",          # every treatment, side by side
    "fig12_FINAL_mean_Rrs.png",         # THE headline R_rs spectrum
    "fig14_all_spectra.png",            # every individual scan, not just the mean
    os.path.join("GIOP", "giop10_final_result.png"),   # THE headline GIOP fit
    os.path.join("GIOP", "giop6_all_fits.png"),         # every scan's own fit
]

#: day-level comparison figures already produced by other scripts this session,
#: gathered here rather than regenerated -- COMPARISON/ is the one place that already
#: has them, this just makes them easy to find alongside the map.
DAY_COMPARISON = [
    "by_location/COMPARISON/all_stations_overplot.png",
    "by_location/COMPARISON/loc1_loc2abc_overplot.png",
    "by_location/COMPARISON/LOC2_site_summary.png",
    "by_location/COMPARISON/LOC2_site_means_overlay.png",
    "by_location/COMPARISON/LOC3_site_summary.png",
    "by_location/COMPARISON/LOC3_site_means_overlay.png",
    "by_location/COMPARISON/giop_cross_station_all.png",
]

#: one representative field photo per DISTINCT target type this field day recorded --
#: not one per station, since several stations are the same target type (open water).
#: (station-relative folder, scan number, what it shows)
FIELD_PHOTOS = [
    ("LOC1_66.89718N_162.60290W/FLENS8_FOV08", "00001", "open water (LOC1, the reference station)"),
    ("LOC1_66.89718N_162.60290W/FLENS8_FOV08", "00014", "LAND target: algae-on-concrete (not R_rs -- reflectance)"),
    ("LOC1_66.89718N_162.60290W/FLENS8_FOV08", "00015", "LAND target: bare concrete beside it (not R_rs -- reflectance)"),
    ("LOC2a_66.89677N_162.57953W/FLENS8_FOV08", "00030", "open water (LOC2a, main population)"),
    ("LOC2b_66.89677N_162.57953W/FLENS8_FOV08", "00027", "disturbed/turbid water (LOC2b, sediment settling)"),
    ("LOC2c_66.89677N_162.57953W/FLENS8_FOV08", "00043", "floating algae mat (LOC2c -- opaque target, reflectance not R_rs)"),
    ("LOC3_66.89235N_162.59149W/FIBR15_FOV15", "00055", "open water (LOC3-FIBR15, clean sub-population)"),
    ("LOC3_66.89235N_162.59149W/FIBR15_FOV15_murky", "00058", "high-sediment water (LOC3-FIBR15, murky sub-population)"),
    ("LOC3_66.89235N_162.59149W/FLENS8_FOV08", "00047", "particulate/bottom-visible water (LOC3-FLENS8 -- see LOC3_BOTTOM_CAVEAT.md, sharpest substrate visibility of the 4 FLENS8 photos)"),
]

R_RS_STATIONS = [
    "LOC1_66.89718N_162.60290W/FLENS8_FOV08",
    "LOC2a_66.89677N_162.57953W/FLENS8_FOV08",
    "LOC2b_66.89677N_162.57953W/FLENS8_FOV08",
    "LOC3_66.89235N_162.59149W/FIBR15_FOV15",
    "LOC3_66.89235N_162.59149W/FIBR15_FOV15_murky",
    "LOC3_66.89235N_162.59149W/FLENS8_FOV08",
]


def _relink(src_abs, dst):
    """A relative symlink at dst pointing at src_abs -- portable within the repo (Dropbox
    syncs symlinks as such on macOS/Linux; falls back to a copy if that ever fails)."""
    if os.path.lexists(dst):
        os.remove(dst)
    rel = os.path.relpath(src_abs, os.path.dirname(dst))
    try:
        os.symlink(rel, dst)
    except OSError:
        import shutil
        shutil.copy2(src_abs, dst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("day_folder")
    a = ap.parse_args()
    day = a.day_folder.rstrip("/")
    byloc = os.path.join(day, "by_location")

    n_station = 0
    for station in R_RS_STATIONS:
        srcdir = os.path.join(byloc, station, "analysis")
        hl = os.path.join(srcdir, "highlights")
        os.makedirs(hl, exist_ok=True)
        for fn in STATION_HIGHLIGHTS:
            src = os.path.join(srcdir, fn)
            if not os.path.exists(src):
                print("  ! missing, skipped: %s" % src)
                continue
            _relink(os.path.abspath(src), os.path.join(hl, os.path.basename(fn)))
        n_station += 1
        print("wrote %s/highlights/ (%d figures)" % (srcdir, len(STATION_HIGHLIGHTS)))

    dayhl = os.path.join(day, "highlights")
    os.makedirs(dayhl, exist_ok=True)
    _relink(os.path.abspath(os.path.join(day, "fig1_map.png")),
            os.path.join(dayhl, "site_map.png"))
    for rel in DAY_COMPARISON:
        src = os.path.join(day, rel)
        if not os.path.exists(src):
            print("  ! missing, skipped: %s" % src)
            continue
        _relink(os.path.abspath(src), os.path.join(dayhl, os.path.basename(rel)))

    photodir = os.path.join(dayhl, "field_photos")
    os.makedirs(photodir, exist_ok=True)
    manifest = []
    for station, n, what in FIELD_PHOTOS:
        src = os.path.join(byloc, station, "NaturaSpecPlus_SN25494G1_%s.jpg" % n)
        if not os.path.exists(src):
            print("  ! missing photo, skipped: %s" % src)
            continue
        dst_name = "%s_%s.jpg" % (station.split("_")[0], n)
        _relink(os.path.abspath(src), os.path.join(photodir, dst_name))
        manifest.append((dst_name, what))
    with open(os.path.join(photodir, "MANIFEST.txt"), "w") as fh:
        for name, what in manifest:
            fh.write("%s -- %s\n" % (name, what))

    print("wrote %s/ (map + %d comparison figures + %d field photos, %d stations)"
          % (dayhl, len(DAY_COMPARISON), len(manifest), n_station))


if __name__ == "__main__":
    main()
