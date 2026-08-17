"""Inventory a day's NaturaSpec .sed files before trying to process them.

    python survey_field_data.py Data_NatureSpec/2026_Aug_16

Prints, for every scan: GPS time, pointing (Y tilt), range, solar elevation, the
radiance level in the visible, and an INFERRED role. Then proposes station groupings.

Nothing here decides anything. Role inference from tilt and brightness is a hypothesis
that the operator has to confirm, because the file format does not record what the
instrument was pointed at. It exists so you can see what you have before you fit it,
which is the step that was missing when 68 files landed with no notes.

Standard library only, like the package.
"""

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fieldrrs.sed import read_sed
from fieldrrs.solar import solar_window_verdict


def read_header_log(folder):
    """Read HeaderLog.csv WITHOUT csv.DictReader, which mis-parses it.

    ``ScanType`` holds "TARGET, ALLOW_SWIR2_UPDATE" -- one logical field containing an
    UNQUOTED COMMA -- so every data row has 23 fields against 22 headers and everything
    after ScanType shifts one place left. DictReader therefore reports
    ``Latitude = -15.6697`` (a detector temperature) and ``Longitude = 66.897`` (the
    real latitude). Silent, and wrong in a way that looks plausible.

    The trailing columns are fixed and unambiguous, so anchor on the RIGHT.
    Everything here is a cross-check anyway: the .sed files carry range, tilt, solar
    elevation, GPS time and position in their own <Metadata> block, and that is what
    this survey actually trusts.
    """
    path = os.path.join(folder, "HeaderLog.csv")
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    for r in rows[1:]:
        if len(r) < 6:
            continue
        gps, rng, xt, yt = r[-4:]
        out[r[0].strip()] = {"ScanType": r[3].strip(), "GPSTime": gps.strip(),
                             "Range": rng, "XTilt": xt, "YTilt": yt}
    return out


def band_mean(spec, col, lo, hi):
    vals = [v for w, v in zip(spec.wavelength, spec.columns[col])
            if lo <= w <= hi and v == v]
    return sum(vals) / len(vals) if vals else float("nan")


def infer_role(ytilt, vis, ref_vis):
    """Hypothesis only.

    Geometry first, because it is recorded rather than inferred: a Y tilt near 40 deg
    from the horizontal is the sky view, near -40 (or 40 past vertical) the water view.
    Brightness breaks ties, since a panel in full sun is far brighter than either.
    """
    if ytilt is None:
        return "?"
    t = float(ytilt)
    if vis > 0.5 * ref_vis:
        return "panel?"
    if t > 20.0:
        return "sky?"
    if t < -20.0:
        return "water?"
    return "?"


def main():
    folder = sys.argv[1] if len(sys.argv) > 1 else "Data_NatureSpec/2026_Aug_16"
    log = read_header_log(folder)
    files = sorted(f for f in os.listdir(folder) if f.lower().endswith(".sed"))
    print("%d .sed files in %s\n" % (len(files), folder))

    rows, refs = [], []
    for f in files:
        spec = read_sed(os.path.join(folder, f))
        stem = os.path.splitext(f)[0]
        h = log.get(stem, {})
        vis = band_mean(spec, "rad_target", 450, 650)
        rvis = band_mean(spec, "rad_ref", 450, 650) if spec.has("rad_ref") else float("nan")
        refs.append(rvis)
        rows.append({
            "n": stem.split("_")[-1], "file": f, "spec": spec,
            "gps": h.get("GPSTime", spec.header.get("GPS Time", "?")),
            "scantype": h.get("ScanType", "?").strip(),
            "ytilt": h.get("YTilt"), "xtilt": h.get("XTilt"),
            "range": h.get("Range"), "vis": vis, "rvis": rvis,
            "sun": spec.solar_elevation_deg,
        })

    ref_vis = max(r for r in refs if r == r) if any(r == r for r in refs) else 1.0

    print("%-6s %-9s %-8s %7s %7s %6s %11s %11s  %s"
          % ("scan", "GPS", "type", "Ytilt", "Xtilt", "range", "L_target", "L_ref",
             "role?"))
    print("-" * 96)
    for r in rows:
        r["role"] = infer_role(r["ytilt"], r["vis"], ref_vis)
        print("%-6s %-9s %-8s %7s %7s %6s %11.3e %11.3e  %s"
              % (r["n"], r["gps"], r["scantype"], r["ytilt"] or "?", r["xtilt"] or "?",
                 r["range"] or "?", r["vis"], r["rvis"], r["role"]))

    # ---- what the reference channel is doing
    uniq = len({tuple(round(v, 10) for v in r["spec"].columns["rad_ref"][:50])
                for r in rows if r["spec"].has("rad_ref")})
    print("\nRad. (Ref.) column: %d distinct reference spectra across %d files"
          % (uniq, len(rows)))
    print("  -> %s" % ("ONE stored reference reused by every scan (the DARWin "
                       "reference-scan workflow: panel=None works)" if uniq == 1
                       else "the reference CHANGES between scans; check which panel "
                            "scan each target belongs to"))

    # ---- geometry and sun
    suns = [r["sun"] for r in rows if r["sun"] is not None]
    if suns:
        tier, msg = solar_window_verdict(sum(suns) / len(suns))
        print("\nSolar elevation %.2f-%.2f deg (mean %.2f) -> %s"
              % (min(suns), max(suns), sum(suns) / len(suns), tier.upper()))
        print("  %s" % msg)

    tilts = sorted({round(float(r["ytilt"]), 1) for r in rows if r["ytilt"]})
    print("\nDistinct Y tilts: %s" % tilts)
    print("Distinct scan types: %s" % sorted({r["scantype"] for r in rows}))
    print("Wavelength grid: %.0f-%.0f nm, %d channels"
          % (rows[0]["spec"].wavelength[0], rows[0]["spec"].wavelength[-1],
             len(rows[0]["spec"].wavelength)))
    print("Columns: %s" % ", ".join(sorted(rows[0]["spec"].raw_columns)))
    print("\nInferred roles: %s"
          % {k: sum(1 for r in rows if r["role"] == k)
             for k in sorted({r["role"] for r in rows})})


if __name__ == "__main__":
    main()
