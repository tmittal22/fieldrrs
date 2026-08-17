"""Turn a day of NaturaSpec scans into R_rs, with figures and a GIOP-ready CSV.

    python process_field_day.py Data_NatureSpec/2026_Aug_16 --out results_20260816

Written against the 2026-08-16 Kotzebue dataset (66.897 N, -162.603 W, 60 scans), whose
structure the instrument does NOT record:

  * There is no role field. Nothing in a .sed says whether it was aimed at the sky, the
    water, or the tundra. Tilt does not separate them either -- Y tilt is 36-50 deg for
    EVERY scan in this dataset, sky and water alike -- so roles are assigned from
    SPECTRAL SHAPE, which is unambiguous (see `classify`).
  * Stations are not labelled. But DARWin stores the panel scan in every file's
    `Rad. (Ref.)` column and rewrites it when you re-reference, so a run of scans
    sharing one reference spectrum IS a station. That is the grouping used here.

Both are inferences. Both are printed, so you can overrule them.

Needs matplotlib for the figures, so it is a development script, not part of the
zero-dependency field package.
"""

import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fieldrrs.rrs import RHO_MOBLEY1999, rrs_three_scan
from fieldrrs.sed import read_sed
from fieldrrs.solar import solar_window_verdict

SKY, WATER, VEG = "sky", "water", "land"
COL = {SKY: "#7fb3d5", WATER: "#1f7a99", VEG: "#2e7d32"}


def band(spec, col, lo, hi):
    v = [x for w, x in zip(spec.wavelength, spec.columns[col]) if lo <= w <= hi]
    return sum(v) / len(v) if v else float("nan")


def classify(spec):
    """Assign a role from spectral shape. Returns (role, diagnostics).

    Three signatures, none of them borderline in this dataset:

    * LAND (an opaque target: concrete, algae, tundra). L(865)/L(450-650) is 1.2-2.0
      against < 0.13 for sky and < 0.10 for water, because water absorbs the NIR and a
      solid surface does not. An order of magnitude of separation.

      `red_edge` = R(750)/R(670) then says whether that surface is VEGETATED. In this
      dataset 00014 is algae on concrete (6.3x, deep chlorophyll absorption at 670) and
      00015 is the bare concrete beside it (1.9x, nearly flat and grey). Both are
      deliberate targets, not stray scans.

      LAND TARGETS DO NOT GO THROUGH THE R_rs PIPELINE, and that is physics rather than
      bookkeeping: rho*L_sky is subtracted because for WATER the reflected sky is a
      contaminant sitting on top of the signal. For an opaque surface the sky is part of
      the ILLUMINATION, so subtracting it would remove real signal. The correct product
      is the hemispherical-directional reflectance factor

          R(lambda) = pi L_target / E_d = R_panel * L_target / L_panel

      which needs no rho and no sky scan at all. See `land_reflectance`.
    * SKY -- Rayleigh. Blue/green is ~2.1, because molecular scattering goes as
      lambda^-4.
    * WATER -- blue/green ~0.55 (the opposite sense) AND the NIR is absorbed, since
      liquid water absorbs strongly beyond 700 nm.

    Deliberately NOT used: tilt. It is recorded, which makes it tempting, but in this
    dataset Y tilt is 36-50 deg for sky AND water scans, so it carries no information
    about what was in the field of view.
    """
    b = band(spec, "rad_target", 440, 490)
    g = band(spec, "rad_target", 540, 580)
    vis = band(spec, "rad_target", 450, 650)
    n865 = band(spec, "rad_target", 850, 880)
    r670 = band(spec, "rad_target", 665, 675)
    r750 = band(spec, "rad_target", 745, 755)
    d = {"blue_green": b / g, "nir_vis": n865 / vis, "L_vis": vis,
         "red_edge": r750 / r670 if r670 else float("nan")}
    if d["nir_vis"] > 0.5:
        return VEG, d
    if d["blue_green"] > 1.3:
        return SKY, d
    return WATER, d


def land_reflectance(spec, panel_reflectance=0.99):
    """Hemispherical-directional reflectance factor of an opaque target.

        R = pi L_target / E_d,   E_d = pi L_panel / R_panel
          = R_panel * L_target / L_panel

    No rho and no sky scan: for a solid surface the sky is illumination, not glint to
    remove. Returns (wavelength, R).
    """
    lt = spec.columns["rad_target"]
    lp = spec.columns["rad_ref"]
    return list(spec.wavelength), [panel_reflectance * t / p if p > 0 else float("nan")
                                   for t, p in zip(lt, lp)]


def station_key(spec):
    """Scans sharing one stored panel reference belong to one station."""
    return round(band(spec, "rad_ref", 450, 650), 6)


def mean_spectrum(specs, col):
    n = len(specs)
    return [sum(s.columns[col][i] for s in specs) / n
            for i in range(len(specs[0].wavelength))]


def load(folder):
    out = []
    for f in sorted(os.listdir(folder)):
        if not f.lower().endswith(".sed"):
            continue
        spec = read_sed(os.path.join(folder, f))
        role, diag = classify(spec)
        out.append({"name": os.path.splitext(f)[0], "n": f.split("_")[-1][:5],
                    "spec": spec, "role": role, "diag": diag,
                    "station": station_key(spec)})
    return out


def process(scans, rho, panel_r, residual):
    """One R_rs per water scan, paired with its station's MEAN sky scan.

    Averaging the sky is deliberate. The sky scans within a station are replicates of a
    slowly varying field, so their mean has less noise than any one of them, and the
    alternative (nearest in time) would inject the noise of a single scan into every
    result. Where the sky is genuinely changing that assumption fails, and the printed
    sky scatter is what tells you.
    """
    stations, results = {}, []
    for s in scans:
        stations.setdefault(s["station"], []).append(s)

    for k in sorted(stations, reverse=True):
        members = stations[k]
        skies = [m for m in members if m["role"] == SKY]
        waters = [m for m in members if m["role"] == WATER]
        if not skies or not waters:
            print("  station ref=%.4f SKIPPED: %d sky, %d water"
                  % (k, len(skies), len(waters)))
            continue
        wl = members[0]["spec"].wavelength
        l_sky = mean_spectrum([m["spec"] for m in skies], "rad_target")
        l_panel = members[0]["spec"].columns["rad_ref"]

        # how much did the sky move across its own replicates?
        vis = [m["diag"]["L_vis"] for m in skies]
        spread = (max(vis) - min(vis)) / (sum(vis) / len(vis))

        st = {"ref": k, "n_sky": len(skies), "n_water": len(waters),
              "sky_spread": spread, "wl": wl, "l_sky": l_sky, "l_panel": l_panel,
              "gps": members[0]["spec"].gps_time,
              "sun": members[0]["spec"].solar_elevation_deg, "rrs": []}
        for w in waters:
            res = rrs_three_scan(wl, w["spec"].columns["rad_target"], l_sky, l_panel,
                                 panel_r, rho, residual)
            st["rrs"].append({"n": w["n"], "res": res, "spec": w["spec"]})
        results.append(st)
    return results


def qc(st):
    """Physical checks that would catch a wrong pairing or a bad scan."""
    out = []
    for r in st["rrs"]:
        wl, v = r["res"].wavelength, r["res"].rrs
        vis = [x for w, x in zip(wl, v) if 400 <= w <= 700]
        neg = sum(1 for x in vis if x < 0)
        peak = wl[max(range(len(v)), key=lambda i: v[i] if 400 <= wl[i] <= 750 else -9)]
        r443 = [x for w, x in zip(wl, v) if 440 <= w <= 446]
        nir = [x for w, x in zip(wl, v) if 800 <= w <= 900]
        out.append({"n": r["n"], "neg_frac": neg / len(vis), "peak": peak,
                    "rrs443": sum(r443) / len(r443),
                    "nir": sum(nir) / len(nir)})
    return out


def hhmm(hours):
    """gps_time is a float in hours; print it as a clock."""
    if hours is None:
        return "?"
    h = int(hours)
    return "%02d:%02d" % (h, int(round((hours - h) * 60)))


def figures(results, scans, outdir):
    os.makedirs(outdir, exist_ok=True)
    paths = []

    # --- 1. every R_rs, by station
    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 4.4), sharey=True)
    axes = [axes] if n == 1 else list(axes)
    # set_xlim does NOT constrain autoscaling in y, and this instrument records to
    # 2500 nm where R_rs is meaningless noise. Without an explicit ylim the SWIR sets
    # the scale and the entire visible signal collapses onto the zero line.
    vis_max = max(v for st in results for r in st["rrs"]
                  for w, v in zip(r["res"].wavelength, r["res"].rrs) if 350 <= w <= 950)
    for ax, st in zip(axes, results):
        for r in st["rrs"]:
            ax.plot(r["res"].wavelength, r["res"].rrs, lw=0.9, alpha=0.65,
                    color=COL[WATER])
        mean = [sum(r["res"].rrs[i] for r in st["rrs"]) / len(st["rrs"])
                for i in range(len(st["wl"]))]
        ax.plot(st["wl"], mean, lw=2.4, color="#c0392b",
                label="station mean (n=%d)" % len(st["rrs"]))
        ax.axhline(0, color="#888", lw=0.8)
        ax.set_xlim(350, 950)
        ax.set_ylim(-0.15 * vis_max, 1.15 * vis_max)
        ax.set_xlabel("wavelength (nm)")
        ax.set_title("station ref=%.4f\n%s UTC · sun %.1f deg"
                     % (st["ref"], hhmm(st["gps"]), st["sun"] or float("nan")),
                     fontsize=9.5)
        ax.legend(fontsize=8); ax.grid(alpha=0.25)
    axes[0].set_ylabel("$R_{rs}$ (sr$^{-1}$)")
    fig.suptitle("Above-water $R_{rs}$ — every water scan, grouped by panel reference",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = os.path.join(outdir, "fig1_rrs_by_station.png")
    fig.savefig(p, dpi=140); plt.close(fig); paths.append(p)

    # --- 2. the raw radiances the classification rests on
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, (title, col, key) in zip(axes, [
            ("Panel  $L_{ref}$  (the irradiance reference)", "#d9534f", "l_panel"),
            ("Sky  $L_{sky}$  (station mean)", COL[SKY], "l_sky"),
            ("Water  $L_t$  (all scans)", COL[WATER], None)]):
        for st in results:
            if key:
                ax.plot(st["wl"], st[key], lw=1.6, label="ref=%.4f" % st["ref"])
            else:
                for r in st["rrs"]:
                    ax.plot(st["wl"], r["spec"].columns["rad_target"], lw=0.7,
                            alpha=0.55, color=COL[WATER])
        ax.set_xlim(350, 1000); ax.set_yscale("log")
        # same trap as fig 1: on a log axis a handful of near-zero SWIR samples drag
        # the lower limit to 1e-17 and flatten every real curve into one line.
        pos = [v for ln in ax.get_lines()
               for w, v in zip(ln.get_xdata(), ln.get_ydata())
               if 350 <= w <= 1000 and v > 0]
        if pos:
            ax.set_ylim(max(min(pos), max(pos) * 1e-4), max(pos) * 1.6)
        ax.set_xlabel("wavelength (nm)"); ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.25)
        if key:
            ax.legend(fontsize=7)
    axes[0].set_ylabel("radiance (W m$^{-2}$ sr$^{-1}$ nm$^{-1}$)")
    fig.suptitle("The three measured quantities. Note the sky is BLUE-rising and the "
                 "water is not: that is what assigns the roles.", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    p = os.path.join(outdir, "fig2_radiances.png")
    fig.savefig(p, dpi=140); plt.close(fig); paths.append(p)

    # --- 3. the classifier, as a scatter
    fig, ax = plt.subplots(figsize=(6.9, 5.4))
    for role in (WATER, SKY, VEG):
        xs = [s["diag"]["blue_green"] for s in scans if s["role"] == role]
        ys = [s["diag"]["nir_vis"] for s in scans if s["role"] == role]
        ax.scatter(xs, ys, s=52, alpha=0.8, color=COL[role], edgecolor="k",
                   linewidth=0.4, label="%s (n=%d)" % (role, len(xs)), zorder=3)
    for s_ in scans:
        if s_["role"] == VEG:
            ax.annotate(s_["n"], (s_["diag"]["blue_green"], s_["diag"]["nir_vis"]),
                        fontsize=7, xytext=(4, 3), textcoords="offset points")
    ax.axhline(0.5, color=COL[VEG], ls="--", lw=1.2)
    ax.axvline(1.3, color=COL[SKY], ls="--", lw=1.2)
    ax.text(0.06, 0.55, "vegetation: NIR/vis > 0.5", fontsize=8.5, color=COL[VEG],
            transform=ax.get_yaxis_transform())
    ax.set_xlabel("blue/green  $L$(440-490) / $L$(540-580)")
    ax.set_ylabel("NIR/vis  $L$(850-880) / $L$(450-650)")
    ax.set_title("Roles separate cleanly; the margin to vegetation is 2.3x\n"
                 "(tilt is NOT used: Y tilt is 36-50 deg for sky and water alike)",
                 fontsize=10.5)
    ax.set_yscale("log"); ax.set_xscale("log")
    ax.grid(alpha=0.25); ax.legend(fontsize=9)
    p = os.path.join(outdir, "fig3_classifier.png")
    fig.savefig(p, dpi=140); plt.close(fig); paths.append(p)
    return paths


def write_csv(results, outdir, wl_lo=350.0, wl_hi=950.0):
    """One row per wavelength, one column per station mean. Feeds giop-workbench."""
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "rrs_stations.csv")
    wl = results[0]["wl"]
    cols, names = [], []
    for st in results:
        cols.append([sum(r["res"].rrs[i] for r in st["rrs"]) / len(st["rrs"])
                     for i in range(len(wl))])
        names.append("Rrs_ref%.4f" % st["ref"])

    # Trim to the range that is actually a measurement. The panel reference is EXACTLY
    # ZERO from 1653 nm up in this dataset (SWIR2 has no calibrated reference scan), so
    # R_rs there is NaN, and beyond ~900 nm the water radiance is at the noise floor
    # anyway. Writing those rows out would hand a downstream inversion 848 NaNs per
    # station and let it decide what to do with them, which is not its job.
    keep = [i for i, lam in enumerate(wl)
            if wl_lo <= lam <= wl_hi and all(c[i] == c[i] for c in cols)]
    dropped = len(wl) - len(keep)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["wavelength"] + names)
        for i in keep:
            w.writerow(["%.1f" % wl[i]] + ["%.8e" % c[i] for c in cols])
    print("  CSV: %d bands written (%.0f-%.0f nm); %d dropped as out-of-range or "
          "non-finite" % (len(keep), wl[keep[0]], wl[keep[-1]], dropped))
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--out", default="results_field")
    ap.add_argument("--rho", type=float, default=RHO_MOBLEY1999)
    ap.add_argument("--panel-reflectance", type=float, default=0.99)
    ap.add_argument("--residual", default="none",
                    choices=["none", "nir_zero", "nir_similarity"])
    a = ap.parse_args()

    scans = load(a.folder)
    print("%d scans in %s" % (len(scans), a.folder))
    roles = {}
    for s in scans:
        roles[s["role"]] = roles.get(s["role"], 0) + 1
    print("roles from spectral shape: %s" % roles)
    veg = [s["n"] for s in scans if s["role"] == VEG]
    if veg:
        print("EXCLUDED as vegetation/land: %s" % ", ".join(veg))

    suns = [s["spec"].solar_elevation_deg for s in scans
            if s["spec"].solar_elevation_deg is not None]
    if suns:
        tier, msg = solar_window_verdict(sum(suns) / len(suns))
        print("solar elevation %.2f-%.2f deg -> %s" % (min(suns), max(suns),
                                                       tier.upper()))
        print("  %s" % msg)

    print("\nstations (grouped by stored panel reference):")
    results = process(scans, a.rho, a.panel_reflectance, a.residual)

    print("\n%-10s %5s %5s %9s %10s %8s %9s %8s"
          % ("station", "sky", "water", "skyspread", "Rrs(443)", "peak", "Rrs(NIR)",
             "neg%"))
    print("-" * 74)
    for st in results:
        q = qc(st)
        r443 = sum(x["rrs443"] for x in q) / len(q)
        nir = sum(x["nir"] for x in q) / len(q)
        neg = sum(x["neg_frac"] for x in q) / len(q)
        peak = sorted(x["peak"] for x in q)[len(q) // 2]
        print("%-10.4f %5d %5d %8.1f%% %10.5f %8.0f %9.5f %7.1f%%"
              % (st["ref"], st["n_sky"], st["n_water"], 100 * st["sky_spread"],
                 r443, peak, nir, 100 * neg))

    paths = figures(results, scans, a.out)
    csv_path = write_csv(results, a.out)
    print("\nwrote %s" % csv_path)
    for p in paths:
        print("wrote %s" % p)


if __name__ == "__main__":
    main()
