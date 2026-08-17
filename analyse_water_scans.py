"""Per-scan QC for a station where the water is NOT one population.

    python analyse_water_scans.py <by_location/LOC*/FOREOPTIC_FOVnn>

`analyse_location.py` treats the water scans as replicates of one spectrum and reports
their spread. That is right at LOC1, where 12 scans agree in shape to 1.7 %. It is WRONG
at a station with sun glint, floating algae, or more than one target in the field of
view, because the "spread" then mixes measurement error with genuinely different things.

This script makes no such assumption. It fits nothing and averages nothing until it has
established, per scan:

  * WHAT WAS IN THE FIELD OF VIEW -- the scan's own photograph, contact-sheeted with its
    spectrum, because that is the only direct evidence and it settles arguments that the
    spectra alone cannot.
  * WHETHER THE SCANS ARE ONE POPULATION -- pairwise spectral angle on shape alone,
    which is amplitude-blind, so a bright scan and a dim scan of the same water are 0 deg
    apart while two different waters are not.
  * WHETHER A DEVIANT SCAN IS GLINT -- and this is a falsifiable test, not a judgement.
    If the deviation is glint, a NIR-based correction (Ruddick et al. 2006) must collapse
    the scan onto the main group. If it does not, the deviation is something else and
    must not be corrected away.
  * WHAT THE RESIDUAL DIFFERENCE IS MADE OF -- the difference spectrum regressed against
    the candidate explanations (sky reflectance, a flat offset, a chlorophyll shape), so
    the attribution is measured rather than asserted.

Outputs go to `analysis/water_scans/`.

Diagnostics and their provenance
--------------------------------
`nir_ratio`   R_rs(780)/R_rs(870). Ruddick et al. (2006) give 1.912 for water whose NIR
              signal is pure-water residual. Departures indicate NIR-scattering material
              OR an imperfect surface correction; the sign does not distinguish them,
              which is why the collapse test above exists.
`peak700`     R_rs(700)/R_rs(675). Above ~1 indicates particulate backscatter and/or
              chlorophyll fluorescence dominating over the 675 nm chlorophyll absorption.
`blue_slope`  R_rs(412)/R_rs(555). Residual sky reflection is blue-heavy, so this is the
              most sensitive single band ratio to a rho error.
`red_edge`    R_rs(750)/R_rs(675). Vegetation and algae mats show a step here; open water
              does not.
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

from fieldrrs.rrs import (SIMILARITY_780_870, rho_at_angle, rrs_three_scan,
                          view_zenith_from_tilt)
from organize_by_location import survey

LO, HI = 400.0, 900.0
#: A scan is called deviant when its shape sits this many times the main group's own
#: internal spread away from that group. Stated, not tuned.
DEVIANT_FACTOR = 1.6


def _ix(wl, x):
    return int(np.argmin(abs(wl - x)))


def build_rrs(water, sky, wl, panel_r, method="none"):
    """R_rs per water scan, each paired to the sky closest in view geometry."""
    out = {}
    for w in water:
        tv = view_zenith_from_tilt(w["spec"].tilt_y_deg)
        rho = rho_at_angle(tv)
        sk = min(sky, key=lambda k: abs(view_zenith_from_tilt(k["spec"].tilt_y_deg) - tv))
        r = rrs_three_scan(wl, w["spec"].columns["rad_target"],
                           sk["spec"].columns["rad_target"],
                           w["spec"].columns["rad_ref"], panel_r, rho, method)
        out[w["n"]] = dict(rrs=np.array(r.rrs), sky=sk["n"], rho=rho, tv=tv)
    return out


def shape_angles(R, names, wl):
    """Pairwise spectral angle in degrees, on unit-norm spectra (amplitude removed)."""
    m = (wl >= LO) & (wl <= HI)
    A = np.array([R[n]["rrs"][m] for n in names])
    U = A / np.linalg.norm(A, axis=1, keepdims=True)
    return np.degrees(np.arccos(np.clip(U @ U.T, -1.0, 1.0)))


def group_scans(D, names):
    """Split into a main group and deviants by median angle to the whole set.

    Deliberately crude: a single-link clustering on a 12-point set invents structure.
    This takes the median-angle ranking, calls the tightest half the seed group, and then
    admits anything within DEVIANT_FACTOR of that group's own internal spread.
    """
    med = np.median(D, axis=1)
    seed = list(np.argsort(med)[:max(3, len(names) // 2)])
    inner = D[np.ix_(seed, seed)]
    spread = inner[inner > 1e-9].max() if (inner > 1e-9).any() else 0.0
    main = [i for i in range(len(names))
            if np.median(D[i, seed]) <= max(spread, 1e-9) * DEVIANT_FACTOR]
    dev = [i for i in range(len(names)) if i not in main]
    return main, dev, spread


def collapse_test(water, sky, wl, panel_r, names, main, dev):
    """Does a NIR glint correction merge the deviants into the main group?

    THE test in this script. Glint is an additive surface term, so a correction that
    removes it must bring a glint-contaminated scan into the main group. A scan that
    stays out after correction is not glint-contaminated, and correcting it further
    would be manufacturing agreement.
    """
    rows = []
    for meth in ("none", "nir_zero", "nir_similarity"):
        try:
            R = build_rrs(water, sky, wl, panel_r, meth)
        except Exception as exc:
            rows.append((meth, np.nan, {}, str(exc)))
            continue
        D = shape_angles(R, names, wl)
        inner = D[np.ix_(main, main)]
        spread = inner[inner > 1e-9].max() if (inner > 1e-9).any() else 0.0
        rows.append((meth, spread,
                     {names[i]: float(np.median(D[i, main])) for i in dev}, ""))
    return rows


def attribute(diff, wl, sky, water, panel_r):
    """What is the difference spectrum made of? Angle to each candidate explanation."""
    m = (wl >= LO) & (wl <= HI)
    Lsky = np.mean([np.array(x["spec"].columns["rad_target"]) for x in sky], axis=0)
    Lpan = np.mean([np.array(x["spec"].columns["rad_ref"]) for x in sky], axis=0)
    Ed = np.pi * Lpan / panel_r
    cand = {"sky reflection (rho error)": Lsky / np.maximum(Ed, 1e-12),
            "flat offset (what NIR methods remove)": np.ones_like(wl),
            "chlorophyll (Gaussian at 675 nm)":
                np.exp(-0.5 * ((wl - 675.0) / 12.0) ** 2)}

    def ang(u, v):
        u = u[m] / np.linalg.norm(u[m])
        v = v[m] / np.linalg.norm(v[m])
        return float(np.degrees(np.arccos(np.clip(u @ v, -1, 1))))

    out = {}
    for k, v in cand.items():
        c = float((diff[m] @ v[m]) / (v[m] @ v[m]))
        res = diff[m] - c * v[m]
        out[k] = dict(angle=ang(diff, v), scale=c,
                      removed=100 * (1 - np.linalg.norm(res) / np.linalg.norm(diff[m])))
    return out


# ------------------------------------------------------------------ figures
def fig_contact(scans, outdir, tag, verdicts):
    """Every scan's photograph, labelled. The direct evidence, all on one page."""
    n = len(scans)
    ncol = 6
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.7 * ncol, 2.9 * nrow))
    axes = np.atleast_1d(axes).ravel()
    col = {"sky": "#7fb3d5", "water": "#1f7a99", "land": "#2e7d32", "panel": "#d9534f"}
    for ax, s in zip(axes, scans):
        jpg = os.path.join(os.path.dirname(s["path"]), s["stem"] + ".jpg")
        if os.path.exists(jpg):
            ax.imshow(mpimg.imread(jpg))
        else:
            ax.text(0.5, 0.5, "no photo", ha="center", va="center",
                    transform=ax.transAxes, fontsize=9, color="#999")
        v = verdicts.get(s["n"], "")
        ax.set_title("%s  %s\n%s" % (s["n"], s["role"], v or "-"), fontsize=7.5,
                     color=col.get(s["role"], "k"),
                     fontweight="bold" if v and v != "clean" else "normal")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_edgecolor("#c0392b" if v and v not in ("clean", "") else "#ccc")
            sp.set_linewidth(2.4 if v and v not in ("clean", "") else 0.8)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("WHAT WAS IN THE FIELD OF VIEW — every scan, %s\n"
                 "red border = flagged by the spectral diagnostics" % tag, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    p = os.path.join(outdir, "w1_contact_sheet.png")
    fig.savefig(p, dpi=125); plt.close(fig)
    return p


def fig_perscan(R, names, wl, groups, outdir, tag, diag):
    """Each water scan on its own axes, with its own numbers."""
    n = len(names)
    ncol = 4
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 3.3 * nrow), sharex=True)
    axes = np.atleast_1d(axes).ravel()
    m = (wl >= LO) & (wl <= HI)
    mainmean = np.mean([R[names[i]]["rrs"] for i in groups["main"]], axis=0)
    for ax, nm in zip(axes, names):
        r = R[nm]["rrs"]
        dv = nm in groups["dev_names"]
        ax.plot(wl[m], mainmean[m], lw=1.2, color="#bbbbbb", label="main-group mean")
        ax.plot(wl[m], r[m], lw=2.0, color="#c0392b" if dv else "#1f7a99")
        d = diag[nm]
        ax.set_title("%s   sky %s   $\\theta_v$=%.1f$^\\circ$  %s\n"
                     "$R_{rs}$555=%.5f  700/675=%.2f  412/555=%.3f\n"
                     "NIR ratio=%.2f (Ruddick %.3f)   angle to main %.2f$^\\circ$"
                     % (nm, R[nm]["sky"], R[nm]["tv"],
                        "** " + groups["verdict"][nm].upper() + " **" if dv else "clean",
                        d["rrs555"], d["peak700"], d["blue_slope"], d["nir_ratio"],
                        SIMILARITY_780_870, d["angle_main"]),
                     fontsize=8, loc="left",
                     color="#c0392b" if dv else "k")
        ax.grid(alpha=0.25); ax.tick_params(labelsize=8)
    for ax in axes[n:]:
        ax.axis("off")
    axes[0].legend(fontsize=8)
    for ax in axes[max(0, n - ncol):n]:
        ax.set_xlabel("wavelength (nm)")
    for i in range(0, n, ncol):
        axes[i].set_ylabel("$R_{rs}$  sr$^{-1}$")
    fig.suptitle("EVERY WATER SCAN SEPARATELY — %s\ngrey = mean of the %d-scan main "
                 "group, so a red curve departing from grey is a real difference, "
                 "not noise" % (tag, len(groups["main"])), fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = os.path.join(outdir, "w2_per_scan.png")
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_groups(R, names, wl, D, groups, collapse, attrib, outdir, tag):
    fig = plt.figure(figsize=(17.5, 10))

    ax = fig.add_subplot(2, 3, 1)
    im = ax.imshow(D, cmap="magma_r")
    ax.set_xticks(range(len(names))); ax.set_yticks(range(len(names)))
    ax.set_xticklabels([n[-3:] for n in names], rotation=90, fontsize=7)
    ax.set_yticklabels([n[-3:] for n in names], fontsize=7)
    plt.colorbar(im, ax=ax, label="spectral angle (deg)")
    ax.set_title("SHAPE-ONLY distance between scans.\nBlocks = separate populations, "
                 "not scatter", fontsize=10, loc="left")

    ax = fig.add_subplot(2, 3, 2)
    m = (wl >= LO) & (wl <= HI)
    i5 = _ix(wl, 555)
    for i, nm in enumerate(names):
        r = R[nm]["rrs"]
        dv = nm in groups["dev_names"]
        ax.plot(wl[m], (r / r[i5])[m], lw=2.0 if dv else 1.0,
                color="#c0392b" if dv else "#bbbbbb",
                label=nm if dv else None)
    mm = np.mean([R[names[i]]["rrs"] for i in groups["main"]], axis=0)
    ax.plot(wl[m], (mm / mm[i5])[m], lw=2.4, color="#1f7a99", label="main group mean")
    ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("$R_{rs}$ / $R_{rs}$(555)")
    ax.legend(fontsize=7.5); ax.grid(alpha=0.25)
    ax.set_title("Normalised at 555 nm: the SHAPE difference", fontsize=10, loc="left")

    ax = fig.add_subplot(2, 3, 3)
    labs = [r[0] for r in collapse]
    x = np.arange(len(labs))
    ax.plot(x, [r[1] for r in collapse], "o-", lw=2.4, color="#1f7a99", ms=9,
            label="main group's OWN spread")
    for nm in groups["dev_names"]:
        ax.plot(x, [r[2].get(nm, np.nan) for r in collapse], "s--", lw=1.8, ms=8,
                label=nm)
    ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=9)
    ax.set_ylabel("angle to the main group (deg)")
    ax.legend(fontsize=8); ax.grid(alpha=0.25)
    ax.set_title("THE COLLAPSE TEST. If a deviation is GLINT, a NIR\ncorrection must "
                 "bring it down to the blue line. If it does\nnot, it is not glint.",
                 fontsize=10, loc="left")

    ax = fig.add_subplot(2, 1, 2)
    ks = list(attrib)
    xs = np.arange(len(ks))
    ax.bar(xs - 0.2, [attrib[k]["angle"] for k in ks], 0.4, color="#2c6f9b",
           label="angle to the difference spectrum (deg, lower = better)")
    ax.bar(xs + 0.2, [attrib[k]["removed"] for k in ks], 0.4, color="#8a6000",
           label="% of the difference it removes")
    ax.set_xticks(xs); ax.set_xticklabels(ks, fontsize=9.5)
    ax.axhline(0, color="k", lw=1)
    ax.legend(fontsize=9); ax.grid(alpha=0.25, axis="y")
    best = min(ks, key=lambda k: attrib[k]["angle"])
    ax.set_title("WHAT IS THE RESIDUAL DIFFERENCE MADE OF?  Best single explanation: "
                 "%s (%.1f$^\\circ$, removes %.0f %%).\nNone of these is a fit — they "
                 "are one-parameter projections, and a partial removal means the "
                 "difference is NOT one clean thing."
                 % (best, attrib[best]["angle"], attrib[best]["removed"]),
                 fontsize=10, loc="left")

    fig.suptitle("IS THIS ONE WATER OR SEVERAL? — %s" % tag, fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = os.path.join(outdir, "w3_groups_and_glint.png")
    fig.savefig(p, dpi=130); plt.close(fig)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--panel-reflectance", type=float, default=0.99)
    a = ap.parse_args()

    scans = survey(a.folder)
    for s in scans:
        s["path"] = os.path.join(a.folder, s["file"])
    loc = os.path.basename(os.path.dirname(a.folder.rstrip("/")))
    fo = os.path.basename(a.folder.rstrip("/"))
    tag = "%s  ·  %s" % (loc, fo)
    outdir = os.path.join(a.folder, "analysis", "water_scans")
    os.makedirs(outdir, exist_ok=True)

    wl = np.array(scans[0]["spec"].wavelength)
    sky = [x for x in scans if x["role"] == "sky"]
    water = sorted([x for x in scans if x["role"] == "water"], key=lambda x: x["n"])
    names = [w["n"] for w in water]
    if len(water) < 3 or not sky:
        msg = ("need >=3 water and >=1 sky scan; got %d/%d" % (len(water), len(sky)))
        print(msg)
        with open(os.path.join(outdir, "REPORT.txt"), "w") as fh:
            fh.write("%s\n\nSKIPPED -- %s\n\n"
                     "This tool's shape-angle clustering and glint-collapse test need "
                     "at least 3 water scans to distinguish 'one deviant scan' from "
                     "'the population itself'; below that there is nothing to cluster "
                     "against. This is not a processing failure -- it is a structural "
                     "limit of the method at this n, the same reason "
                     "analyse_location.py's fig_sensitivity/fig_variability skip below "
                     "n=4/n=3 (see their own docstrings).\n\n"
                     "Consequence: this sub-station's homogeneity does NOT rest on this "
                     "tool. It rests on whatever established it as its own population in "
                     "the first place (spectral-angle distance from the sibling group it "
                     "was split from, and/or visual/field evidence) -- check the "
                     "station's own split-rationale document, not this file.\n"
                     % (tag, msg))
        print("wrote %s/REPORT.txt (skip reason)" % outdir)
        return

    R = build_rrs(water, sky, wl, a.panel_reflectance, "none")
    D = shape_angles(R, names, wl)
    main, dev, spread = group_scans(D, names)
    groups = dict(main=main, dev=dev,
                  dev_names=[names[i] for i in dev], spread=spread)

    collapse = collapse_test(water, sky, wl, a.panel_reflectance, names, main, dev)
    base = dict(collapse[0][2]); best = dict(collapse[-1][2])
    post_spread = collapse[-1][1]
    verdict = {}
    for nm in names:
        if nm not in groups["dev_names"]:
            verdict[nm] = "clean"
        elif best.get(nm, np.inf) <= post_spread * DEVIANT_FACTOR:
            verdict[nm] = "glint (correctable)"
        else:
            verdict[nm] = "deviant, NOT glint"
    groups["verdict"] = verdict

    diag = {}
    for nm in names:
        r = R[nm]["rrs"]
        diag[nm] = dict(
            rrs555=float(r[_ix(wl, 555)]),
            rrs443=float(r[_ix(wl, 443)]),
            peak700=float(r[_ix(wl, 700)] / r[_ix(wl, 675)]),
            blue_slope=float(r[_ix(wl, 412)] / r[_ix(wl, 555)]),
            nir_ratio=float(r[_ix(wl, 780)] / r[_ix(wl, 870)]),
            red_edge=float(r[_ix(wl, 750)] / r[_ix(wl, 675)]),
            angle_main=float(np.median(D[names.index(nm), main])),
            tilt=R[nm]["tv"], sky=R[nm]["sky"], verdict=verdict[nm])

    notglint = [n for n in groups["dev_names"] if verdict[n] == "deviant, NOT glint"]
    attrib = {}
    if notglint:
        Rc = build_rrs(water, sky, wl, a.panel_reflectance, "nir_similarity")
        d = (np.mean([Rc[n]["rrs"] for n in notglint], axis=0)
             - np.mean([Rc[names[i]]["rrs"] for i in main], axis=0))
        attrib = attribute(d, wl, sky, water, a.panel_reflectance)

    ps = [fig_contact(scans, outdir, tag, verdict),
          fig_perscan(R, names, wl, groups, outdir, tag, diag)]
    if attrib:
        ps.append(fig_groups(R, names, wl, D, groups, collapse, attrib, outdir, tag))

    L = []
    P = L.append
    P("=" * 78)
    P("PER-SCAN WATER QC  ·  %s" % tag)
    P("=" * 78)
    P("%d water, %d sky. Nothing is averaged until the scans are shown to be one"
      % (len(water), len(sky)))
    P("population, because a 'spread' over different targets is not an uncertainty.")
    P("")
    P("IS IT ONE POPULATION?")
    P("  main group: %d scans, internal shape spread %.2f deg" % (len(main), spread))
    P("  deviant   : %s" % (", ".join(groups["dev_names"]) or "none"))
    P("")
    P("THE COLLAPSE TEST -- glint is additive, so a NIR correction must remove it")
    P("  %-16s %14s %s" % ("method", "main spread", "   ".join(groups["dev_names"])))
    for meth, sp, dd, err in collapse:
        P("  %-16s %13.2f  %s%s"
          % (meth, sp, "  ".join("%7.2f" % dd.get(n, float("nan"))
                                 for n in groups["dev_names"]), "  " + err if err else ""))
    P("")
    P("VERDICT PER SCAN")
    P("  %-7s %-6s %-6s %9s %8s %8s %8s %8s  %s"
      % ("scan", "sky", "theta", "Rrs555", "700/675", "412/555", "NIRrat", "ang", "verdict"))
    for nm in names:
        d = diag[nm]
        P("  %-7s %-6s %6.1f %9.5f %8.2f %8.3f %8.2f %8.2f  %s"
          % (nm, d["sky"], d["tilt"], d["rrs555"], d["peak700"], d["blue_slope"],
             d["nir_ratio"], d["angle_main"], d["verdict"]))
    if attrib:
        P("")
        P("WHAT IS THE NON-GLINT DIFFERENCE MADE OF?  (%s vs the main group)"
          % ", ".join(notglint))
        for k, v in sorted(attrib.items(), key=lambda kv: kv[1]["angle"]):
            P("  %-40s angle %6.2f deg   removes %5.1f %%" % (k, v["angle"], v["removed"]))
        P("  A partial removal means the difference is not one clean thing.")
    P("")
    P("Ruddick et al. (2006) NIR similarity ratio for pure-water residual: %.3f"
      % SIMILARITY_780_870)

    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(outdir, "REPORT.txt"), "w") as fh:
        fh.write(txt + "\n")
    with open(os.path.join(outdir, "per_scan_diagnostics.csv"), "w", newline="") as fh:
        w_ = csv.DictWriter(fh, fieldnames=["scan"] + list(next(iter(diag.values()))))
        w_.writeheader()
        for nm in names:
            w_.writerow(dict(scan=nm, **diag[nm]))
    print("\nwrote %s/REPORT.txt" % outdir)
    print("wrote %s/per_scan_diagnostics.csv" % outdir)
    for p in ps:
        print("wrote %s" % p)


if __name__ == "__main__":
    main()
