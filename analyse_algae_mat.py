"""LOC2c — the floating algae mat, a target the water R_rs pipeline does not fit.

    python analyse_algae_mat.py <by_location/LOC2c_*/FLENS8_FOV08>

00043 and 00044 were classified "land" by `organize_by_location.classify()` on the basis
of their spectral shape (a red-edge step that open water does not show). The photographs
(`analysis/water_scans/w1_contact_sheet.png` at LOC2a, or this station's own .jpg files)
show what that classification actually found: a floating clump of algae, one mid-water
and one against a dock, NOT dry land. "Land" here means "opaque enough that the sky is
illumination rather than a glint contaminant to subtract" -- the physical test the
pipeline actually applies -- not a claim about what the target is.

WHY THIS IS ITS OWN STATION, NOT PART OF LOC2a's water pipeline: `rrs_three_scan`
subtracts rho*L_sky because for a water column that term is reflected sky sitting on top
of the water-leaving signal. For a mat that is opaque enough to visually block the water
beneath it, subtracting rho*L_sky would remove real reflected light the mat is supposed
to have, not a contaminant. The right product for an opaque/emergent target is
reflectance, not R_rs (see `process_field_day.land_reflectance`, and LOC1's own
`fig7_land_targets.png` for the same target-type distinction there).

CAVEAT, stated once here rather than left implicit: a wet algae mat is not perfectly
opaque. Some fraction of the signal may be transmitted water-leaving radiance from
beneath or beside the mat, mixed into the "target" beam by the instrument's finite field
of view. The reflectance product below is the right FIRST-ORDER treatment (same one
LOC1 used for its own algae-on-concrete target), not a claim that the mat is opaque to
the metre.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from organize_by_location import survey
from process_field_day import land_reflectance


def clip(wl, v, lo=380.0, hi=950.0):
    p = [(w, x) for w, x in zip(wl, v) if lo <= w <= hi]
    return [a for a, _ in p], [b for _, b in p]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--panel-reflectance", type=float, default=0.99)
    a = ap.parse_args()

    scans = survey(a.folder)
    loc = os.path.basename(os.path.dirname(a.folder.rstrip("/")))
    fo = os.path.basename(a.folder.rstrip("/"))
    tag = "%s  ·  %s" % (loc, fo)
    outdir = os.path.join(a.folder, "analysis")
    os.makedirs(outdir, exist_ok=True)

    targets = sorted([s for s in scans if s["role"] == "land"], key=lambda x: x["n"])
    if not targets:
        print("no 'land'-role scans found at %s" % a.folder)
        return

    def band_mean(w, r, lo, hi):
        v = [x for ww, x in zip(w, r) if lo <= ww <= hi]
        return sum(v) / len(v) if v else float("nan")

    R, edge = {}, {}
    for s in targets:
        w, r = land_reflectance(s["spec"], a.panel_reflectance)
        w, r = np.array(w), np.array(r)
        R[s["n"]] = (w, r)
        # Self-computed from the REFLECTANCE this figure plots, same band convention
        # (745-755 / 665-675) as process_field_day.classify()'s diagnostic -- but that
        # diagnostic is a RAW TARGET RADIANCE ratio, used only to gate the sky/water/
        # land classification, not a property of R. Reusing its printed value here
        # would label this reflectance curve with a number computed from a different
        # quantity (radiance vs reflectance differ by the illumination spectrum's own
        # red/blue slope). Recomputed from R itself so the label matches the curve.
        edge[s["n"]] = band_mean(w, r, 745, 755) / band_mean(w, r, 665, 675)

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.4))
    ax = axes[0]
    for s in targets:
        w, r = R[s["n"]]
        wc, rc = clip(w, r)
        ax.plot(wc, rc, lw=2.2, label="%s  red-edge %.2fx" % (s["n"], edge[s["n"]]))
    ax.axvspan(670, 680, color="#c0392b", alpha=0.12)
    ax.text(682, ax.get_ylim()[1] * 0.55, "chlorophyll\nabsorption", fontsize=8.5,
            color="#c0392b")
    ax.axvspan(700, 760, color="#2e7d32", alpha=0.10)
    ax.text(706, ax.get_ylim()[1] * 0.12, "red edge", fontsize=8.5, color="#2e7d32")
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("reflectance factor  $R$  (NOT $R_{rs}$)")
    ax.legend(fontsize=9); ax.grid(alpha=0.25)
    ax.set_title("$R = R_{panel}\\,L_t/L_{panel}$ -- no $\\rho$, no sky subtraction.\n"
                 "A floating mat, treated as opaque, like LOC1's land targets.",
                 fontsize=10, loc="left")

    ax = axes[1]
    ax.axis("off")
    lines = ["WHAT THESE SCANS ARE", ""]
    for s in targets:
        w, r = R[s["n"]]
        i550, i670, i750 = (np.argmin(abs(w - x)) for x in (550, 670, 750))
        lines.append("%s  R(550)=%.4f  R(670)=%.4f  R(750)=%.4f" %
                     (s["n"], r[i550], r[i670], r[i750]))
        lines.append("   red edge (from R, 745-755/665-675) = %.2f  ->  %s"
                     % (edge[s["n"]],
                        "VEGETATED signature (chlorophyll)" if edge[s["n"]] > 3 else
                        "a real, visible red-edge rise, below the >3x threshold "
                        "classify() uses for its OWN raw-radiance diagnostic"))
        lines.append("")
    lines += ["Both scans classified 'land' by SHAPE (an 850-880/450-650 nm radiance",
             "ratio open water does not show -- see process_field_day.classify()),",
             "not by the field notes. The photographs (see the folder's .jpg files,",
             "or LOC2a's contact sheet) show floating algae clumps, one mid-water and",
             "one against a dock -- not dry ground. 'Land' in this pipeline means",
             "'opaque enough that sky is illumination, not glint to remove', which is",
             "a physical test this figure applies correctly even though the English",
             "word is misleading here.",
             "",
             "NOTE ON THE RED-EDGE NUMBER: classify()'s own 'red_edge' diagnostic is a",
             "RAW TARGET RADIANCE ratio (no panel normalisation), used only to print a",
             "vegetated/bare label in reports -- it does not gate the classification",
             "itself (nir_vis and blue_green do). It is NOT a property of the",
             "reflectance R plotted here (radiance and reflectance red-edge ratios",
             "differ by the illumination spectrum's own red/blue slope). The number",
             "above is recomputed FROM R with the same band convention, so it matches",
             "the curve it is printed next to.",
             "",
             "This is REFLECTANCE, not R_rs -- do not compare its magnitude directly",
             "to a water R_rs spectrum without accounting for the different physical",
             "quantity (see the module docstring)."]
    ax.text(0.0, 0.98, "\n".join(lines), transform=ax.transAxes, va="top",
            fontsize=9, family="monospace")

    fig.suptitle("LOC2c — THE FLOATING ALGAE MAT, %s" % tag, fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    p = os.path.join(outdir, "c1_algae_reflectance.png")
    fig.savefig(p, dpi=140); plt.close(fig)

    with open(os.path.join(outdir, "REPORT.txt"), "w") as fh:
        fh.write("LOC2c -- floating algae mat, %s\n" % tag)
        fh.write("%d scans, classified 'land' by spectral shape (see module "
                 "docstring for what that means here).\n\n" % len(targets))
        for s in targets:
            w, r = R[s["n"]]
            i550, i670, i750 = (np.argmin(abs(w - x)) for x in (550, 670, 750))
            fh.write("%s  R(550)=%.4f  R(670)=%.4f  R(750)=%.4f  red-edge=%.2fx "
                     "(from R, not classify()'s raw-radiance diagnostic)\n"
                     % (s["n"], r[i550], r[i670], r[i750], edge[s["n"]]))
        fh.write("\nReflectance, NOT R_rs. See module docstring for why an opaque-\n"
                 "target treatment was used and its limits.\n")
    print("wrote %s" % p)
    print("wrote %s/REPORT.txt" % outdir)
    return R, targets


if __name__ == "__main__":
    main()
