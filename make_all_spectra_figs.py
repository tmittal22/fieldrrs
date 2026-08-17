"""Every individual water spectrum, full range, real amplitude -- not just the mean.

    python make_all_spectra_figs.py <location folder>          # one sub-case
    python make_all_spectra_figs.py --site LOC2                # a site summary
    python make_all_spectra_figs.py --site LOC3

Every figure in the existing pipeline that shows "all the data" does it as a small inset
inside a bigger figure making a different point: fig4's third panel (angle-matched
spectra, but next to two other panels about pairing), fig9's first panel (SHAPE only --
normalised at 555, amplitude removed, 420-750 nm only). Neither is a clear, dedicated,
full-range, real-amplitude view of every individual scan. This fills that gap.

Uses the SAME R_rs (angle-matched pairing, per-scan rho, whatever --glint that station
was processed with -- read from the station's own REPORT.txt) as the rest of the
pipeline, so nothing here can disagree with FINAL_Rrs.csv by construction.
"""

import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyse_location import match_by_angle
from fieldrrs.rrs import rho_at_angle, rrs_three_scan, view_zenith_from_tilt
from organize_by_location import survey

LO, HI = 400.0, 900.0


def _glint_of(loc):
    """Read back which --glint this station's FINAL_Rrs.csv was built with."""
    rpt = os.path.join(loc, "analysis", "REPORT.txt")
    if not os.path.exists(rpt):
        return "none"
    m = re.search(r"glint correction applied to the final product: (\S+)",
                 open(rpt).read())
    return m.group(1) if m else "none"


def _excluded_of(loc):
    """Read back which scans --exclude-water dropped from FINAL_Rrs.csv, if any.

    per_scan_rrs() deliberately re-derives R_rs for EVERY water scan present on disk,
    ignoring --exclude-water -- that is the point of "all the data, not just means".
    But a figure that silently blends a known-different population back in without
    marking it would misrepresent what FINAL_Rrs.csv actually used. Read the same
    REPORT.txt line analyse_location.py itself prints for this.
    """
    rpt = os.path.join(loc, "analysis", "REPORT.txt")
    if not os.path.exists(rpt):
        return set()
    m = re.search(r"EXCLUDED FROM THE FINAL PRODUCT.*?\n\s*(\S.*?)\s*-- see",
                 open(rpt).read())
    if not m:
        return set()
    return {s.strip() for s in m.group(1).split(",")}


def per_scan_rrs(loc, glint="none", panel_r=0.99):
    scans = survey(loc)
    sky = [x for x in scans if x["role"] == "sky"]
    water = sorted([x for x in scans if x["role"] == "water"], key=lambda x: x["n"])
    wl = np.array(scans[0]["spec"].wavelength)
    if not water or not sky:
        return wl, {}
    out = {}
    for w, sk, dm, _ in match_by_angle(water, sky):
        rho = rho_at_angle(view_zenith_from_tilt(w["spec"].tilt_y_deg))
        r = rrs_three_scan(wl, w["spec"].columns["rad_target"],
                           sk["spec"].columns["rad_target"],
                           w["spec"].columns["rad_ref"], panel_r, rho, glint)
        out[w["n"]] = np.array(r.rrs)
    return wl, out


def fig_all_spectra(loc, out=None, title=None):
    """One sub-case: every individual water spectrum, coloured and labelled, plus the
    scaled mean for reference. This is the figure this module exists to provide."""
    glint = _glint_of(loc)
    excl = _excluded_of(loc)
    wl, R = per_scan_rrs(loc, glint)
    if not R:
        print("no water/sky scans at %s" % loc)
        return None
    out = out or os.path.join(loc, "analysis")
    os.makedirs(out, exist_ok=True)
    m = (wl >= LO) & (wl <= HI)
    names = sorted(R)
    kept = [n for n in names if n not in excl]
    cmap = plt.get_cmap("turbo" if len(names) > 6 else "tab10")

    fig, ax = plt.subplots(figsize=(12, 6.5))
    for i, n in enumerate(names):
        if n in excl:
            ax.plot(wl[m], R[n][m], lw=1.8, ls=":", color="#999999",
                    label="%s (EXCLUDED, not in FINAL_Rrs)" % n, alpha=0.85, zorder=4)
        else:
            ax.plot(wl[m], R[n][m], lw=1.6,
                    color=cmap(i / max(1, len(names) - 1)) if len(names) > 6 else cmap(i),
                    label=n, alpha=0.9)
    mean = np.mean([R[n] for n in kept], axis=0)
    ax.plot(wl[m], mean[m], lw=3.2, color="k", ls="--",
            label="plain mean (n=%d kept)" % len(kept), zorder=5)
    ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("$R_{rs}$  sr$^{-1}$")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.5, ncol=2 if len(names) > 8 else 1, loc="upper right")
    excl_note = ("  (%d shown dotted grey = excluded from FINAL_Rrs, see REPORT.txt)"
                 % len(excl) if excl else "")
    ax.set_title((title or os.path.basename(os.path.dirname(loc.rstrip("/"))) + "  ·  "
                 + os.path.basename(loc.rstrip("/"))) +
                "\nEVERY individual water scan on disk (n=%d), angle-matched pairing, "
                "glint=%s -- not just the mean%s"
                % (len(names), glint, excl_note), fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    p = os.path.join(out, "fig14_all_spectra.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    print("wrote", p, " n=%d scans (%d excluded): %s"
          % (len(names), len(excl), ", ".join(names)))
    return p


#: (location relative path, short label, colour) per site, so the site-summary figure
#: has one entry point to maintain instead of being re-derived from folder scanning.
SITES = {
    "LOC2": [
        ("Data_NatureSpec/2026_Aug_16/by_location/LOC2a_66.89677N_162.57953W/FLENS8_FOV08",
         "LOC2a (main)", "#2e7d32"),
        ("Data_NatureSpec/2026_Aug_16/by_location/LOC2b_66.89677N_162.57953W/FLENS8_FOV08",
         "LOC2b (disturbed)", "#c0392b"),
    ],
    "LOC3": [
        ("Data_NatureSpec/2026_Aug_16/by_location/LOC3_66.89235N_162.59149W/FIBR15_FOV15",
         "FIBR15 main (15°)", "#8a6000"),
        ("Data_NatureSpec/2026_Aug_16/by_location/LOC3_66.89235N_162.59149W/FIBR15_FOV15_murky",
         "FIBR15 murky (15°)", "#d4a017"),
        ("Data_NatureSpec/2026_Aug_16/by_location/LOC3_66.89235N_162.59149W/FLENS8_FOV08",
         "FLENS8 (8°)", "#6a3d9a"),
    ],
}


def fig_site_summary(site):
    """Every sub-case's every individual scan, on one page, grouped by sub-case."""
    subs = SITES[site]
    ncol = len(subs)
    fig, axes = plt.subplots(1, ncol, figsize=(6.3 * ncol, 6), sharey=False)
    axes = np.atleast_1d(axes)
    all_R = {}
    for ax, (loc, lab, col) in zip(axes, subs):
        glint = _glint_of(loc)
        excl = _excluded_of(loc)
        wl, R = per_scan_rrs(loc, glint)
        if not R:
            ax.axis("off"); ax.set_title("%s -- no data" % lab); continue
        m = (wl >= LO) & (wl <= HI)
        names = sorted(R)
        kept = [n for n in names if n not in excl]
        for i, n in enumerate(names):
            if n in excl:
                ax.plot(wl[m], R[n][m], lw=1.2, ls=":", color="#999999", alpha=0.8)
            else:
                ax.plot(wl[m], R[n][m], lw=1.3, color=col,
                        alpha=0.35 + 0.5 * (i + 1) / len(names))
        mean = np.mean([R[n] for n in kept], axis=0)
        ax.plot(wl[m], mean[m], lw=2.8, color="k", ls="--")
        ax.set_xlabel("wavelength (nm)")
        ax.grid(alpha=0.25)
        ax.set_title("%s\nn=%d kept%s, glint=%s"
                     % (lab, len(kept), " (+%d excluded, dotted)" % len(excl)
                        if excl else "", glint), fontsize=10.5)
        all_R[lab] = (wl, mean, col)
    axes[0].set_ylabel("$R_{rs}$  sr$^{-1}$")
    fig.suptitle("%s -- every sub-case, every individual scan" % site, fontsize=13.5,
                weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    outdir = "Data_NatureSpec/2026_Aug_16/by_location/COMPARISON"
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, "%s_site_summary.png" % site)
    fig.savefig(p, dpi=140); plt.close(fig)
    print("wrote", p)

    # second panel: all sub-case MEANS overlaid, so the site tells one story too
    fig, ax = plt.subplots(figsize=(9, 6))
    for lab, (wl, mean, col) in all_R.items():
        m = (wl >= LO) & (wl <= HI)
        ax.plot(wl[m], mean[m], lw=2.4, color=col, label=lab)
    ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("$R_{rs}$  sr$^{-1}$")
    ax.legend(fontsize=9); ax.grid(alpha=0.25)
    ax.set_title("%s -- sub-case means, overlaid" % site, fontsize=12.5)
    fig.tight_layout()
    p2 = os.path.join(outdir, "%s_site_means_overlay.png" % site)
    fig.savefig(p2, dpi=140); plt.close(fig)
    print("wrote", p2)
    return p, p2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("location", nargs="?", help="one sub-case folder")
    ap.add_argument("--site", choices=list(SITES), help="site summary instead")
    a = ap.parse_args()
    if a.site:
        fig_site_summary(a.site)
    elif a.location:
        fig_all_spectra(a.location.rstrip("/"))
    else:
        ap.error("give a location folder, or --site LOC2/LOC3")


if __name__ == "__main__":
    main()
