"""Cross-station comparison figures -- R_rs and GIOP composition, every validated station.

    python make_global_comparison_figs.py Data_NatureSpec/2026_Aug_16

Writes into `by_location/COMPARISON/`: `all_stations_overplot.png` (every station's final
R_rs, one axis), `giop_cross_station_all.png` (a_dg/b_bp/M_phi bar chart, free vs max
freedom, every station), `giop_cross_station.png` (the same, LOC1+LOC2 only -- the two
stations without the LOC3 shallow-water caveat), and `loc1_loc2abc_overplot.png` (R_rs,
LOC1 and the LOC2 split only).

Written because these four figures existed in `by_location/COMPARISON/` with no script
that produced them -- built once, ad hoc, earlier in the session. That is a real gap: when
LOC1's FINAL_Rrs.csv changed (the 2026-08-17 glint correction), there was no way to refresh
them short of hand-editing PNGs. This is that script, reconstructed to match what the
originals showed (verified against GLOBAL_COMPARISON.md's own description of each), so it
is reusable for every station added or corrected from here on, not just this one fix.
"""

import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

#: (station-relative path, legend label, color, linestyle, is_loc3)
STATIONS = [
    ("LOC1_66.89718N_162.60290W/FLENS8_FOV08", "LOC1 (n=12)", "#1f77b4", "-", False),
    ("LOC2a_66.89677N_162.57953W/FLENS8_FOV08", "LOC2a (n=9, main)", "#2e7d32", "-", False),
    ("LOC2b_66.89677N_162.57953W/FLENS8_FOV08", "LOC2b (n=3, disturbed)", "#c0392b", "-", False),
    ("LOC3_66.89235N_162.59149W/FIBR15_FOV15", "LOC3 FIBR15 (n=3, 15deg)", "#8a6000", "--", True),
    ("LOC3_66.89235N_162.59149W/FIBR15_FOV15_murky", "LOC3 FIBR15 murky (n=2)", "#d4a017", "--", True),
    ("LOC3_66.89235N_162.59149W/FLENS8_FOV08", "LOC3 FLENS8 (n=4, 8deg)", "#6a3d9a", "--", True),
]


def _read_final_rrs(path):
    with open(path) as fh:
        rows = list(csv.reader(fh))
    hdr_row = next(i for i, r in enumerate(rows) if r and r[0] == "wavelength_nm")
    data = np.array([[float(x) for x in r] for r in rows[hdr_row + 1:]])
    return data[:, 0], data[:, 1], data[:, 2]   # wl, mean, shape_sd


def _read_giop_final(path):
    """Pull the 'free' and 'max free' mean_spectrum rows: M_phi, a_dg(443), b_bp(443)."""
    out = {}
    with open(path) as fh:
        for row in csv.reader(fh):
            if not row or row[0].startswith("#"):
                continue
            if row[0] == "config":
                continue
            if row[0] == "free" and row[1] == "mean_spectrum":
                out["free"] = dict(M_phi=float(row[2]), adg443=float(row[3]), bbp443=float(row[4]))
            if row[0] == "max free" and row[1].startswith("mean_spectrum"):
                out["max_free"] = dict(M_phi=float(row[2]), adg443=float(row[3]), bbp443=float(row[4]))
    return out


def fig_all_stations_overplot(byloc, outdir, stations, title, fname):
    fig, ax = plt.subplots(figsize=(9, 6))
    for rel, label, color, ls, is_loc3 in stations:
        path = os.path.join(byloc, rel, "analysis", "FINAL_Rrs.csv")
        if not os.path.exists(path):
            print("  ! missing, skipped:", path); continue
        wl, mean, sd = _read_final_rrs(path)
        ax.plot(wl, mean, color=color, ls=ls, lw=1.8, label=label)
        ax.fill_between(wl, mean - sd, mean + sd, color=color, alpha=0.15, lw=0)
    ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("$R_{rs}$  sr$^{-1}$")
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=8.5); ax.grid(alpha=0.25)
    fig.tight_layout()
    p = os.path.join(outdir, fname)
    fig.savefig(p, dpi=140); plt.close(fig)
    print("wrote", p)


def fig_giop_cross_station(byloc, outdir, stations, title, fname):
    labels, free_vals, maxfree_vals = [], {"M_phi": [], "adg443": [], "bbp443": []}, \
                                       {"M_phi": [], "adg443": [], "bbp443": []}
    shaded = []
    for rel, label, color, ls, is_loc3 in stations:
        path = os.path.join(byloc, rel, "analysis", "GIOP", "giop_FINAL.csv")
        if not os.path.exists(path):
            print("  ! missing, skipped:", path); continue
        g = _read_giop_final(path)
        if "free" not in g or "max_free" not in g:
            print("  ! incomplete GIOP table, skipped:", path); continue
        labels.append(label.split(" (")[0])
        shaded.append(is_loc3)
        for k in ("M_phi", "adg443", "bbp443"):
            free_vals[k].append(g["free"][k])
            maxfree_vals[k].append(g["max_free"][k])

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    panels = [("adg443", "$a_{dg}(443)$  m$^{-1}$"), ("bbp443", "$b_{bp}(443)$  m$^{-1}$"),
             ("M_phi", "$M_\\phi$")]
    x = np.arange(len(labels))
    for ax, (key, ylab) in zip(axes, panels):
        ax.bar(x - 0.19, free_vals[key], width=0.38, label="free (OC4-seeded)",
              color="#1f77b4", hatch="//", edgecolor="k", linewidth=0.4)
        ax.bar(x + 0.19, maxfree_vals[key], width=0.38, label="max freedom",
              color="#ff7f0e", hatch="\\\\", edgecolor="k", linewidth=0.4)
        for i, is_loc3 in enumerate(shaded):
            if is_loc3:
                ax.axvspan(i - 0.5, i + 0.5, color="grey", alpha=0.12, zorder=0)
        ax.set_xticks(x); ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
        ax.set_ylabel(ylab); ax.grid(alpha=0.2, axis="y")
    axes[0].legend(fontsize=8, loc="upper left")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = os.path.join(outdir, fname)
    fig.savefig(p, dpi=140); plt.close(fig)
    print("wrote", p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("day_folder")
    a = ap.parse_args()
    byloc = os.path.join(a.day_folder.rstrip("/"), "by_location")
    outdir = os.path.join(byloc, "COMPARISON")
    os.makedirs(outdir, exist_ok=True)

    fig_all_stations_overplot(
        byloc, outdir, STATIONS,
        "Every validated station, one field day\n"
        "dashed = LOC3 (shallow-water caveat, see LOC3_BOTTOM_CAVEAT.md); band = shape uncertainty",
        "all_stations_overplot.png")

    loc12_only = [s for s in STATIONS if not s[4]]
    fig_all_stations_overplot(
        byloc, outdir, loc12_only, "LOC1 and the LOC2 split, R_rs",
        "loc1_loc2abc_overplot.png")

    fig_giop_cross_station(
        byloc, outdir, STATIONS,
        "GIOP retrievals across every validated station\n"
        "shaded = LOC3, conditional on optical depth (LOC3_BOTTOM_CAVEAT.md)",
        "giop_cross_station_all.png")

    fig_giop_cross_station(
        byloc, outdir, loc12_only,
        "GIOP retrievals -- LOC1 and LOC2 only (no shallow-water caveat)",
        "giop_cross_station.png")


if __name__ == "__main__":
    main()
