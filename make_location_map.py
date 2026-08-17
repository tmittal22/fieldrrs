"""Station map on satellite imagery.

    python make_location_map.py Data_NatureSpec/2026_Aug_16 --out <dir> --place "Kotzebue, Alaska"

Imagery is fetched once and cached under .tilecache/, so re-runs are offline.

Everything on the map is derived from the scans themselves (date, coordinates, the
context-inset bounding box) EXCEPT the place name, which nothing in a .sed file records
-- pass `--place` or the title/context annotation just omit it rather than guess or
carry over a stale placeholder from a previous day's run.
"""

import argparse
import datetime
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import basemap as bm
from organize_by_location import cluster, fov_deg, survey

COL = {"water": "#00d5ff", "sky": "#ffd166", "land": "#7CFC00"}
SITE = ["#ff2d55", "#ffffff", "#ff9f0a"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--out", default=".")
    ap.add_argument("--zoom", type=int, default=16)
    ap.add_argument("--provider", default="esri_imagery")
    ap.add_argument("--place", default="",
                    help="place name for the title/inset, e.g. 'Kotzebue, Alaska'. "
                         "Nothing in a .sed file records this, so it is not guessed -- "
                         "omit it and the map just shows date + coordinates.")
    a = ap.parse_args()

    scans = survey(a.folder)
    locs = cluster(scans)
    # ground truth, not the folder name (which can be named anything) or a hardcoded
    # date left over from whichever day this script was last run on.
    dates = sorted({s["spec"].header.get("Date", "").split(",")[0].strip()
                    for s in scans if s["spec"].header.get("Date")})
    if dates:
        try:
            mm, dd, yy = dates[0].split("/")
            date_str = datetime.date(int(yy), int(mm), int(dd)).isoformat()
        except ValueError:
            date_str = dates[0]
        if len(dates) > 1:
            date_str += " to " + dates[-1]
    else:
        date_str = "date unknown"
    lats = [s["lat"] for s in scans]
    lons = [s["lon"] for s in scans]
    padx = max(0.006, (max(lons) - min(lons)) * 0.45)
    pady = max(0.0016, (max(lats) - min(lats)) * 0.45)
    box = (min(lats) - pady, max(lats) + pady, min(lons) - padx, max(lons) + padx)

    img, ext = bm.mosaic(*box, zoom=a.zoom, provider=a.provider)
    err = bm.check_distortion(ext, img.size[1])

    fig = plt.figure(figsize=(15.5, 8.2))
    ax = fig.add_axes([0.05, 0.08, 0.60, 0.84])
    ax.imshow(img, extent=ext, origin="upper", interpolation="bilinear")
    ax.set_xlim(box[2], box[3])
    ax.set_ylim(box[0], box[1])
    latm = sum(lats) / len(lats)
    ax.set_aspect(1.0 / math.cos(math.radians(latm)))

    for s in scans:
        ax.scatter(s["lon"], s["lat"], s=30, color=COL[s["role"]], alpha=0.9,
                   edgecolor="k", linewidth=0.4, zorder=4)

    # Label boxes are positioned in AXES FRACTION, not data coordinates, so they can
    # never leave the canvas however the sites happen to fall. Leader lines still point
    # at the true positions.
    place = [(0.055, 0.93), (0.055, 0.75), (0.055, 0.57)]
    for i, c in enumerate(locs):
        gps = [x["gps"] for x in c["scans"] if x["gps"] is not None]
        fos = sorted({x["fo"] for x in c["scans"]})
        n_w = sum(1 for x in c["scans"] if x["role"] == "water")
        n_s = sum(1 for x in c["scans"] if x["role"] == "sky")
        ax.scatter(c["lon"], c["lat"], s=420, marker="o", facecolor="none",
                   edgecolor=SITE[i % 3], linewidth=2.8, zorder=6)
        ax.annotate(
            "LOC%d   %.5f$^\\circ$N  %.5f$^\\circ$W\n"
            "%d scans  ·  %d water / %d sky\n"
            "%s UTC  ·  %s"
            % (i + 1, c["lat"], -c["lon"], len(c["scans"]), n_w, n_s,
               "%02d:%02d" % (int(min(gps)), int(round((min(gps) % 1) * 60)))
               if gps else "?",
               " + ".join("%s (%.0f$^\\circ$ FOV)" % (f, fov_deg(f)) for f in fos)),
            xy=(c["lon"], c["lat"]), xycoords="data",
            xytext=place[i % len(place)], textcoords="axes fraction",
            fontsize=9.5, weight="bold", color="k", zorder=7, ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=SITE[i % 3], lw=2,
                      alpha=0.93),
            arrowprops=dict(arrowstyle="-", color=SITE[i % 3], lw=2))

    # scale bar
    kx = 111320.0 * math.cos(math.radians(latm))
    d = 250.0 / kx
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    sx = x1 - 0.05 * (x1 - x0) - d
    sy = y0 + 0.06 * (y1 - y0)
    ax.plot([sx, sx + d], [sy, sy], color="w", lw=6, solid_capstyle="butt", zorder=8)
    ax.plot([sx, sx + d], [sy, sy], color="k", lw=3, solid_capstyle="butt", zorder=9)
    ax.text(sx + d / 2, sy + 0.012 * (y1 - y0), "250 m", ha="center", fontsize=10,
            weight="bold", color="w", zorder=9,
            path_effects=None)
    ax.set_xlabel("longitude ($^\\circ$E)"); ax.set_ylabel("latitude ($^\\circ$N)")
    ax.set_title("%s%s%d scans at %d locations"
                 % (date_str, "  ·  %s  ·  " % a.place if a.place else "  ·  ",
                    len(scans), len(locs)), fontsize=13, weight="bold")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=COL[k],
                              markeredgecolor="k", label=k)
                       for k in ("water", "sky", "land")],
              fontsize=9.5, loc="lower left", framealpha=0.92)
    # matplotlib renders these longitudes with a '+6.689e1' offset by default, which is
    # unreadable on a map; force plain ticks.
    for axis in (ax.xaxis, ax.yaxis):
        axis.get_major_formatter().set_useOffset(False)
        axis.get_major_formatter().set_scientific(False)
    ax.text(0.0, -0.115, bm.ATTRIBUTION[a.provider] + "   |   Mercator-on-linear axis "
            "error %.2f px" % err, transform=ax.transAxes, fontsize=7.5, color="#444")

    # ---- right column: regional context + per-location detail
    # Computed from the data's own centroid, not hardcoded to any one region -- a fixed
    # Alaska/Kotzebue-Sound box here would silently mislabel any other field day's map.
    axc = fig.add_axes([0.685, 0.56, 0.29, 0.36])
    lonm = sum(lons) / len(lons)
    ctx_lat = (latm - 8, latm + 8)
    ctx_lon = (lonm - 16, lonm + 16)
    axc.set_xlim(*ctx_lon); axc.set_ylim(*ctx_lat)
    axc.add_patch(Rectangle((box[2], box[0]), box[3] - box[2], box[1] - box[0],
                            fc="none", ec="#ff2d55", lw=2, zorder=3))
    axc.scatter([lonm], [latm], s=150, marker="*", color="#ff2d55", zorder=4)
    if ctx_lat[0] <= 66.5622 <= ctx_lat[1]:
        axc.axhline(66.5622, color="#2c6f9b", ls="--", lw=1.3)
        axc.text(0.02, 0.9, "Arctic Circle", fontsize=8.5, color="#2c6f9b",
                 transform=axc.transAxes)
    if a.place:
        axc.annotate(a.place, (lonm, latm), xytext=(0.55, 0.12),
                     textcoords="axes fraction",
                     fontsize=10, weight="bold", color="#ff2d55",
                     arrowprops=dict(arrowstyle="->", color="#ff2d55", lw=1.5))
    axc.set_title("Regional context" + (" — %s" % a.place if a.place else ""),
                 fontsize=10.5, loc="left")
    axc.grid(alpha=0.3); axc.tick_params(labelsize=8)

    axt = fig.add_axes([0.685, 0.08, 0.29, 0.42]); axt.axis("off")
    lines = ["OCCUPATIONS, in time order", ""]
    for i, c in enumerate(locs):
        gps = [x["gps"] for x in c["scans"] if x["gps"] is not None]
        lines.append("LOC%d   %.5f N  %.5f W" % (i + 1, c["lat"], -c["lon"]))
        lines.append("   %s-%s UTC   sun %.1f-%.1f deg"
                     % ("%02d:%02d" % (int(min(gps)), int(round((min(gps) % 1) * 60))),
                        "%02d:%02d" % (int(max(gps)), int(round((max(gps) % 1) * 60))),
                        min(x["sun"] for x in c["scans"]),
                        max(x["sun"] for x in c["scans"])))
        by = {}
        for s in c["scans"]:
            by.setdefault(s["fo"], []).append(s)
        for fo in sorted(by):
            g = by[fo]
            lines.append("   %-7s FOV %2.0f deg  %2d scans: %d water, %d sky, %d land"
                         % (fo, fov_deg(fo), len(g),
                            sum(1 for x in g if x["role"] == "water"),
                            sum(1 for x in g if x["role"] == "sky"),
                            sum(1 for x in g if x["role"] == "land")))
        lines.append("")
    mixed = [i + 1 for i, c in enumerate(locs)
             if len({x["fo"] for x in c["scans"]}) > 1]
    if mixed:
        lines.append("LOC%s uses TWO foreoptics." % ", ".join(map(str, mixed)))
        lines.append("Treat each FOV as a separate dataset:")
        lines.append("8 and 15 deg average over different")
        lines.append("footprints, so they do not share rho.")
    axt.text(0, 1, "\n".join(lines), va="top", ha="left", fontsize=9,
             family="monospace", transform=axt.transAxes)

    os.makedirs(a.out, exist_ok=True)
    p = os.path.join(a.out, "fig1_map.png")
    fig.savefig(p, dpi=145)
    plt.close(fig)
    print("wrote %s" % p)
    print("distortion %.2f px (negligible below 1)" % err)


if __name__ == "__main__":
    main()
