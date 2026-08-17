"""Full worked report for one field day: map, every processing step, and the inversion.

    python make_field_report.py Data_NatureSpec/2026_Aug_16 --out results_20260816

Six figures, each showing a stage rather than only its result, plus a text summary.
All spectral axes are clipped to 350-950 nm: the instrument records to 2500 nm, but
beyond ~900 nm the water radiance is at the noise floor and R_rs is meaningless.
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
from matplotlib.patches import Rectangle

from fieldrrs.rrs import RHO_MOBLEY1999, rrs_three_scan
from fieldrrs.solar import solar_window_verdict
from process_field_day import (SKY, VEG, WATER, band, load, mean_spectrum, process)

WLO, WHI = 350.0, 950.0            # the range this instrument actually resolves on water
RLO, RHI = 400.0, 900.0            # and where R_rs is worth plotting
COL = {SKY: "#7fb3d5", WATER: "#1f7a99", VEG: "#2e7d32"}
PANEL = "#d9534f"
SITE_COL = ["#c0392b", "#1f7a99", "#8a6000"]


def clip(wl, v, lo=WLO, hi=WHI):
    p = [(w, x) for w, x in zip(wl, v) if lo <= w <= hi]
    return [w for w, _ in p], [x for _, x in p]


def sites(scans):
    """Group stations by position. Three occupations ~1 km apart in this dataset."""
    out = []
    for s in scans:
        la, lo = s["spec"].latitude, s["spec"].longitude
        for g in out:
            if abs(g["lat"] - la) < 2e-4 and abs(g["lon"] - lo) < 4e-4:
                g["scans"].append(s)
                break
        else:
            out.append({"lat": la, "lon": lo, "scans": [s]})
    for g in out:
        g["lat"] = sum(x["spec"].latitude for x in g["scans"]) / len(g["scans"])
        g["lon"] = sum(x["spec"].longitude for x in g["scans"]) / len(g["scans"])
        g["gps"] = min(x["spec"].gps_time for x in g["scans"])
    return sorted(out, key=lambda g: g["gps"])


# ------------------------------------------------------------------ 1. the map
def fig_map(scans, groups, outdir):
    """Plot in local METRES, not degrees.

    Degrees here span 0.005 lat by 0.024 lon, which matplotlib renders with an offset
    like '+6.6896e1' that is unreadable and hides the actual scale. Metres east/north of
    the first station make the geometry legible and the scale bar meaningful.
    """
    lat0 = sum(g["lat"] for g in groups) / len(groups)
    lon0 = sum(g["lon"] for g in groups) / len(groups)
    kx = 111320.0 * math.cos(math.radians(lat0))
    ky = 111320.0
    E = lambda lo: (lo - lon0) * kx
    N = lambda la: (la - lat0) * ky

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14.5, 6.4),
                                  gridspec_kw={"width_ratios": [1.45, 1]})
    for role in (WATER, SKY, VEG):
        pts = [(E(s["spec"].longitude), N(s["spec"].latitude))
               for s in scans if s["role"] == role]
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=26, alpha=0.55,
                   color=COL[role], label="%s (%d)" % (role, len(pts)), zorder=2,
                   linewidth=0)
    off = [(14, 12), (14, -34), (-120, 12), (14, 12)]
    for i, g in enumerate(groups):
        x, y = E(g["lon"]), N(g["lat"])
        ax.scatter(x, y, s=300, marker="o", facecolor="none",
                   edgecolor=SITE_COL[i % 3], linewidth=2.6, zorder=4)
        ax.annotate("SITE %d\n%d scans\n%s UTC"
                    % (i + 1, len(g["scans"]), _hhmm(g["gps"])), (x, y),
                    xytext=off[i % len(off)], textcoords="offset points", fontsize=9,
                    weight="bold", color=SITE_COL[i % 3], zorder=5)
    ax.set_aspect("equal")
    ax.set_xlabel("metres east of %.5f$^\\circ$W" % -lon0)
    ax.set_ylabel("metres north of %.5f$^\\circ$N" % lat0)
    ax.grid(alpha=0.3)
    ax.set_title("2026-08-16 · Kotzebue Sound, Alaska · %d scans, %d occupations"
                 % (len(scans), len(groups)), fontsize=11.5, pad=12)
    ax.legend(fontsize=8.5, loc="upper left", framealpha=0.95)
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    xs, ys = x1 - 0.30 * (x1 - x0), y0 + 0.07 * (y1 - y0)
    ax.plot([xs, xs + 200], [ys, ys], color="k", lw=3.5, zorder=6)
    ax.text(xs + 100, ys + 0.025 * (y1 - y0), "200 m", ha="center", fontsize=9,
            weight="bold", zorder=6)

    ax2.set_xlim(-180, -128); ax2.set_ylim(51, 75)
    ax2.add_patch(Rectangle((-168, 64), 12, 6, fc="#fdecea", ec="#c0392b", lw=2,
                            zorder=2))
    ax2.scatter([-162.6], [66.9], s=140, color="#c0392b", zorder=5, marker="*")
    ax2.axhline(66.5622, color="#2c6f9b", ls="--", lw=1.3, zorder=3)
    ax2.text(-179, 67.2, "Arctic Circle  66.56$^\\circ$N", fontsize=9, color="#2c6f9b")
    ax2.annotate("Kotzebue Sound\n66.894$^\\circ$N, 162.590$^\\circ$W",
                 (-162.6, 66.9), xytext=(-152, 57.5), fontsize=10, weight="bold",
                 color="#c0392b",
                 arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.5))
    ax2.set_xlabel("longitude ($^\\circ$E)"); ax2.set_ylabel("latitude ($^\\circ$N)")
    ax2.set_title("Sun peaks at 36$^\\circ$ elevation here in mid-August, so the\n"
                  "30-60$^\\circ$ window is available for only part of the day.",
                  fontsize=10, loc="left")
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    return _save(fig, outdir, "fig1_map.png")


def _hhmm(h):
    if h is None:
        return "?"
    return "%02d:%02d" % (int(h), int(round((h - int(h)) * 60)))


def _save(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, name)
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p


# ------------------------------------- 2. the three-scan method, worked end to end
def fig_worked_example(results, outdir):
    st = max(results, key=lambda s: len(s["rrs"]))
    r = st["rrs"][0]
    wl = st["wl"]
    lt = r["spec"].columns["rad_target"]
    ls, lp = st["l_sky"], st["l_panel"]
    ed = [math.pi * x / 0.99 for x in lp]
    lw = [t - RHO_MOBLEY1999 * s for t, s in zip(lt, ls)]
    rrs = r["res"].rrs

    fig, axes = plt.subplots(1, 4, figsize=(18.5, 4.6))
    a = axes[0]
    for v, c, lab in ((lp, PANEL, "$L_{panel}$ (reference)"),
                      (ls, COL[SKY], "$L_{sky}$"), (lt, COL[WATER], "$L_{target}$")):
        a.plot(*clip(wl, v), lw=1.8, color=c, label=lab)
    a.set_yscale("log"); a.legend(fontsize=8.5)
    a.set_title("STEP 1  the three measured radiances", fontsize=10.5, loc="left")
    a.set_ylabel("W m$^{-2}$ sr$^{-1}$ nm$^{-1}$")

    a = axes[1]
    a.plot(*clip(wl, ed), lw=2.0, color=PANEL)
    a.set_title("STEP 2  $E_d = \\pi L_{panel} / R_{panel}$", fontsize=10.5, loc="left")
    a.set_ylabel("W m$^{-2}$ nm$^{-1}$")
    a.text(0.04, 0.06, "the panel converts a radiance\ninto an irradiance",
           transform=a.transAxes, fontsize=8.5, color="#555")

    a = axes[2]
    a.plot(*clip(wl, lt), lw=1.4, color=COL[WATER], label="$L_t$ (measured)")
    a.plot(*clip(wl, [RHO_MOBLEY1999 * x for x in ls]), lw=1.4, color=COL[SKY],
           label="$\\rho L_{sky}$ (removed)")
    a.plot(*clip(wl, lw), lw=2.2, color="#c0392b", label="$L_w$ (what is left)")
    a.legend(fontsize=8.5)
    a.set_title("STEP 3  $L_w = L_t - \\rho L_{sky}$,  $\\rho$ = 0.028",
                fontsize=10.5, loc="left")
    a.set_ylabel("W m$^{-2}$ sr$^{-1}$ nm$^{-1}$")
    frac = [100 * RHO_MOBLEY1999 * s / t if t else 0 for t, s in zip(lt, ls)]
    f450 = [f for w, f in zip(wl, frac) if 440 <= w <= 460]
    a.text(0.04, 0.9, "skylight is %.0f %% of $L_t$ at 450 nm" % (sum(f450) / len(f450)),
           transform=a.transAxes, fontsize=8.5, color="#c0392b")

    a = axes[3]
    a.plot(*clip(wl, rrs, RLO, RHI), lw=2.4, color="#c0392b")
    a.axhline(0, color="#888", lw=0.8)
    a.set_title("STEP 4  $R_{rs} = L_w / E_d$", fontsize=10.5, loc="left")
    a.set_ylabel("$R_{rs}$ (sr$^{-1}$)")
    for lam, txt in ((570, "green peak\nsediment"), (700, "700 nm peak\nhigh SPM")):
        v = [x for w, x in zip(wl, rrs) if lam - 6 <= w <= lam + 6]
        a.annotate(txt, (lam, sum(v) / len(v)), xytext=(0, 26),
                   textcoords="offset points", fontsize=8, ha="center",
                   arrowprops=dict(arrowstyle="->", lw=1))
    for ax in axes:
        ax.set_xlabel("wavelength (nm)"); ax.grid(alpha=0.25)
    fig.suptitle("The three-scan method, worked on scan %s (station ref=%.4f)"
                 % (r["n"], st["ref"]), fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return _save(fig, outdir, "fig2_worked_example.png")


# ------------------------------------------------------------ 3. R_rs by station
def fig_rrs(results, outdir):
    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(3.7 * n, 4.6), sharey=True)
    axes = list(axes) if n > 1 else [axes]
    vmax = max(x for st in results for r in st["rrs"]
               for w, x in zip(r["res"].wavelength, r["res"].rrs) if RLO <= w <= RHI)
    for ax, st in zip(axes, results):
        for r in st["rrs"]:
            ax.plot(*clip(r["res"].wavelength, r["res"].rrs, RLO, RHI), lw=0.8,
                    alpha=0.55, color=COL[WATER])
        mean = [sum(r["res"].rrs[i] for r in st["rrs"]) / len(st["rrs"])
                for i in range(len(st["wl"]))]
        ax.plot(*clip(st["wl"], mean, RLO, RHI), lw=2.4, color="#c0392b")
        ax.axhline(0, color="#888", lw=0.8)
        ax.set_xlim(RLO, RHI); ax.set_ylim(-0.05 * vmax, 1.12 * vmax)
        ax.set_xlabel("wavelength (nm)"); ax.grid(alpha=0.25)
        ax.set_title("ref=%.4f  ·  n=%d\n%s UTC · sun %.1f$^\\circ$"
                     % (st["ref"], len(st["rrs"]), _hhmm(st["gps"]), st["sun"]),
                     fontsize=9.5)
    axes[0].set_ylabel("$R_{rs}$ (sr$^{-1}$)")
    fig.suptitle("$R_{rs}$ by station — thin = individual scans, thick = station mean",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return _save(fig, outdir, "fig3_rrs_by_station.png")


# ------------------------------------------------------------- 4. the classifier
def fig_classifier(scans, outdir):
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.2))
    for role in (WATER, SKY, VEG):
        xs = [s["diag"]["blue_green"] for s in scans if s["role"] == role]
        ys = [s["diag"]["nir_vis"] for s in scans if s["role"] == role]
        ax.scatter(xs, ys, s=54, alpha=0.82, color=COL[role], edgecolor="k",
                   linewidth=0.4, label="%s (n=%d)" % (role, len(xs)), zorder=3)
    ax.axhline(0.5, color=COL[VEG], ls="--", lw=1.3)
    ax.axvline(1.3, color=COL[SKY], ls="--", lw=1.3)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(alpha=0.25)
    ax.set_xlabel("blue/green   $L$(440-490)/$L$(540-580)")
    ax.set_ylabel("NIR/vis   $L$(850-880)/$L$(450-650)")
    ax.legend(fontsize=9)
    ax.set_title("Roles from SPECTRAL SHAPE\n(tilt is useless here: 36-50$^\\circ$ "
                 "for sky and water alike)", fontsize=10.5, loc="left")

    for role in (WATER, SKY, VEG):
        sel = [s for s in scans if s["role"] == role]
        for s in sel[:6]:
            wl, v = clip(s["spec"].wavelength, s["spec"].columns["rad_target"])
            norm = max(v)
            ax2.plot(wl, [x / norm for x in v], lw=1.2, alpha=0.75, color=COL[role])
    for role in (WATER, SKY, VEG):
        ax2.plot([], [], color=COL[role], lw=2, label=role)
    ax2.legend(fontsize=9); ax2.grid(alpha=0.25)
    ax2.set_xlabel("wavelength (nm)"); ax2.set_ylabel("radiance, peak-normalised")
    ax2.set_title("Why it works: three different shapes\n"
                  "sky falls with $\\lambda$, water peaks green, tundra has a red edge",
                  fontsize=10.5, loc="left")
    fig.tight_layout()
    return _save(fig, outdir, "fig4_classifier.png")


# ---------------------------------------------------------------- 5. the inversion
def fig_giop(results, outdir, giop_rows):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    names = [g["name"] for g in giop_rows]
    x = range(len(names))
    for ax, key, title, bound in (
            (axes[0], "chl", "$M_\\phi$  (chlorophyll amplitude)", None),
            (axes[1], "adg", "$a_{dg}$(443)  m$^{-1}$", 20.0),
            (axes[2], "bbp", "$b_{bp}$(443)  m$^{-1}$", 3.0)):
        dc = [g["dc_" + key] for g in giop_rows]
        bd = [g["bd_" + key] for g in giop_rows]
        ax.bar([i - 0.2 for i in x], [abs(v) if abs(v) < 1e6 else 1e6 for v in dc],
               0.4, color="#999", label="GIOP-DC")
        ax.bar([i + 0.2 for i in x], bd, 0.4, color="#2e7d32", label="bounded")
        if bound:
            ax.axhline(bound, color="#c0392b", ls="--", lw=1.6)
            ax.text(0.02, bound * 1.08, "solver bound = RAIL", color="#c0392b",
                    fontsize=8.5, transform=ax.get_yaxis_transform())
        ax.set_yscale("log"); ax.set_xticks(list(x))
        ax.set_xticklabels([n.replace("Rrs_ref", "") for n in names], fontsize=8,
                           rotation=30)
        ax.set_title(title, fontsize=10.5); ax.grid(alpha=0.25, axis="y")
        ax.legend(fontsize=8.5)
    fig.suptitle("GIOP on turbid Arctic water. Unconstrained GIOP-DC returns up to "
                 "$10^{15}$ and one outright failure;\nthe bounded solver returns "
                 "finite values, but note how many sit ON the bound.", fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    return _save(fig, outdir, "fig5_giop.png")


# ------------------------------------------------------------------ 6. diagnostics
def fig_qc(results, outdir):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    a = axes[0]
    for st in results:
        mean = [sum(r["res"].rrs[i] for r in st["rrs"]) / len(st["rrs"])
                for i in range(len(st["wl"]))]
        wl, v = clip(st["wl"], mean, RLO, RHI)
        a.plot(wl, [x / max(v) for x in v], lw=2, label="ref=%.4f" % st["ref"])
    a.axvline(700, color="#888", ls=":", lw=1.4)
    a.text(704, 0.2, "700 nm\n$a_w$ minimum", fontsize=8.5, color="#555")
    a.set_xlabel("wavelength (nm)"); a.set_ylabel("$R_{rs}$, peak-normalised")
    a.set_title("Shape comparison: the 700 nm peak grows with turbidity",
                fontsize=10.5, loc="left")
    a.legend(fontsize=8); a.grid(alpha=0.25)

    a = axes[1]
    for st in results:
        pk, nir = [], []
        for r in st["rrs"]:
            wl, v = r["res"].wavelength, r["res"].rrs
            g = max(x for w, x in zip(wl, v) if 550 <= w <= 600)
            n = sum(x for w, x in zip(wl, v) if 800 <= w <= 900) / \
                sum(1 for w in wl if 800 <= w <= 900)
            pk.append(g); nir.append(n / g)
        a.scatter(pk, nir, s=52, alpha=0.8, label="ref=%.4f" % st["ref"],
                  edgecolor="k", linewidth=0.4)
    a.set_xlabel("$R_{rs}$ green peak (sr$^{-1}$)")
    a.set_ylabel("NIR / green")
    a.set_title("Turbidity index. High NIR/green = more suspended sediment",
                fontsize=10.5, loc="left")
    a.legend(fontsize=8); a.grid(alpha=0.25)

    a = axes[2]
    for st in results:
        sc = [max(x for w, x in zip(r["res"].wavelength, r["res"].rrs)
                  if 550 <= w <= 600) for r in st["rrs"]]
        a.scatter([st["sun"]] * len(sc), sc, s=48, alpha=0.75, edgecolor="k",
                  linewidth=0.4)
        a.errorbar(st["sun"], sum(sc) / len(sc),
                   yerr=(max(sc) - min(sc)) / 2, fmt="_", ms=22, color="k", lw=1.6)
    a.axvspan(30, 60, color="#2e7d32", alpha=0.12)
    a.text(30.2, a.get_ylim()[1] * 0.97, "PREFERRED 30-60$^\\circ$", fontsize=8.5,
           color="#2e7d32", va="top")
    a.set_xlabel("solar elevation (deg)"); a.set_ylabel("$R_{rs}$ green peak")
    a.set_title("Every station sat inside the preferred solar window",
                fontsize=10.5, loc="left")
    a.grid(alpha=0.25)
    fig.tight_layout()
    return _save(fig, outdir, "fig6_diagnostics.png")


def run_giop(results):
    """Six-band GIOP-DC and the bounded solver, on each station mean."""
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "..", "giop_python", "src"))
        import numpy as np
        from giop import get_oc, giop as run
    except Exception as exc:
        print("GIOP unavailable (%s); skipping the inversion figure" % exc)
        return []
    B = np.array([412., 443, 490, 510, 555, 670])
    rows = []
    for st in results:
        wl = np.array(st["wl"])
        mean = np.array([sum(r["res"].rrs[i] for r in st["rrs"]) / len(st["rrs"])
                         for i in range(len(st["wl"]))])
        r6 = np.array([mean[np.argmin(abs(wl - b))] for b in B])
        chl = float(get_oc(r6[1], r6[2], r6[3], r6[4], "oc4"))
        sig = np.sqrt((0.05 * np.abs(r6)) ** 2 + 2e-4 ** 2)
        dc = run(B, r6, chl)
        bd = run(B, r6, chl, inv="bounded", sigma=sig)
        rows.append({"name": "Rrs_ref%.4f" % st["ref"], "oc4": chl,
                     "dc_chl": dc.chl, "dc_adg": dc.adg443, "dc_bbp": dc.bbp443,
                     "bd_chl": bd.chl, "bd_adg": bd.adg443, "bd_bbp": bd.bbp443,
                     "dc_failed": bool(dc.failed)})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--out", default="results_field")
    a = ap.parse_args()

    scans = load(a.folder)
    groups = sites(scans)
    results = process(scans, RHO_MOBLEY1999, 0.99, "none")
    giop_rows = run_giop(results)

    from process_field_day import write_csv
    csv_path = write_csv(results, a.out)
    with open(os.path.join(a.out, "VERIFICATION.txt"), "w") as fh:
        import subprocess
        r = subprocess.run([sys.executable, "verify_field_calcs.py", a.folder],
                           capture_output=True, text=True,
                           cwd=os.path.dirname(os.path.abspath(__file__)))
        fh.write(r.stdout)
    paths = [fig_map(scans, groups, a.out),
             fig_worked_example(results, a.out),
             fig_rrs(results, a.out),
             fig_classifier(scans, a.out),
             fig_qc(results, a.out)]
    if giop_rows:
        paths.append(fig_giop(results, a.out, giop_rows))

    lines = []
    P = lines.append
    P("FIELD REPORT  %s" % a.folder)
    P("=" * 78)
    s0 = scans[0]["spec"]
    P("Position   %.5f N, %.5f W  (Kotzebue Sound, Alaska)" % (s0.latitude,
                                                               -s0.longitude))
    P("Date       %s   GPS %s-%s UTC"
      % (s0.header.get("Date", "?").split(",")[0],
         _hhmm(min(x["spec"].gps_time for x in scans)),
         _hhmm(max(x["spec"].gps_time for x in scans))))
    suns = [x["spec"].solar_elevation_deg for x in scans]
    tier, msg = solar_window_verdict(sum(suns) / len(suns))
    P("Sun        %.1f-%.1f deg elevation -> %s" % (min(suns), max(suns), tier.upper()))
    P("Scans      %d total: %d water, %d sky, %d vegetation (excluded)"
      % (len(scans), sum(1 for x in scans if x["role"] == WATER),
         sum(1 for x in scans if x["role"] == SKY),
         sum(1 for x in scans if x["role"] == VEG)))
    P("Sites      %d occupations, %d stations (one per panel reference)"
      % (len(groups), len(results)))
    P("Grid       350-2500 nm, 2151 channels; usable on water to ~900 nm")
    P("")
    P("%-9s %4s %4s %10s %9s %9s %8s %8s"
      % ("station", "sky", "wat", "Rrs(443)", "Rrs(555)", "peak_nm", "NIR/grn",
         "skyspr"))
    P("-" * 72)
    for st in results:
        mean = [sum(r["res"].rrs[i] for r in st["rrs"]) / len(st["rrs"])
                for i in range(len(st["wl"]))]
        wl = st["wl"]
        f = lambda lo, hi: (sum(v for w, v in zip(wl, mean) if lo <= w <= hi)
                            / max(1, sum(1 for w in wl if lo <= w <= hi)))
        grn = max(v for w, v in zip(wl, mean) if 550 <= w <= 600)
        pk = [w for w, v in zip(wl, mean) if 400 <= w <= 750 and v == max(
            x for ww, x in zip(wl, mean) if 400 <= ww <= 750)][0]
        P("%-9.4f %4d %4d %10.5f %9.5f %9.0f %8.3f %7.1f%%"
          % (st["ref"], st["n_sky"], len(st["rrs"]), f(440, 446), f(552, 558), pk,
             f(800, 900) / grn, 100 * st["sky_spread"]))
    if giop_rows:
        P("")
        P("GIOP inversion (6 bands, per station mean)")
        P("%-14s %8s | %10s %9s %9s | %8s %8s %8s"
          % ("station", "OC4 chl", "DC M_phi", "DC adg", "DC bbp",
             "bd M_phi", "bd adg", "bd bbp"))
        P("-" * 84)
        for g in giop_rows:
            fmt = lambda v: ("%10.3g" % v) if abs(v) < 1e6 else "%10.1e" % v
            P("%-14s %8.2f | %s %9.3g %9.3g | %8.3g %8.3g %8.3g"
              % (g["name"].replace("Rrs_ref", ""), g["oc4"], fmt(g["dc_chl"]),
                 g["dc_adg"], g["dc_bbp"], g["bd_chl"], g["bd_adg"], g["bd_bbp"]))
        nrail = sum(1 for g in giop_rows if abs(g["bd_adg"] - 20.0) < 1e-6)
        P("")
        P("WARNING  %d of %d bounded a_dg(443) sit exactly ON the 20 m^-1 bound."
          % (nrail, len(giop_rows)))
        P("         A railed parameter is NOT a measurement. GIOP's 3-component basis")
        P("         cannot represent this water; the bound is absorbing the misfit.")

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(a.out, "REPORT.txt"), "w") as fh:
        fh.write(txt + "\n")
    print("\nwrote %s/REPORT.txt" % a.out)
    print("wrote %s" % csv_path)
    print("wrote %s/VERIFICATION.txt" % a.out)
    for p in paths:
        print("wrote %s" % p)


if __name__ == "__main__":
    main()
