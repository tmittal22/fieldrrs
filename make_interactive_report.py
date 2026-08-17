"""Interactive (zoomable) HTML report for one location/foreoptic dataset.

    python make_interactive_report.py <by_location/LOC*/FOREOPTIC_FOVxx>

Writes analysis/REPORT.html: one self-contained file, no internet needed, with plotly.js
embedded. Every scan is its own named trace, so you can click the legend to isolate a
scan, hover to read exact values, and box-zoom into any feature.

The static PNGs from analyse_location.py stay: they are what goes in a paper. This is
for looking.
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from analyse_location import (RHO_MOBLEY1999, assert_same_dataset, hhmm,
                              match_by_angle, rrs_three_scan, stats)
from fieldrrs.rrs import rho_at_angle, view_zenith_from_tilt
from organize_by_location import fov_deg, survey
from process_field_day import band, land_reflectance

WLO, WHI, RLO, RHI = 350.0, 950.0, 400.0, 900.0
C = {"sky": "#7fb3d5", "water": "#1f7a99", "panel": "#d9534f", "land": "#2e7d32",
     "mean": "#c0392b"}
PANEL_R = 0.99


def sub(wl, v, lo=WLO, hi=WHI):
    m = [(w, x) for w, x in zip(wl, v) if lo <= w <= hi and x == x]
    return [a for a, _ in m], [b for _, b in m]


def fig_html(fig, first):
    return fig.to_html(full_html=False,
                       include_plotlyjs=("inline" if first else False),
                       config={"displaylogo": False,
                               "toImageButtonOptions": {"scale": 3}})


def f_pooled(scans, wl):
    groups = [("CALIBRATION PANEL", "rad_ref", "panel",
               sorted({round(band(s["spec"], "rad_ref", 450, 650), 6): s
                       for s in scans}.items())),
              ("SKY", "rad_target", "sky",
               [(i, s) for i, s in enumerate(x for x in scans if x["role"] == "sky")]),
              ("WATER", "rad_target", "water",
               [(i, s) for i, s in enumerate(x for x in scans if x["role"] == "water")]),
              ("LAND", "rad_target", "land",
               [(i, s) for i, s in enumerate(x for x in scans if x["role"] == "land")])]
    groups = [g for g in groups if g[3]]
    fig = make_subplots(rows=1, cols=len(groups),
                        subplot_titles=[g[0] for g in groups],
                        horizontal_spacing=0.045)
    for c, (title, key, role, items) in enumerate(groups, start=1):
        for _, s in items:
            x, y = sub(wl, s["spec"].columns[key])
            fig.add_trace(go.Scatter(
                x=x, y=y, name="%s %s" % (title[:3], s["n"]), legendgroup=title,
                line=dict(width=1.1, color=C[role]), opacity=0.65,
                hovertemplate="%s<br>%%{x:.1f} nm<br>%%{y:.4e}<extra></extra>" % s["n"]),
                row=1, col=c)
        cur = [s["spec"].columns[key] for _, s in items]
        m, _, _ = stats(cur)
        x, y = sub(wl, m)
        fig.add_trace(go.Scatter(x=x, y=y, name="%s mean" % title[:3],
                                 legendgroup=title,
                                 line=dict(width=3, color="black")), row=1, col=c)
        fig.update_yaxes(type="log", row=1, col=c)
        fig.update_xaxes(title_text="nm", row=1, col=c)
    fig.update_yaxes(title_text="W m⁻² sr⁻¹ nm⁻¹", row=1, col=1)
    fig.update_layout(height=470, title="1 · The measured radiances, every scan",
                      hovermode="closest", template="plotly_white")
    return fig


def f_ed(scans, panels, wl, sun):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "giop_python", "src"))
    from giop.water import f0_solar
    wl_a = np.array(wl)
    eds = [np.array([math.pi * p / PANEL_R for p in s["spec"].columns["rad_ref"]])
           for s in panels]
    ed = sum(eds) / len(eds)
    mu = math.cos(math.radians(90.0 - sun))
    T = ed / (f0_solar(wl_a) * mu)
    m = (wl_a >= 380) & (wl_a <= 950)
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.09,
                        subplot_titles=("E_d vs top of atmosphere",
                                        "total transmittance T"))
    fig.add_trace(go.Scatter(x=wl_a[m], y=ed[m], name="E_d from panel",
                             line=dict(color=C["panel"], width=2.5)), row=1, col=1)
    fig.add_trace(go.Scatter(x=wl_a[m], y=(f0_solar(wl_a) * mu)[m], name="F₀·cosθs",
                             line=dict(color="#8a6000", width=2, dash="dash")),
                  row=1, col=1)
    fig.add_trace(go.Scatter(x=wl_a[m], y=T[m], name="T", line=dict(color="#2e7d32",
                             width=2.5)), row=1, col=2)
    fig.add_hline(y=1.0, line=dict(color="#c0392b", width=1.5), row=1, col=2)
    for lam, lab in ((762, "O₂-A"), (940, "H₂O")):
        fig.add_vline(x=lam, line=dict(color="#999", dash="dot"), row=1, col=2,
                      annotation_text=lab)
    fig.update_yaxes(title_text="W m⁻² nm⁻¹", row=1, col=1)
    fig.update_yaxes(title_text="T", range=[0, 1.15], row=1, col=2)
    fig.update_xaxes(title_text="nm")
    fig.update_layout(height=430, template="plotly_white",
                      title="3 · Is E_d reasonable?  median T(450-650) = %.3f, "
                            "sun %.1f°" % (float(np.median(T[(wl_a >= 450) &
                                                             (wl_a <= 650)])), sun))
    return fig


def f_rrs(water, sky, wl):
    def rrs_of(w, sk, angle=True):
        r = rho_at_angle(view_zenith_from_tilt(w["spec"].tilt_y_deg)) if angle \
            else RHO_MOBLEY1999
        return rrs_three_scan(wl, w["spec"].columns["rad_target"],
                              sk["spec"].columns["rad_target"],
                              w["spec"].columns["rad_ref"], PANEL_R, r, "none").rrs

    pairs = match_by_angle(water, sky)
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.09,
                        subplot_titles=("every pairing (blind, fixed ρ)",
                                        "angle-matched, ρ(θv), one sky per water"))
    for w in water:
        for sk in sky:
            x, y = sub(wl, rrs_of(w, sk, False), RLO, RHI)
            fig.add_trace(go.Scatter(x=x, y=y, name="w%s·s%s" % (w["n"], sk["n"]),
                                     line=dict(width=0.7, color="#bbb"),
                                     showlegend=False, hoverinfo="skip"),
                          row=1, col=1)
    A = [rrs_of(w, sk, False) for w in water for sk in sky]
    mA, _, _ = stats(A)
    x, y = sub(wl, mA, RLO, RHI)
    fig.add_trace(go.Scatter(x=x, y=y, name="blind mean",
                             line=dict(width=3, color=C["mean"])), row=1, col=1)
    Cc = []
    for w, sk, dm, _ in pairs:
        v = rrs_of(w, sk, True)
        Cc.append(v)
        x, y = sub(wl, v, RLO, RHI)
        fig.add_trace(go.Scatter(
            x=x, y=y, name="w%s + s%s" % (w["n"], sk["n"]),
            line=dict(width=1.2, color=C["water"]), opacity=0.8,
            hovertemplate="water %s θv=%.1f°<br>sky %s θv=%.1f°<br>Δ=%.1f°<br>"
                          "%%{x:.1f} nm  %%{y:.5f}<extra></extra>"
                          % (w["n"], w["spec"].tilt_y_deg, sk["n"],
                             sk["spec"].tilt_y_deg, dm)), row=1, col=2)
    mC, _, _ = stats(Cc)
    x, y = sub(wl, mC, RLO, RHI)
    fig.add_trace(go.Scatter(x=x, y=y, name="angle-matched mean",
                             line=dict(width=3.5, color="#2e7d32")), row=1, col=2)
    fig.update_xaxes(title_text="nm", range=[RLO, RHI])
    fig.update_yaxes(title_text="R_rs  sr⁻¹")
    i443 = min(range(len(wl)), key=lambda k: abs(wl[k] - 443))
    sd = np.std([c[i443] for c in Cc]) / np.mean([c[i443] for c in Cc]) * 100
    fig.update_layout(height=470, template="plotly_white",
                      title="4 · R_rs — R_rs(443) = %.5f sr⁻¹, sd %.1f %% (n=%d)"
                            % (mC[i443], sd, len(Cc)))
    return fig, mC


def f_geometry(water, sky, wl):
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.10,
                        subplot_titles=("ρ used, per scan",
                                        "R_rs(443) vs view angle"))
    tt = np.linspace(35, 55, 60)
    fig.add_trace(go.Scatter(x=tt, y=[rho_at_angle(t) for t in tt],
                             name="ρ(θv) Fresnel-scaled",
                             line=dict(color="#2e7d32", width=2.5)), row=1, col=1)
    fig.add_hline(y=RHO_MOBLEY1999, line=dict(color="#c0392b", dash="dash"),
                  row=1, col=1, annotation_text="fixed 0.028")
    tw = [w["spec"].tilt_y_deg for w in water]
    fig.add_trace(go.Scatter(x=tw, y=[rho_at_angle(view_zenith_from_tilt(t))
                                      for t in tw], mode="markers", name="water scans",
                             marker=dict(size=11, color=C["water"],
                                         line=dict(width=1, color="black")),
                             text=[w["n"] for w in water],
                             hovertemplate="%{text}<br>θv=%{x:.1f}°<br>ρ=%{y:.5f}"
                                           "<extra></extra>"), row=1, col=1)
    l_sky = [sum(s["spec"].columns["rad_target"][i] for s in sky) / len(sky)
             for i in range(len(wl))]
    i443 = min(range(len(wl)), key=lambda k: abs(wl[k] - 443))
    for lab, angle, col in (("fixed ρ", False, "#c0392b"),
                            ("ρ(θv)", True, "#2e7d32")):
        y = []
        for w in water:
            r = rho_at_angle(view_zenith_from_tilt(w["spec"].tilt_y_deg)) if angle \
                else RHO_MOBLEY1999
            y.append(rrs_three_scan(wl, w["spec"].columns["rad_target"], l_sky,
                                    w["spec"].columns["rad_ref"], PANEL_R, r,
                                    "none").rrs[i443])
        from scipy import stats as sps
        r_, p_ = sps.pearsonr(tw, y)
        fig.add_trace(go.Scatter(x=tw, y=y, mode="markers",
                                 name="%s  r=%+.2f p=%.3f" % (lab, r_, p_),
                                 marker=dict(size=11, color=col,
                                             line=dict(width=1, color="black")),
                                 text=[w["n"] for w in water],
                                 hovertemplate="%{text}<br>θv=%{x:.1f}°<br>"
                                               "R_rs=%{y:.5f}<extra></extra>"),
                      row=1, col=2)
    fig.update_xaxes(title_text="view zenith angle (deg)")
    fig.update_yaxes(title_text="ρ", row=1, col=1)
    fig.update_yaxes(title_text="R_rs(443)  sr⁻¹", row=1, col=2)
    fig.update_layout(height=440, template="plotly_white",
                      title="5 · Geometry: per-scan ρ and the residual angle trend")
    return fig


def f_land(land, wl):
    fig = go.Figure()
    for s in land:
        w, r = land_reflectance(s["spec"], PANEL_R)
        x, y = sub(w, r, 380, 950)
        fig.add_trace(go.Scatter(x=x, y=y, name="%s  red edge %.1f×"
                                 % (s["n"], s["diag"]["red_edge"]),
                                 line=dict(width=2.5),
                                 hovertemplate="%{x:.1f} nm  R=%{y:.4f}<extra></extra>"))
    fig.add_vrect(x0=665, x1=680, fillcolor="#c0392b", opacity=0.12, line_width=0,
                  annotation_text="chlorophyll")
    fig.add_vrect(x0=700, x1=760, fillcolor="#2e7d32", opacity=0.10, line_width=0,
                  annotation_text="red edge")
    fig.update_layout(height=440, template="plotly_white",
                      xaxis_title="nm", yaxis_title="reflectance factor R",
                      title="6 · Land targets — reflectance, NOT R_rs "
                            "(no ρ, no sky scan)")
    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    a = ap.parse_args()
    scans = survey(a.folder)
    span = assert_same_dataset(scans)
    wl = scans[0]["spec"].wavelength
    sky = sorted([s for s in scans if s["role"] == "sky"], key=lambda x: x["n"])
    water = sorted([s for s in scans if s["role"] == "water"], key=lambda x: x["n"])
    land = sorted([s for s in scans if s["role"] == "land"], key=lambda x: x["n"])
    seen, panels = set(), []
    for s in scans:
        k = round(band(s["spec"], "rad_ref", 450, 650), 6)
        if k not in seen:
            seen.add(k); panels.append(s)
    sun = sum(s["sun"] for s in scans) / len(scans)
    loc = os.path.basename(os.path.dirname(a.folder.rstrip("/")))
    fo = os.path.basename(a.folder.rstrip("/"))
    outdir = os.path.join(a.folder, "analysis")
    os.makedirs(outdir, exist_ok=True)

    figs = [f_pooled(scans, wl), f_ed(scans, panels, wl, sun)]
    fr, mC = f_rrs(water, sky, wl)
    figs.append(fr)
    figs.append(f_geometry(water, sky, wl))
    if land:
        figs.append(f_land(land, wl))

    report = os.path.join(outdir, "REPORT.txt")
    notes = open(report).read() if os.path.exists(report) else ""

    parts = ["""<!doctype html><html><head><meta charset="utf-8">
<title>%s · %s</title><style>
body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:0 auto;max-width:1500px;
padding:26px;color:#1a1a1a;background:#fff}
h1{font-size:23px;margin-bottom:4px}h2{font-size:17px;margin-top:34px;color:#33506e}
.meta{color:#555;font-size:14px;margin-bottom:18px}
pre{background:#f7f7f7;border-left:4px solid #8aa8c8;padding:14px 16px;font-size:12.5px;
overflow-x:auto;line-height:1.45}
.tip{background:#eef5ff;border-left:4px solid #2c6f9b;padding:10px 14px;font-size:13.5px;
margin:14px 0}
</style></head><body>
<h1>%s &middot; %s</h1>
<div class="meta">%d scans &middot; %d water, %d sky, %d land &middot; %s&ndash;%s UTC
&middot; sun %.1f&ndash;%.1f&deg; &middot; FOV %.0f&deg; &middot; all within %.0f m</div>
<div class="tip"><b>Interactive.</b> Box-zoom to inspect a feature; double-click to reset.
Click a legend entry to hide one trace, double-click it to isolate it. Hover reads exact
values, and on the R<sub>rs</sub> panel it shows which sky was paired with which water and
the angle mismatch. Camera icon saves a PNG.</div>
""" % (loc, fo, loc, fo, len(scans), len(water), len(sky), len(land),
       hhmm(min(s["gps"] for s in scans)), hhmm(max(s["gps"] for s in scans)),
       min(s["sun"] for s in scans), max(s["sun"] for s in scans),
       fov_deg(fo.split("_")[0]), span)]

    for i, f in enumerate(figs):
        parts.append(fig_html(f, i == 0))
    if notes:
        parts.append("<h2>Full numeric report</h2><pre>%s</pre>" % notes)
    parts.append("</body></html>")

    out = os.path.join(outdir, "REPORT.html")
    with open(out, "w") as fh:
        fh.write("\n".join(parts))
    print("wrote %s  (%.1f MB, self-contained)"
          % (out, os.path.getsize(out) / 1e6))


if __name__ == "__main__":
    main()
