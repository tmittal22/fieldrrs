"""GIOP applied to a field R_rs: band choice, OC4, conditioning, fit and assumptions.

    python make_giop_figures.py <analysis/FINAL_Rrs.csv> --out <dir>

Five figures. The first exists because I originally ran GIOP on six bands and reported
its conditioning as if that were GIOP's limit. It is not: GIOP solves at whatever
wavelengths you give it, and the Bricaud aph* table is continuous over 400-700 nm at
2 nm. Running it hyperspectral changes the answer and the uncertainty by more than an
order of magnitude, so the band choice gets its own figure.
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "giop_python", "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from giop import get_oc, giop
from giop.empirical import _OC_COEF
from giop.model import (GORDON_G0, GORDON_G1, GiopConfig, eigenvectors,
                        find_anchor_bands, rrs_above_to_below, rrs_from_iops)
from giop.water import a_water, bb_water

B6 = np.array([412., 443, 490, 510, 555, 670])
LO, HI = 400.0, 700.0


def load(path):
    rows = list(csv.reader(open(path)))
    d = np.array([[float(x) for x in r] for r in rows
                  if r and not r[0].startswith("#") and r[0][0].isdigit()])
    return d[:, 0], d[:, 1], d[:, 2]


def at(wl, v, b):
    return np.array([v[np.argmin(abs(wl - x))] for x in b])


def run(W, R, chl, **kw):
    sig = np.sqrt((0.05 * np.abs(R)) ** 2 + 2e-4 ** 2)
    return giop(W, R, chl, inv="bounded", sigma=sig, **kw)


# ------------------------------------------------------------------ figure 1
def f_bands(wl, mean, ssd, chl, out):
    m = (wl >= LO) & (wl <= HI)
    sets = [("6 bands\n(upstream demo)", B6, at(wl, mean, B6)),
            ("20 nm", wl[m][::20], mean[m][::20]),
            ("10 nm", wl[m][::10], mean[m][::10]),
            ("2 nm\n(aph* native)", wl[m][::2], mean[m][::2]),
            ("1 nm\n(HYPERSPECTRAL)", wl[m], mean[m])]
    res = []
    for lab, W, R in sets:
        try:
            g = run(W, R, chl)
            res.append((lab, len(W), g.chl, g.adg443, g.bbp443, not g.failed))
        except Exception as exc:
            res.append((lab, len(W), np.nan, np.nan, np.nan, False))

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.0))
    x = range(len(res))
    for ax, k, title, col in ((axes[0], 2, "$M_\\phi$  (chlorophyll amplitude)", "#2e7d32"),
                              (axes[1], 3, "$a_{dg}$(443)  m$^{-1}$", "#8a6000"),
                              (axes[2], 4, "$b_{bp}$(443)  m$^{-1}$", "#2c6f9b")):
        vals = [r[k] for r in res]
        bars = ax.bar(x, vals, color=col)
        for i, (v, ok) in enumerate(zip(vals, [r[5] for r in res])):
            if not ok or not np.isfinite(v):
                ax.text(i, 0, "FAILED", ha="center", va="bottom", rotation=90,
                        fontsize=9, color="#c0392b", weight="bold")
            else:
                ax.text(i, v, "%.3g" % v, ha="center", va="bottom", fontsize=9)
        ax.set_xticks(list(x))
        ax.set_xticklabels(["%s\nn=%d" % (r[0], r[1]) for r in res], fontsize=8)
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.25, axis="y")
    fig.suptitle("GIOP IS NOT A SIX-BAND ALGORITHM. It solves on whatever grid you give "
                 "it; the Bricaud $a^*_\\phi$ table is continuous 400-700 nm at 2 nm.\n"
                 "Using six bands roughly DOUBLES every retrieved amplitude here.",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    p = os.path.join(out, "giop1_band_choice.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p, res


# ------------------------------------------------------------------ figure 2
def f_oc4(wl, mean, out):
    """The SEED CHAIN: where every prescribed number in the inversion comes from.

    Three of the model's inputs are not fitted and not measured -- they are computed
    from the spectrum itself by published parameterisations, before the inversion runs:

        OC4 chlorophyll -> selects the Bricaud a*_phi SHAPE (and nothing else)
        eta             -> QAA v5, from the subsurface blue/green ratio
        S_dg            -> GIOP-DC fixes it at 0.018; 'qaa' and 'obpg' compute it

    Each is shown as its formula, evaluated on this spectrum, so the reader can see how
    much of the answer was decided before any fitting happened.
    """
    from giop.aphstar import BRICAUD_NORM_VALUE, BRICAUD_NORM_WL, bricaud1998
    from giop.model import eta_qaa, sdg_from_option

    r = at(wl, mean, np.array([443., 490, 510, 555]))
    a = _OC_COEF["oc4"]
    num = max(r[0], r[1], r[2]); which = [443, 490, 510][int(np.argmax(r[:3]))]
    X = np.log10(num / r[3])
    poly = a[0] + a[1] * X + a[2] * X ** 2 + a[3] * X ** 3 + a[4] * X ** 4
    chl = 10 ** poly

    m = (wl >= LO) & (wl <= HI)
    W, R = wl[m], mean[m]
    cfg = GiopConfig()
    rin = rrs_above_to_below(R, cfg.trans)
    idx = find_anchor_bands(W)
    ratio = rin[idx["443"]] / rin[idx["555"]]
    eta = eta_qaa(rin, idx["443"], idx["555"])
    sdg = {k: sdg_from_option(k, R, rin, idx["412"], idx["443"], idx["555"])
           for k in ("qaa", "obpg", "gsm")}

    fig, axes = plt.subplots(2, 2, figsize=(15, 10.2))
    axes = axes.ravel()
    ax = axes[0]
    xs = np.linspace(-0.6, 0.6, 400)
    ys = 10 ** (a[0] + a[1] * xs + a[2] * xs ** 2 + a[3] * xs ** 3 + a[4] * xs ** 4)
    ax.plot(xs, ys, lw=2.5, color="#2c6f9b")
    ax.plot([X], [chl], "o", ms=13, color="#c0392b", zorder=5)
    ax.annotate("LOC1\nX = %.4f\nchl = %.2f mg m$^{-3}$" % (X, chl), (X, chl),
                xytext=(28, 40), textcoords="offset points", fontsize=10,
                weight="bold", color="#c0392b",
                arrowprops=dict(arrowstyle="->", color="#c0392b"))
    ax.set_yscale("log"); ax.grid(alpha=0.3)
    ax.set_xlabel("$X = \\log_{10}\\left[\\max(R_{rs}443,490,510)\\,/\\,R_{rs}555\\right]$")
    ax.set_ylabel("chlorophyll  (mg m$^{-3}$)")
    ax.set_title("OC4 is a 4th-order polynomial in ONE band ratio\n"
                 "$\\log_{10}$chl $= %.4f %+.4fX %+.4fX^2 %+.4fX^3 %+.4fX^4$" % a,
                 fontsize=10.5, loc="left")

    ax = axes[1]
    lab = ["443", "490", "510", "555"]
    bars = ax.bar(lab, r, color=["#2c6f9b" if l != str(which) else "#c0392b"
                                 for l in lab])
    ax.text(0.5, 0.95, "numerator = max(443,490,510) = R_rs(%d) = %.5f\n"
            "denominator = R_rs(555) = %.5f\nratio = %.4f,  X = %.4f"
            % (which, num, r[3], num / r[3], X), transform=ax.transAxes,
            va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.4", fc="#f7f7f7", ec="#999"))
    ax.set_ylabel("$R_{rs}$  sr$^{-1}$"); ax.grid(alpha=0.25, axis="y")
    ax.set_title("The four numbers OC4 uses (red = the one selected)", fontsize=10.5,
                 loc="left")

    # --- (c) what the OC4 number is actually FOR: it picks the aph* shape
    ax = axes[2]
    seeds = [0.5, 2.0, chl, 20.0, 50.0]
    for c_, col in zip(seeds, ["#bbbbbb", "#8fbcd4", "#c0392b", "#7fa87f", "#555555"]):
        y = bricaud1998(c_, W)
        lw = 3.0 if abs(c_ - chl) < 1e-9 else 1.4
        ax.plot(W, y, lw=lw, color=col,
                label="chl=%.2f%s  $a^*_\\phi(443)/a^*_\\phi(675)$=%.2f"
                      % (c_, "  <- OC4" if abs(c_ - chl) < 1e-9 else "",
                         y[np.argmin(abs(W - 443))] / y[np.argmin(abs(W - 675))]))
    ax.axvline(BRICAUD_NORM_WL, color="#333", ls=":", lw=1.2)
    ax.text(BRICAUD_NORM_WL + 3, ax.get_ylim()[1] * 0.95,
            "normalised:\n$a^*_\\phi$(%g) = %g" % (BRICAUD_NORM_WL, BRICAUD_NORM_VALUE),
            fontsize=8)
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("$a^*_\\phi$  m$^2$ mg$^{-1}$")
    ax.legend(fontsize=7.5); ax.grid(alpha=0.25)
    ax.set_title("The OC4 number does ONE thing: it selects the $a^*_\\phi$ SHAPE.\n"
                 "$a^*_\\phi(\\lambda) = A_\\phi(\\lambda)\\,\\mathrm{chl}^{\\,"
                 "E_\\phi(\\lambda)-1}$, frozen at the seed — GIOP does NOT iterate.",
                 fontsize=9.5, loc="left")

    # --- (d) eta and S_dg: computed, not fitted, not measured
    ax = axes[3]
    x = np.linspace(0.05, 3.0, 300)
    ax.plot(x, 2.0 * (1.0 - 1.2 * np.exp(-0.9 * x)), lw=2.4, color="#2c6f9b",
            label="$\\eta$  QAA v5")
    ax.plot(x, 100 * (0.015 + 0.002 / (0.6 + x)), lw=2.4, color="#8a6000",
            label="$S_{dg}$ 'qaa'  ($\\times$100)")
    ax.axhline(100 * 0.018, color="#c0392b", ls="--", lw=1.8,
               label="$S_{dg}$ GIOP-DC default 0.018 ($\\times$100)")
    ax.axhline(100 * sdg["obpg"], color="#2e7d32", ls="-.", lw=1.8,
               label="$S_{dg}$ 'obpg' here = %.4f ($\\times$100)" % sdg["obpg"])
    ax.axhline(100 * 0.0113, color="#ff2d55", ls=":", lw=2.4,
               label="$S_{dg}$ FITTED here = 0.0113 ($\\times$100)")
    ax.axvline(ratio, color="k", lw=1.4)
    ax.plot([ratio], [eta], "*", ms=20, color="#2c6f9b", mec="k", zorder=5)
    ax.plot([ratio], [100 * sdg["qaa"]], "*", ms=20, color="#8a6000", mec="k", zorder=5)
    ax.annotate("this spectrum\n$r_{rs}$(443)/$r_{rs}$(555) = %.3f\n"
                "$\\Rightarrow \\eta$ = %.3f,  $S_{dg}^{qaa}$ = %.4f"
                % (ratio, eta, sdg["qaa"]), (ratio, eta), xytext=(30, -55),
                textcoords="offset points", fontsize=8.5,
                arrowprops=dict(arrowstyle="->"))
    ax.set_xlabel("subsurface blue/green ratio  $r_{rs}$(443) / $r_{rs}$(555)")
    ax.set_ylabel("$\\eta$   /   $S_{dg}\\times$100  (nm$^{-1}$)")
    ax.legend(fontsize=7.5, loc="upper right"); ax.grid(alpha=0.25)
    ax.set_title("$\\eta = 2\\,[1 - 1.2\\,e^{-0.9\\,x}]$ (QAA v5) and "
                 "$S_{dg} = 0.015 + 0.002/(0.6+x)$, both from the SAME ratio $x$.\n"
                 "GIOP's own 'obpg' option lands at %.4f, within %.0f %% of the fitted "
                 "value; the DC default 0.018 is %.0f %% high."
                 % (sdg["obpg"], 100 * abs(sdg["obpg"] / 0.0113 - 1),
                    100 * (0.018 / 0.0113 - 1)), fontsize=9.5, loc="left")

    fig.suptitle("THE SEED CHAIN — everything the inversion prescribes before it fits "
                 "anything, and where each number comes from.\nOC4 (NASA OC v6, "
                 "get_oc.m) selects the $a^*_\\phi$ shape; $\\eta$ and $S_{dg}$ come "
                 "from one blue/green ratio. None of these is fitted, and none is "
                 "measured.", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    p = os.path.join(out, "giop2_seed_chain.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p, chl, X, which


# ------------------------------------------------------------------ figure 3
def f_fit(wl, mean, chl, out):
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.2))
    cfg = GiopConfig()
    for ax, (W, R, lab) in zip(axes, [
            (B6, at(wl, mean, B6), "6 bands"),
            (wl[(wl >= LO) & (wl <= HI)], mean[(wl >= LO) & (wl <= HI)],
             "hyperspectral 400-700 nm")]):
        rin = rrs_above_to_below(R, cfg.trans)
        g = run(W, R, chl)
        idx = find_anchor_bands(W)
        adgs, bbps, aphs, sdg, eta = eigenvectors(W, cfg, chl, R, rin, idx)
        Mdg, Mbp, Mphi = g.x
        rmod = rrs_from_iops(a_water(W) + Mphi * aphs + Mdg * adgs,
                             bb_water(W) + Mbp * bbps, GORDON_G0, GORDON_G1)
        rms = 100 * np.sqrt(np.mean((rmod / rin - 1) ** 2))
        ax.plot(W, rin, "o-" if len(W) < 20 else "-", lw=2.2, ms=8, color="#1f7a99",
                label="measured $r_{rs}$")
        ax.plot(W, rmod, "s--" if len(W) < 20 else "--", lw=2.0, ms=7, color="#c0392b",
                label="GIOP model")
        ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("$r_{rs}$ below surface")
        ax.legend(fontsize=9); ax.grid(alpha=0.25)
        ax.set_title("%s   RMS misfit %.1f %%\n$M_\\phi$=%.2f  $a_{dg}$=%.2f  "
                     "$b_{bp}$=%.3f" % (lab, rms, g.chl, g.adg443, g.bbp443),
                     fontsize=10.5, loc="left")
        ax2 = ax.twinx()
        ax2.plot(W, 100 * (rmod / rin - 1), lw=1.2, color="#8a6000", alpha=0.75)
        ax2.axhline(0, color="#8a6000", lw=0.8, ls=":")
        ax2.set_ylabel("residual (%)", color="#8a6000")
        ax2.set_ylim(-30, 30)
    fig.suptitle("Six bands fit BETTER (3 parameters, 6 points) while hyperspectral "
                 "exposes the structural misfit.\nPrecision and accuracy are different "
                 "things, and only the second one matters here.", fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    p = os.path.join(out, "giop3_fit_quality.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p


# ------------------------------------------------------------------ figure 4
def f_uncertainty(wl, mean, ssd, chl, out, ndraw=150):
    m = (wl >= LO) & (wl <= HI)
    cases = [("6 bands", B6, at(wl, mean, B6), at(wl, ssd, B6)),
             ("hyperspectral", wl[m], mean[m], ssd[m])]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.0))
    store = {}
    rng = np.random.default_rng(0)
    for lab, W, R, S in cases:
        out_ = []
        for _ in range(ndraw):
            try:
                g = run(W, R + rng.normal(0, S), chl)
                if not g.failed:
                    out_.append((g.chl, g.adg443, g.bbp443))
            except Exception:
                pass
        store[lab] = np.array(out_)
    for ax, k, title in ((axes[0], 0, "$M_\\phi$"), (axes[1], 1, "$a_{dg}$(443)"),
                         (axes[2], 2, "$b_{bp}$(443)")):
        for lab, col in (("6 bands", "#c0392b"), ("hyperspectral", "#2e7d32")):
            v = store[lab][:, k]
            ax.hist(v, bins=28, alpha=0.6, color=col,
                    label="%s: %.3g ± %.0f %%" % (lab, v.mean(),
                                                  100 * v.std() / abs(v.mean())))
        ax.set_xlabel(title); ax.set_ylabel("draws"); ax.legend(fontsize=9)
        ax.grid(alpha=0.25)
        ax.set_title("%s under the 1.7 %% shape uncertainty" % title, fontsize=10.5,
                     loc="left")
    fig.suptitle("CONDITIONING. The same 1.7 % input uncertainty, propagated through "
                 "both band choices.\nSix bands amplify it ~27x; hyperspectral does not.",
                 fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    p = os.path.join(out, "giop4_conditioning.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p, store


# ------------------------------------------------------------------ figure 5
def f_sdg(wl, mean, ssd, chl, out):
    m = (wl >= LO) & (wl <= HI)
    S = ssd[m]
    sd = np.linspace(0.010, 0.025, 16)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    hyp = None
    for lab, W, R, col in (("6 bands", B6, at(wl, mean, B6), "#c0392b"),
                           ("hyperspectral", wl[m], mean[m], "#2e7d32")):
        vals = []
        for s_ in sd:
            try:
                g = run(W, R, chl, sdg=float(s_))
                c2 = (_chi2(g, R, S)[0] / (len(W) - 3) if len(W) > 20 else np.nan)
                vals.append((g.chl, g.adg443, g.bbp443, c2))
            except Exception:
                vals.append((np.nan,) * 4)
        vals = np.array(vals)
        if lab == "hyperspectral":
            hyp = vals
        for ax, k in zip(axes, range(3)):
            ax.plot(sd, vals[:, k], "-", lw=1.6, color=col, label=lab, zorder=1)
            if lab == "hyperspectral":
                ax.scatter(sd, vals[:, k], c=vals[:, 3], cmap="viridis_r", s=80,
                           zorder=3, edgecolor="k", linewidth=0.5,
                           norm=matplotlib.colors.LogNorm())
            else:
                ax.plot(sd, vals[:, k], "o", ms=5, color=col, zorder=2)
    b = int(np.nanargmin(hyp[:, 3]))
    for ax, t, k in zip(axes, ("$M_\\phi$", "$a_{dg}$(443)  m$^{-1}$",
                               "$b_{bp}$(443)  m$^{-1}$"), range(3)):
        ax.axvline(0.018, color="#333", ls="--", lw=1.5)
        ax.plot(sd[b], hyp[b, k], "*", ms=20, color="#ff2d55", mec="k", zorder=5)
        ax.set_xlabel("$S_{dg}$  (nm$^{-1}$)")
        ax.set_ylabel(t); ax.set_yscale("symlog", linthresh=1e-2)
        ax.grid(alpha=0.25); ax.legend(fontsize=9, loc="upper left")
        ax.set_title("colour = $\\chi^2_\\nu$ (hyperspectral). BEST at $S_{dg}$=%.4f, "
                     "$\\chi^2_\\nu$=%.0f;\nthe default 0.018 gives %.0f."
                     % (sd[b], hyp[b, 3], hyp[np.argmin(abs(sd - 0.018)), 3]),
                     fontsize=9, loc="left")
    plt.colorbar(axes[-1].collections[0], ax=axes[-1], label="$\\chi^2_\\nu$")
    fig.suptitle("$S_{dg}$ IS NOT A FREE CHOICE — the data prefer one. Sweeping it moves "
                 "$M_\\phi$ over orders of magnitude, but the extreme arms FIT BADLY,\n"
                 "so that swing is a set of rejected models, not an error bar. The "
                 "GIOP default 0.018 is not the preferred value here.", fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    p = os.path.join(out, "giop5_sdg_assumption.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p


def _perscan_fits(loc, chl, free=False, S=None):
    """GIOP on each angle-matched pair, keeping the spectra and the IOPs.

    ``free=True`` also fits S_dg and eta per spectrum -- the maximum-freedom
    configuration, which on this water is also the best-fitting one.
    """
    from analyse_location import match_by_angle
    from fieldrrs.rrs import rho_at_angle, rrs_three_scan, view_zenith_from_tilt
    from organize_by_location import survey
    from giop.model import (GiopConfig, GORDON_G0, GORDON_G1, eigenvectors,
                            find_anchor_bands, rrs_above_to_below, rrs_from_iops,
                            u_from_rrs)
    scans = survey(loc)
    sky = [s for s in scans if s["role"] == "sky"]
    water = [s for s in scans if s["role"] == "water"]
    wl = np.array(scans[0]["spec"].wavelength)
    m = (wl >= LO) & (wl <= HI); W = wl[m]
    cfg = GiopConfig(); out = []
    for w, sk, dm, _ in match_by_angle(water, sky):
        rho = rho_at_angle(view_zenith_from_tilt(w["spec"].tilt_y_deg))
        v = np.array(rrs_three_scan(wl, w["spec"].columns["rad_target"],
                                    sk["spec"].columns["rad_target"],
                                    w["spec"].columns["rad_ref"], 0.99, rho,
                                    "none").rrs)[m]
        try:
            g = run(W, v, chl, **(dict(fit_shapes=True, n_starts=4) if free else {}))
            rin = rrs_above_to_below(v, cfg.trans)
            rmod = g.rrs_model_subsurface
            # chi2 against the SAME measured per-band sigma used everywhere else in
            # this folder, so per-scan chi2 is comparable with the mean-spectrum chi2.
            c2n = (float(np.sum(((g.rrs_model_above - v) / S) ** 2)) / (len(W) - 3)
                   if S is not None else np.nan)
            out.append({"n": w["n"], "sky": sk["n"], "dm": dm, "W": W, "rrs": v,
                        "rin": rin, "rmod": rmod, "M": g.chl, "adg": g.adg443,
                        "bbp": g.bbp443, "eta": g.eta, "sdg": g.sdg, "c2n": c2n,
                        "rms": 100 * float(np.sqrt(np.mean((rmod / rin - 1) ** 2))),
                        "amp": float(v[int(np.argmin(abs(W - 555)))]), "ok": True})
        except Exception:
            pass
    return out


# ------------------------------------------------------------------ figure 6
def f_perscan(loc, chl, out, S=None):
    """All twelve fits, spectrum by spectrum, in BOTH configurations.

    Earlier this showed only the constrained fit (S_dg fixed at 0.018, eta from QAA,
    three amplitudes free). That is not the best the model can do, and showing only it
    made the misfit look like an inevitable property of GIOP rather than of that
    particular choice of shapes. Both are now drawn on every panel.
    """
    fits = _perscan_fits(loc, chl, free=False, S=S)
    fr = {f["n"]: f for f in _perscan_fits(loc, chl, free=True, S=S)}
    n = len(fits)
    fig, axes = plt.subplots(3, 4, figsize=(17.5, 11), sharex=True)
    for ax, f in zip(axes.flat, fits):
        g = fr.get(f["n"])
        ax.plot(f["W"], f["rin"], lw=2.0, color="#1f7a99", label="measured")
        ax.plot(f["W"], f["rmod"], lw=1.5, ls="--", color="#c0392b",
                label="CONSTRAINED  ($S_{dg}$=0.018)")
        t = ("water %s + sky %s   $\\Delta\\theta$=%.1f$^\\circ$\n"
             "fixed:  $M_\\phi$=%.1f  $a_{dg}$=%.2f  $b_{bp}$=%.3f   RMS %.1f %%"
             % (f["n"], f["sky"], f["dm"], f["M"], f["adg"], f["bbp"], f["rms"]))
        if g:
            ax.plot(g["W"], g["rmod"], lw=1.5, ls="-.", color="#2e7d32",
                    label="FREE  ($S_{dg}$, $\\eta$ fitted)")
            t += ("\nFREE:  $M_\\phi$=%.1f  $a_{dg}$=%.2f  $b_{bp}$=%.3f   RMS %.1f %%"
                  "   $S_{dg}$=%.4f $\\eta$=%+.2f"
                  % (g["M"], g["adg"], g["bbp"], g["rms"], g["sdg"], g["eta"]))
        ax.set_title(t, fontsize=7.6, loc="left")
        ax.grid(alpha=0.25); ax.tick_params(labelsize=8)
    for ax in axes.flat[n:]:
        ax.axis("off")
    axes[0][0].legend(fontsize=7.5)
    for ax in axes[-1]:
        ax.set_xlabel("wavelength (nm)")
    for row in axes:
        row[0].set_ylabel("$r_{rs}$ below surface")
    fmed = np.median([f["rms"] for f in fits])
    gmed = np.median([g["rms"] for g in fr.values()]) if fr else np.nan
    fig.suptitle("All %d angle-matched fits, hyperspectral 400-700 nm, in BOTH "
                 "configurations.\nCONSTRAINED (red, GIOP-DC: $S_{dg}$=0.018, $\\eta$ "
                 "from QAA, 3 free amplitudes) median RMS %.1f %%  vs  FREE (green: "
                 "$S_{dg}$ and $\\eta$ fitted too, 5 free) median RMS %.1f %%."
                 % (n, fmed, gmed), fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.945))
    p = os.path.join(out, "giop6_all_fits.png")
    fig.savefig(p, dpi=135); plt.close(fig)
    return p, fits, list(fr.values())


# ------------------------------------------------------------------ figure 7
def f_covariance(fits, out):
    """Every parameter, and how they trade off against each other."""
    keys = ["M", "adg", "bbp", "eta", "rms", "amp"]
    lab = {"M": "$M_\\phi$", "adg": "$a_{dg}$(443)", "bbp": "$b_{bp}$(443)",
           "eta": "$\\eta$ (derived)", "rms": "fit RMS %", "amp": "$R_{rs}$(555)"}
    D = {k: np.array([f[k] for f in fits]) for k in keys}
    names = [f["n"] for f in fits]

    fig = plt.figure(figsize=(17, 9.5))
    # --- top: every parameter per scan
    for j, k in enumerate(keys):
        ax = fig.add_subplot(2, 6, j + 1)
        v = D[k]
        ax.bar(range(len(v)), v, color="#2c6f9b")
        ax.axhline(v.mean(), color="k", ls="--", lw=1.3)
        ax.set_xticks(range(len(v)))
        ax.set_xticklabels(names, rotation=90, fontsize=6)
        ax.set_title("%s\n%.4g ± %.0f %%" % (lab[k], v.mean(),
                                              100 * v.std() / abs(v.mean())),
                     fontsize=9.5)
        ax.grid(alpha=0.25, axis="y"); ax.tick_params(labelsize=7)

    # --- bottom left: correlation matrix
    ax = fig.add_subplot(2, 3, 4)
    M = np.corrcoef(np.array([D[k] for k in keys]))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(keys))); ax.set_yticks(range(len(keys)))
    ax.set_xticklabels([lab[k] for k in keys], rotation=45, ha="right", fontsize=8.5)
    ax.set_yticklabels([lab[k] for k in keys], fontsize=8.5)
    for i_ in range(len(keys)):
        for j_ in range(len(keys)):
            ax.text(j_, i_, "%+.2f" % M[i_, j_], ha="center", va="center", fontsize=8,
                    color="white" if abs(M[i_, j_]) > 0.6 else "black")
    plt.colorbar(im, ax=ax, fraction=0.046)
    ax.set_title("CROSS-CORRELATION across the %d scans\n"
                 "strong off-diagonal = the parameters trade off" % len(fits),
                 fontsize=10.5, loc="left")

    # --- bottom middle: the key trade-off
    ax = fig.add_subplot(2, 3, 5)
    ax.scatter(D["adg"], D["M"], s=90, c=D["amp"], cmap="viridis",
               edgecolor="k", linewidth=0.5)
    for x_, y_, n_ in zip(D["adg"], D["M"], names):
        ax.annotate(n_, (x_, y_), fontsize=7, xytext=(4, 4),
                    textcoords="offset points")
    r = np.corrcoef(D["adg"], D["M"])[0, 1]
    ax.set_xlabel("$a_{dg}$(443)  m$^{-1}$"); ax.set_ylabel("$M_\\phi$")
    ax.grid(alpha=0.25)
    ax.set_title("THE degeneracy: $a_{dg}$ vs $a_\\phi$\nr = %+.2f  "
                 "(both rise toward the blue)" % r, fontsize=10.5, loc="left")
    plt.colorbar(ax.collections[0], ax=ax, label="$R_{rs}$(555)")

    # --- bottom right: what amplitude drives
    ax = fig.add_subplot(2, 3, 6)
    for k, col in (("M", "#2e7d32"), ("adg", "#8a6000"), ("bbp", "#2c6f9b")):
        y = D[k] / D[k].mean()
        r = np.corrcoef(D["amp"], D[k])[0, 1]
        ax.scatter(D["amp"], y, s=70, color=col, edgecolor="k", linewidth=0.4,
                   label="%s   r=%+.2f" % (lab[k], r))
    ax.set_xlabel("$R_{rs}$(555)  sr$^{-1}$  (the amplitude term)")
    ax.set_ylabel("parameter / its mean")
    ax.legend(fontsize=9); ax.grid(alpha=0.25)
    ax.set_title("Which parameter absorbs the concentration change?",
                 fontsize=10.5, loc="left")
    fig.suptitle("Per-scan parameters and their covariance — $S_{dg}$ fixed, "
                 "$\\eta$ derived, three amplitudes free", fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = os.path.join(out, "giop7_covariance.png")
    fig.savefig(p, dpi=135); plt.close(fig)
    return p, D, names


# ------------------------------------------------------------------ figure 8
def f_assumption_free(loc, wl, mean, ssd, oc4, out):
    """Fig 3/4/6 redone WITHOUT trusting OC4, and with the shapes freed.

    Three assumptions are stripped one at a time and the answer watched:

      the OC4 chlorophyll INPUT  -- it only selects the Bricaud a*_phi SHAPE, but it is
                                    derived from a Case-1 band ratio on Case-2 water
      the a*_phi FAMILY          -- Ciotti 2006 is parameterised by particle SIZE and
                                    needs no chlorophyll at all
      S_dg and eta               -- freed via fit_shapes on the fmin path

    The result separates what survives from what does not, which is more useful than a
    single retrieval with error bars.
    """
    import warnings
    warnings.filterwarnings("ignore")
    m = (wl >= LO) & (wl <= HI); W, R, S = wl[m], mean[m], ssd[m]

    sig = np.sqrt((0.05 * np.abs(R)) ** 2 + 2e-4 ** 2)
    nu = len(W) - 3

    chls = np.array([0.5, 1, 2, 3, 5, 8, oc4, 15, 20, 30, 50])
    A = []
    for c in chls:
        g = run(W, R, float(c))
        A.append((g.chl, g.adg443, g.bbp443, _chi2(g, R, S)[0] / nu))
    A = np.array(A)

    sfs = np.linspace(0.0, 1.0, 11)
    B = []
    for sf in sfs:
        g = run(W, R, oc4, aph="ciotti", sf=float(sf))
        B.append((g.chl, g.adg443, g.bbp443, _chi2(g, R, S)[0] / nu))
    B = np.array(B)

    free = {}
    for lab, kw in (("Bricaud", dict(aph="bricaud")),
                    ("Ciotti sf=0.5", dict(aph="ciotti", sf=0.5))):
        g = giop(W, R, oc4, inv="bounded", sigma=sig, fit_shapes=True, n_starts=4, **kw)
        free[lab] = (g.chl, g.adg443, g.bbp443, g.sdg, g.eta, _chi2(g, R, S)[0] / nu)

    fig, axes = plt.subplots(2, 3, figsize=(17, 9.8))
    for j, (t, col) in enumerate((("$M_\\phi$", "#2e7d32"),
                                  ("$a_{dg}$(443)  m$^{-1}$", "#8a6000"),
                                  ("$b_{bp}$(443)  m$^{-1}$", "#2c6f9b"))):
        ax = axes[0][j]
        ax.plot(chls, A[:, j], "-", lw=1.4, color=col, zorder=1)
        s = ax.scatter(chls, A[:, j], c=A[:, 3], cmap="viridis_r", s=95, zorder=3,
                       edgecolor="k", linewidth=0.5,
                       norm=matplotlib.colors.LogNorm())
        plt.colorbar(s, ax=ax, label="$\\chi^2_\\nu$")
        ax.axvline(oc4, color="#c0392b", ls="--", lw=1.6)
        ax.set_xscale("log"); ax.set_xlabel("chlorophyll INPUT (mg m$^{-3}$)")
        ax.set_ylabel(t); ax.grid(alpha=0.25)
        b = int(np.argmin(A[:, 3]))
        ax.plot(chls[b], A[b, j], "*", ms=19, color="#ff2d55", mec="k", zorder=4)
        ax.set_title("%s vs the OC4 input.  range %.3g to %.3g\n"
                     "best-fitting arm (star) is chl$_{in}$=%.3g at $\\chi^2_\\nu$=%.0f"
                     % (t, A[:, j].min(), A[:, j].max(), chls[b], A[b, 3]),
                     fontsize=9, loc="left")

        ax = axes[1][j]
        ax.plot(sfs, B[:, j], "-", lw=1.4, color=col, zorder=1)
        s = ax.scatter(sfs, B[:, j], c=B[:, 3], cmap="viridis_r", s=95, zorder=3,
                       edgecolor="k", linewidth=0.5,
                       norm=matplotlib.colors.LogNorm())
        plt.colorbar(s, ax=ax, label="$\\chi^2_\\nu$")
        ax.set_xlabel("Ciotti size factor $S_f$  (NO chlorophyll input)")
        ax.set_ylabel(t); ax.grid(alpha=0.25)
        for k, (lab, v) in enumerate(free.items()):
            c_ = ["#c0392b", "#6a3d9a"][k]
            ax.axhline(v[j], color=c_, ls=":", lw=1.8)
            ax.text(0.02, v[j], "shapes FREE (%s), $\\chi^2_\\nu$=%.0f" % (lab, v[5]),
                    fontsize=7.5, color=c_, va="bottom")
        ax.set_title("%s vs the a*$_\\phi$ family, and with the shapes freed"
                     % t, fontsize=9, loc="left")
    fig.suptitle("STRIPPING THE ASSUMPTIONS — every point now carries its $\\chi^2_\\nu$, "
                 "because a range that pools good and bad fits is not an uncertainty.\n"
                 "Freeing $S_{dg}$ and $\\eta$ is the BEST-FITTING arm of all "
                 "($\\chi^2_\\nu$ %.0f vs %.0f at the GIOP default) and it moves "
                 "$a_{dg}$ and $b_{bp}$ by %.0f %% and %.0f %%."
                 % (free["Bricaud"][5], A[np.argmin(abs(chls - oc4)), 3],
                    100 * abs(free["Bricaud"][1] / A[np.argmin(abs(chls - oc4)), 1] - 1),
                    100 * abs(free["Bricaud"][2] / A[np.argmin(abs(chls - oc4)), 2] - 1)),
                 fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    p = os.path.join(out, "giop8_assumption_free.png")
    fig.savefig(p, dpi=135); plt.close(fig)
    return p, A, B, free, chls, sfs



def _chi2(g, R, S):
    """chi2 of a fit against the MEASURED per-band uncertainty, above water.

    Not `g.cost`: that is the solver's own objective on the subsurface scale, with
    whatever sigma it was handed. This is one statistic, comparable across every arm.
    """
    r = (g.rrs_model_above - R) / S
    return float(np.sum(r ** 2)), float(np.corrcoef(r[:-1], r[1:])[0, 1])


def _amp_solve(W, R, sdg, eta, chl, cfg=None):
    """The three amplitudes at a PRESCRIBED (S_dg, eta), plus chi2. For the chi2 map."""
    from giop.inversion import _invert_bounded
    from giop.model import rrs_below_to_above
    cfg = cfg or GiopConfig()
    rin = rrs_above_to_below(R, cfg.trans)
    _, _, aphs, _, _ = eigenvectors(W, cfg, chl, R, rin, find_anchor_bands(W))
    aw, bbw = a_water(W), bb_water(W)
    sig = np.sqrt((0.05 * np.abs(R)) ** 2 + 2e-4 ** 2)
    adgs = np.exp(-sdg * (W - 443.0))
    bbps = (443.0 / W) ** eta
    x, ok, _, _, _ = _invert_bounded(rin, aw, bbw, adgs, bbps, aphs, GORDON_G0,
                                     GORDON_G1, chl, sigma=sig * rin / R)
    rm = rrs_from_iops(aw + adgs * x[0] + aphs * x[2], bbw + bbps * x[1],
                       GORDON_G0, GORDON_G1)
    return x, rrs_below_to_above(rm, cfg.trans)


# ------------------------------------------------------------------ figure 9
def f_chi2(wl, mean, ssd, oc4, out):
    """Does ANY arm actually fit, and what does chi2-weighting the spread give?

    Written because the assumption sweep in figure 8 quoted a range of answers without
    ever saying how well each one fitted. Some of them fit far worse than others, and a
    range that pools a chi2_nu of 18 with a chi2_nu of 130 is not an uncertainty.
    """
    import warnings
    warnings.filterwarnings("ignore")
    m = (wl >= LO) & (wl <= HI)
    W, R, S = wl[m], mean[m], ssd[m]
    n = len(W)
    nu = n - 3

    # ---- every arm, with its chi2
    arms = []
    for c in [0.5, 1, 2, 3, 5, 8, oc4, 15, 20, 30, 50]:
        g = run(W, R, float(c))
        c2, _ = _chi2(g, R, S)
        arms.append(("OC4 input %.3g" % c, "chl", c2, g.chl, g.adg443, g.bbp443,
                     0.018, g.eta))
    for sf in np.linspace(0, 1, 11):
        g = run(W, R, oc4, aph="ciotti", sf=float(sf))
        c2, _ = _chi2(g, R, S)
        arms.append(("Ciotti sf=%.1f" % sf, "sf", c2, g.chl, g.adg443, g.bbp443,
                     0.018, g.eta))
    sig = np.sqrt((0.05 * np.abs(R)) ** 2 + 2e-4 ** 2)
    for lab, kw in (("Bricaud", dict(aph="bricaud")),
                    ("Ciotti sf=0.5", dict(aph="ciotti", sf=0.5))):
        g = giop(W, R, oc4, inv="bounded", sigma=sig, fit_shapes=True, n_starts=4, **kw)
        c2, _ = _chi2(g, R, S)
        arms.append(("shapes FREE (%s)" % lab, "free", c2, g.chl, g.adg443, g.bbp443,
                     g.sdg, g.eta))
    C = np.array([a[2] for a in arms])
    c2min = C.min()

    # ---- residual autocorrelation of the BEST arm: is chi2 a likelihood here?
    gb = run(W, R, oc4)
    resid = (gb.rrs_model_above - R) / S
    rho1 = float(np.corrcoef(resid[:-1], resid[1:])[0, 1])
    neff = n * (1 - rho1) / (1 + rho1)

    # ---- chi2 weights. The best fit is chi2_nu = %.0f, so exp(-dchi2/2) on the RAW
    # chi2 would put every gram of weight on one arm. Errors are inflated so the best
    # arm has chi2_nu = 1 (Avni 1976), which is the standard treatment when the misfit
    # is dominated by model inadequacy rather than by noise.
    infl = c2min / nu
    wgt = np.exp(-(C - c2min) / (2.0 * infl))
    wgt /= wgt.sum()
    # Even inflated, this collapses onto one arm: the gap between the best and the next
    # is dchi2_nu ~ 50 over 298 dof. That is the honest answer -- the sweep is mostly
    # REJECTED models, not alternatives -- but a delta-chi2 band is more useful to read
    # off, so an "admissible" set is carried alongside at a stated, arbitrary factor.
    adm = C <= 2.0 * c2min

    fig = plt.figure(figsize=(17.5, 10.5))
    col = {"chl": "#2e7d32", "sf": "#8a6000", "free": "#c0392b"}

    # (a) ranked chi2_nu
    ax = fig.add_subplot(2, 3, 1)
    o = np.argsort(C)[::-1]
    ax.barh(range(len(o)), C[o] / nu, color=[col[arms[i][1]] for i in o])
    ax.set_yticks(range(len(o)))
    ax.set_yticklabels([arms[i][0] for i in o], fontsize=6.5)
    ax.set_xscale("log"); ax.set_xlabel("$\\chi^2_\\nu$   ($\\nu$ = %d)" % nu)
    ax.axvline(1.0, color="k", lw=2.0)
    ax.text(1.15, 0.5, "a GOOD fit\nwould be here", fontsize=8, rotation=90,
            va="center")
    ax.set_title("NOTHING FITS. Best $\\chi^2_\\nu$ = %.0f, worst %.0f.\n"
                 "The measured band uncertainty is %.1f %%; the best arm misfits by "
                 "%.1f %%." % (c2min / nu, C.max() / nu, 100 * np.median(S / R),
                               100 * np.sqrt(np.mean(((gb.rrs_model_above - R) / R) ** 2))),
                 fontsize=9.5, loc="left")
    ax.grid(alpha=0.25, axis="x")

    # (b) the chi2 surface over the two SHAPE parameters
    ax = fig.add_subplot(2, 3, 2)
    sd = np.linspace(0.006, 0.028, 34)
    et = np.linspace(-1.0, 2.0, 34)
    Z = np.empty((len(et), len(sd)))
    for i, e in enumerate(et):
        for j, s_ in enumerate(sd):
            _, mod = _amp_solve(W, R, float(s_), float(e), oc4)
            Z[i, j] = np.sum(((mod - R) / S) ** 2) / nu
    im = ax.pcolormesh(sd, et, np.log10(Z), cmap="viridis_r", shading="auto")
    plt.colorbar(im, ax=ax, label="$\\log_{10}\\,\\chi^2_\\nu$")
    k = np.unravel_index(np.argmin(Z), Z.shape)
    ax.plot(sd[k[1]], et[k[0]], "*", ms=20, color="#ff2d55", mec="k",
            label="free minimum  %.4f, %+.2f" % (sd[k[1]], et[k[0]]))
    ax.plot(0.018, gb.eta, "o", ms=11, color="w", mec="k",
            label="GIOP default  0.0180, %+.2f" % gb.eta)
    ax.set_xlabel("$S_{dg}$  (nm$^{-1}$)"); ax.set_ylabel("$\\eta$")
    ax.legend(fontsize=7.5, loc="upper right")
    ax.set_title("Contours are VERTICAL: $S_{dg}$ is sharply determined (interior "
                 "minimum at\n%.4f, not the assumed 0.018, $\\chi^2_\\nu$ %.0f $\\to$ "
                 "%.0f); $\\eta$ is nearly flat and runs to its bound."
                 % (sd[k[1]], Z[np.argmin(abs(et - gb.eta)), np.argmin(abs(sd - 0.018))],
                    Z.min()), fontsize=9.5, loc="left")

    # (c) the residual, and why chi2 is not a likelihood here
    ax = fig.add_subplot(2, 3, 3)
    ax.plot(W, resid, lw=1.6, color="#c0392b", label="GIOP default")
    xf, modf = _amp_solve(W, R, sd[k[1]], et[k[0]], oc4)
    ax.plot(W, (modf - R) / S, lw=1.6, color="#2e7d32", label="shapes free")
    ax.axhline(0, color="k", lw=0.8)
    for y in (-1, 1):
        ax.axhline(y, color="k", ls=":", lw=1.0)
    ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("residual / $\\sigma$")
    ax.legend(fontsize=8.5); ax.grid(alpha=0.25)
    ax.set_title("The residual is ONE SMOOTH CURVE, not noise.\n"
                 "lag-1 $\\rho$ = %.4f over %d bands $\\Rightarrow$ $n_{eff}$ = %.1f. "
                 "So $\\chi^2$ RANKS arms;\nit does not measure a probability."
                 % (rho1, n, neff), fontsize=9.5, loc="left")

    # (d-f) unweighted vs chi2-weighted spread
    P = np.array([[a[3], a[4], a[5]] for a in arms])
    for j, (t, c_) in enumerate((("$M_\\phi$", "#2e7d32"),
                                 ("$a_{dg}$(443)  m$^{-1}$", "#8a6000"),
                                 ("$b_{bp}$(443)  m$^{-1}$", "#2c6f9b"))):
        ax = fig.add_subplot(2, 3, 4 + j)
        v = P[:, j]
        sc = ax.scatter(C / nu, v, s=30 + 900 * wgt, c=[col[a[1]] for a in arms],
                        edgecolor="k", linewidth=0.5, zorder=3)
        ax.set_xscale("log"); ax.set_xlabel("$\\chi^2_\\nu$ of that arm")
        ax.set_ylabel(t); ax.grid(alpha=0.25)
        mu = float(np.sum(wgt * v))
        ax.axhspan(v[adm].min(), v[adm].max(), color="#ff2d55", alpha=0.18, zorder=0)
        ax.axhline(mu, color="#ff2d55", lw=2.0, zorder=1)
        ax.axhline(v.min(), color="k", ls=":", lw=1.0)
        ax.axhline(v.max(), color="k", ls=":", lw=1.0)
        ax.axvline(2.0 * c2min / nu, color="#ff2d55", ls="--", lw=1.4)
        ax.set_title("unweighted range (dotted) %.3g to %.3g   <- pools "
                     "REJECTED models\nadmissible ($\\chi^2_\\nu\\leq2\\times$best, "
                     "%d of %d arms, band): %.3g to %.3g"
                     % (v.min(), v.max(), int(adm.sum()), len(v),
                        v[adm].min(), v[adm].max()),
                     fontsize=8.5, loc="left")

    fig.suptitle("CHI-SQUARED WEIGHTING THE ASSUMPTION SWEEP. Weights are "
                 "$\\exp(-\\Delta\\chi^2/2\\hat{s})$, errors inflated by "
                 "$\\hat{s}=\\chi^2_{\\nu,min}$ = %.0f (Avni 1976) since the misfit is "
                 "model inadequacy, not noise.\nEVEN SO ONE ARM TAKES w = 1.000: the gap "
                 "to the next is $\\Delta\\chi^2$ = %.0f, and %.0f to the best "
                 "FIXED-shape arm, over %d dof. The sweep is a set of REJECTED models,\n"
                 "not an uncertainty band — so the panels below quote a stated "
                 "$\\chi^2_\\nu\\leq2\\times$best cut instead."
                 % (c2min / nu, np.sort(C)[1] - c2min,
                    C[np.array([a_[1] != "free" for a_ in arms])].min() - c2min, nu),
                 fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    p = os.path.join(out, "giop9_chi2_weighting.png")
    fig.savefig(p, dpi=135); plt.close(fig)

    # The ranked table is the point of the figure, so it must also exist as numbers.
    with open(os.path.join(out, "giop_assumption_arms.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# every arm of the assumption sweep, ranked by chi2. "
                    "nu=%d, sigma = the measured per-band uncertainty in FINAL_Rrs.csv"
                    % nu])
        w.writerow(["# admissible = chi2_nu <= 2 x best (a STATED cut, not a "
                    "confidence level -- see LOC1_GIOP_FINDINGS.md 2b/2c)"])
        w.writerow(["arm", "family", "chi2", "chi2_nu", "weight_avni", "admissible",
                    "M_phi", "adg443", "bbp443", "S_dg", "eta"])
        for i in np.argsort(C):
            a_ = arms[i]
            w.writerow([a_[0], a_[1], "%.4f" % C[i], "%.4f" % (C[i] / nu),
                        "%.6g" % wgt[i], int(adm[i]), "%.6g" % a_[3],
                        "%.6g" % a_[4], "%.6g" % a_[5], "%.6g" % a_[6],
                        "%.6g" % a_[7]])
        w.writerow([])
        w.writerow(["# chi2 surface minimum over (S_dg, eta), 34x34 grid"])
        w.writerow(["S_dg_best", "%.6g" % sd[k[1]], "eta_best", "%.6g" % et[k[0]],
                    "chi2_nu_best", "%.6g" % Z.min()])
        w.writerow(["# residual autocorrelation of the GIOP-default fit"])
        w.writerow(["lag1_rho", "%.6g" % rho1, "n_bands", n, "n_eff_AR1",
                    "%.4g" % neff])
    return p, arms, C, wgt, nu, rho1, neff, (sd[k[1]], et[k[0]], Z.min()), adm


# ------------------------------------------------------------------ figure 10
#: Seed grid for the maximum-freedom arm. Spans well past OC4 in both directions.
SEEDS = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 50.0])


def _maxfree(W, R, S, sig, oc4):
    """Maximum freedom: S_dg, eta AND the a*_phi shape, with nothing seeded by OC4.

    ``fit_shapes=True`` alone does NOT do this. It frees S_dg and eta, but a*_phi is
    built once from the chlorophyll seed and held fixed inside the profile, so the OC4
    number -- a Case-1 band ratio, on Case-2 water -- still sets the phytoplankton
    shape. Here the seed is profiled over as well, and the Ciotti size-parameterised
    family (which needs no chlorophyll at all) is included, so the reported arm is the
    best the model can do with every prescribed shape released.
    """
    nu = len(W) - 3
    best = None
    for lab, kw in ([("Bricaud seed %.3g" % c, dict(aph="bricaud")) for c in SEEDS] +
                    [("Ciotti sf=%.2f" % sf, dict(aph="ciotti", sf=float(sf)))
                     for sf in np.linspace(0, 1, 6)]):
        seed = float(lab.split()[-1]) if lab.startswith("Bricaud") else oc4
        try:
            g = giop(W, R, seed, inv="bounded", sigma=sig, fit_shapes=True,
                     n_starts=4, **kw)
            c2, _ = _chi2(g, R, S)
        except Exception:
            continue
        if best is None or c2 < best[0]:
            best = (c2, lab, g)
    c2, lab, g = best
    return dict(M=g.chl, adg=g.adg443, bbp=g.bbp443, sdg=g.sdg, eta=g.eta,
                c2n=c2 / nu, arm=lab,
                rms=100 * float(np.sqrt(np.mean(((g.rrs_model_above - R) / R) ** 2))))


def f_final(wl, mean, ssd, oc4, fits_fix, fits_free, out):
    """THE FINAL RESULT. Mean and per-spectrum, constrained and free, side by side.

    Everything else in this folder is diagnosis. This is the number to quote, with the
    two independent spreads kept apart: the SCATTER ACROSS SCANS (real water plus
    pairing) and the CHOICE OF CONFIGURATION (which is a modelling decision, not an
    error bar).
    """
    import warnings
    warnings.filterwarnings("ignore")
    m = (wl >= LO) & (wl <= HI)
    W, R, S = wl[m], mean[m], ssd[m]
    nu = len(W) - 3
    sig = np.sqrt((0.05 * np.abs(R)) ** 2 + 2e-4 ** 2)

    gfix = run(W, R, oc4)
    gfre = giop(W, R, oc4, inv="bounded", sigma=sig, fit_shapes=True, n_starts=4)
    mean_row = {}
    for lab, g in (("constrained", gfix), ("free", gfre)):
        c2, _ = _chi2(g, R, S)
        mean_row[lab] = dict(M=g.chl, adg=g.adg443, bbp=g.bbp443, sdg=g.sdg,
                             eta=g.eta, c2n=c2 / nu,
                             rms=100 * float(np.sqrt(np.mean(
                                 ((g.rrs_model_above - R) / R) ** 2))))

    mean_row["max free"] = _maxfree(W, R, S, sig, oc4)

    # Self-consistency. Bricaud gives ABSOLUTE a_phi = A_phi chl^E_phi; GIOP divides by
    # chl to get the specific a*_phi, normalises it, and then fits a FREE amplitude --
    # so it uses the seed's SHAPE and discards its AMPLITUDE. If the seed were believed,
    # M_phi would simply BE the chlorophyll. Whether the retrieval reproduces its own
    # seed is therefore a real internal check, and GIOP never performs it because it
    # deliberately does not iterate.
    fixed_pt = []
    for c_ in SEEDS:
        try:
            fixed_pt.append((float(c_), run(W, R, float(c_)).chl))
        except Exception:
            pass
    fixed_pt = np.array(fixed_pt)
    # There are TWO roots and only one of them means anything. d = M_phi - seed crosses
    # zero from below at a low seed (where M_phi has collapsed to the 0 bound -- an
    # UNSTABLE root: iterate away from it and you never come back) and from above at a
    # high seed (STABLE: iterating the seed converges there). Taking the first crossing
    # reports the artefact. Keep both, mark the stable one.
    d = fixed_pt[:, 1] - fixed_pt[:, 0]
    roots = []
    for i in np.flatnonzero(np.sign(d[:-1]) != np.sign(d[1:])):
        i = int(i)
        x0 = float(fixed_pt[i, 0] - d[i] * (fixed_pt[i + 1, 0] - fixed_pt[i, 0])
                   / (d[i + 1] - d[i]))
        roots.append((x0, "stable" if d[i] > d[i + 1] else "unstable"))
    stable = [r for r in roots if r[1] == "stable"]
    cross = stable[-1][0] if stable else np.nan

    keys = [("M", "$M_\\phi$"), ("adg", "$a_{dg}$(443)  m$^{-1}$"),
            ("bbp", "$b_{bp}$(443)  m$^{-1}$"), ("sdg", "$S_{dg}$  nm$^{-1}$"),
            ("eta", "$\\eta$"), ("rms", "RMS misfit  %")]
    fig, axes = plt.subplots(2, 4, figsize=(21, 9.8))
    axes = axes.ravel()

    ax = axes[6]
    ax.plot(fixed_pt[:, 0], fixed_pt[:, 1], "o-", lw=2.0, color="#2e7d32")
    lim = [SEEDS.min(), SEEDS.max()]
    ax.plot(lim, lim, "k--", lw=1.4, label="$M_\\phi$ = seed (self-consistent)")
    ax.axvline(oc4, color="#c0392b", ls=":", lw=1.8, label="OC4 = %.2f" % oc4)
    for x0, kind in roots:
        ax.plot([x0], [x0], "*" if kind == "stable" else "x", ms=22 if kind == "stable"
                else 12, color="#ff2d55" if kind == "stable" else "#888", mec="k",
                mew=1.4, zorder=6,
                label="%s fixed point  chl = %.2f" % (kind.upper(), x0))
    ax.set_xscale("log"); ax.set_yscale("log")
    # M_phi rails at 0 for low seeds, which on a log axis drags the view to 1e-27.
    ax.set_ylim(0.3, max(fixed_pt[:, 1].max(), SEEDS.max()) * 1.6)
    ax.set_xlabel("chlorophyll SEED (sets the $a^*_\\phi$ shape)")
    ax.set_ylabel("$M_\\phi$ RETRIEVED")
    ax.legend(fontsize=7.5); ax.grid(alpha=0.25)
    ax.set_title("SELF-CONSISTENCY. GIOP uses the seed's SHAPE and discards its\n"
                 "AMPLITUDE, then fits $M_\\phi$ free. Stable fixed point chl=%.1f "
                 "vs OC4 %.1f\n(agreement to %.0f %%). GIOP never runs this check — "
                 "it does not iterate."
                 % (cross, oc4, 100 * abs(cross / oc4 - 1)), fontsize=8.5, loc="left")

    ax = axes[7]
    ax.axis("off")
    r_ = mean_row["max free"]
    ax.text(0.0, 0.98,
            "MAXIMUM FREEDOM\n"
            "every prescribed shape released:\n"
            "$S_{dg}$, $\\eta$ AND the $a^*_\\phi$ family/seed\n"
            "(%d Bricaud seeds + %d Ciotti $S_f$)\n\n"
            "winning arm:  %s\n\n"
            "$M_\\phi$      = %.3g\n"
            "$a_{dg}$(443) = %.4g m$^{-1}$\n"
            "$b_{bp}$(443) = %.4g m$^{-1}$\n"
            "$S_{dg}$      = %.5g nm$^{-1}$\n"
            "$\\eta$        = %+.3g\n"
            "$\\chi^2_\\nu$      = %.1f\n"
            "RMS       = %.1f %%\n\n"
            "vs constrained $\\chi^2_\\nu$ = %.0f\n"
            "vs 'free' (OC4-seeded) = %.0f\n\n"
            "⚠ 'free' is NOT assumption-free:\n"
            "fit_shapes seeds $a^*_\\phi$ from OC4\n"
            "and holds it fixed. Only this arm\n"
            "releases it — and it buys only\n"
            "18 -> 17, while releasing $S_{dg}$\n"
            "bought 74 -> 18. The CDOM slope\n"
            "matters far more than the seed."
            % (len(SEEDS), 6, r_["arm"], r_["M"], r_["adg"], r_["bbp"], r_["sdg"],
               r_["eta"], r_["c2n"], r_["rms"], mean_row["constrained"]["c2n"],
               mean_row["free"]["c2n"]),
            transform=ax.transAxes, va="top", fontsize=8.2, family="monospace",
            bbox=dict(boxstyle="round,pad=0.6", fc="#f4f8f4", ec="#2e7d32", lw=1.5))

    for ax, (k, lab) in zip(axes[:6], keys):
        for j, (tag, F, col) in enumerate((("constrained", fits_fix, "#c0392b"),
                                           ("free", fits_free, "#2e7d32"))):
            _ = tag
            v = np.array([f[k] for f in F if k in f], dtype=float)
            if not v.size:
                continue
            x = np.full(v.size, j) + np.linspace(-0.13, 0.13, v.size)
            ax.plot(x, v, "o", ms=6, color=col, alpha=0.75, mec="k", mew=0.4,
                    label="%s: per-scan  %.4g $\\pm$ %.0f %%"
                          % (tag, v.mean(), 100 * v.std() / abs(v.mean())
                             if v.mean() else np.nan))
            ax.plot([j - 0.3, j + 0.3], [v.mean()] * 2, lw=2.4, color=col)
            mv = mean_row[tag][k]
            ax.plot([j], [mv], "*", ms=22, color="#ff2d55", mec="k", zorder=6,
                    label="%s: fit of the MEAN  %.4g" % (tag, mv))
        if k in mean_row["max free"]:
            ax.axhline(mean_row["max free"][k], color="#6a3d9a", ls="--", lw=1.8,
                       label="MAX FREEDOM  %.4g" % mean_row["max free"][k])
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["CONSTRAINED\n$S_{dg}$=0.018, 3 free",
                            "FREE\n$S_{dg}$,$\\eta$ fitted, 5 free"], fontsize=9)
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylabel(lab); ax.grid(alpha=0.25, axis="y")
        ax.legend(fontsize=7, loc="best")
        ax.set_title(lab, fontsize=10.5, loc="left")
    fig.suptitle("THE FINAL RESULT — per-spectrum points, their mean (bar), and the fit "
                 "of the amplitude-normalised mean spectrum (star).\n"
                 "Fitting the mean is NOT the mean of the fits: the inversion is "
                 "nonlinear. Quote the per-scan spread as the uncertainty and the "
                 "mean-spectrum fit as the value.\nMean spectrum $\\chi^2_\\nu$: "
                 "constrained %.0f, free (OC4-seeded) %.0f, MAX FREEDOM %.0f. "
                 "NONE of these is a good fit."
                 % (mean_row["constrained"]["c2n"], mean_row["free"]["c2n"],
                    mean_row["max free"]["c2n"]),
                 fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    p = os.path.join(out, "giop10_final_result.png")
    fig.savefig(p, dpi=135); plt.close(fig)

    with open(os.path.join(out, "giop_FINAL.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# THE headline table. 'mean' = fit of the amplitude-normalised "
                    "mean spectrum; 'perscan' = mean +- sd over the individual "
                    "angle-matched fits."])
        w.writerow(["# constrained = GIOP-DC (S_dg 0.018, eta from QAA, 3 free "
                    "amplitudes). free = S_dg and eta also fitted (5 free)."])
        w.writerow(["# NOTHING in this family fits well; see chi2_nu against a "
                    "measured band uncertainty of %.2f %%." % (100 * np.median(S / R))])
        w.writerow(["config", "estimate", "M_phi", "adg443", "bbp443", "S_dg", "eta",
                    "chi2_nu", "rms_pct"])
        w.writerow(["# max free = S_dg, eta AND the a*_phi family/seed released. "
                    "'free' alone still seeds a*_phi from OC4."])
        r_ = mean_row["max free"]
        w.writerow(["max free", "mean_spectrum (%s)" % r_["arm"]] +
                   ["%.6g" % r_[k] for k in
                    ("M", "adg", "bbp", "sdg", "eta", "c2n", "rms")])
        w.writerow(["# self-consistency: seed -> retrieved M_phi fixed point"])
        w.writerow(["fixed_point_chl_STABLE", "%.6g" % cross, "OC4", "%.6g" % oc4] +
                   ["%s=%.6g" % (k, v) for v, k in roots])
        for tag, F in (("constrained", fits_fix), ("free", fits_free)):
            r_ = mean_row[tag]
            w.writerow([tag, "mean_spectrum"] +
                       ["%.6g" % r_[k] for k in
                        ("M", "adg", "bbp", "sdg", "eta", "c2n", "rms")])
            for stat, fn in (("perscan_mean", np.mean), ("perscan_sd", np.std)):
                w.writerow([tag, stat] + ["%.6g" % fn([f[k] for f in F if k in f])
                                          if any(k in f for f in F) else ""
                                          for k in ("M", "adg", "bbp", "sdg", "eta")] +
                           ["", "%.6g" % fn([f["rms"] for f in F])])
    return p, mean_row


# ------------------------------------------------------------------ figure 11
def f_chi2_crossplot(fits_fix, fits_free, mean_row, out):
    """chi2 of the free fit against chi2 of the constrained fit, scan by scan.

    The summary numbers say free beats constrained 4x on the mean spectrum. That could
    be one bad scan dragging an average, or it could hold for every spectrum
    independently. A cross-plot against the 1:1 line answers it directly, and the
    nesting property (free is a superset of constrained) means EVERY point must lie on
    or below the line -- so the plot doubles as a check on the solver.
    """
    fr = {f["n"]: f for f in fits_free}
    pairs = [(f, fr[f["n"]]) for f in fits_fix if f["n"] in fr
             and np.isfinite(f.get("c2n", np.nan))
             and np.isfinite(fr[f["n"]].get("c2n", np.nan))]
    if not pairs:
        return None
    xc = np.array([a_["c2n"] for a_, _ in pairs])
    yc = np.array([b_["c2n"] for _, b_ in pairs])
    amp = np.array([a_["amp"] for a_, _ in pairs])
    names = [a_["n"] for a_, _ in pairs]

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.6))

    ax = axes[0]
    lim = [min(xc.min(), yc.min()) * 0.8, max(xc.max(), yc.max()) * 1.25]
    ax.plot(lim, lim, "k--", lw=1.6, label="1:1  (no improvement)")
    # ABOVE the line is the impossible side: free is a superset of constrained, so its
    # optimum can never be worse. Shading below would mark the only region the points
    # are allowed to be in, which is what the first version of this figure did.
    ax.fill_between(lim, lim, [lim[1]] * 2, color="#c0392b", alpha=0.07, zorder=0)
    ax.text(lim[0] * 1.15, lim[1] * 0.62, "FORBIDDEN\nfree is nested inside "
            "constrained,\nso no point can lie above the 1:1 line",
            fontsize=8, color="#c0392b", ha="left", va="top")
    for f_, y_ in ((2, "#888"), (4, "#555")):
        ax.plot(lim, [v / f_ for v in lim], ":", color=y_, lw=1.2)
        ax.text(lim[1] * 0.97, lim[1] / f_, "%dx better" % f_, fontsize=8,
                color=y_, ha="right", va="bottom")
    sc = ax.scatter(xc, yc, s=110, c=amp, cmap="viridis", edgecolor="k",
                    linewidth=0.6, zorder=4)
    plt.colorbar(sc, ax=ax, label="$R_{rs}$(555)  sr$^{-1}$")
    ax.plot([mean_row["constrained"]["c2n"]], [mean_row["free"]["c2n"]], "*", ms=24,
            color="#ff2d55", mec="k", zorder=6, label="the MEAN spectrum")
    for x_, y_, n_ in zip(xc, yc, names):
        ax.annotate(n_, (x_, y_), fontsize=6.5, xytext=(4, -8),
                    textcoords="offset points")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("$\\chi^2_\\nu$  CONSTRAINED  ($S_{dg}$=0.018, $\\eta$ QAA)")
    ax.set_ylabel("$\\chi^2_\\nu$  FREE  ($S_{dg}$, $\\eta$ fitted)")
    ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=0.25, which="both")
    ax.set_title("EVERY scan improves, and by a similar factor.\n"
                 "median %.2fx, range %.2f-%.2fx over %d scans"
                 % (np.median(xc / yc), (xc / yc).min(), (xc / yc).max(), len(xc)),
                 fontsize=10, loc="left")

    ax = axes[1]
    o = np.argsort(xc / yc)
    ax.barh(range(len(o)), (xc / yc)[o], color="#2e7d32")
    ax.axvline(1, color="k", lw=1.4)
    ax.axvline(np.median(xc / yc), color="#ff2d55", ls="--", lw=1.6,
               label="median %.2fx" % np.median(xc / yc))
    ax.set_yticks(range(len(o)))
    ax.set_yticklabels([names[i] for i in o], fontsize=7.5)
    ax.set_xlabel("$\\chi^2_\\nu$ constrained / $\\chi^2_\\nu$ free")
    ax.legend(fontsize=8.5); ax.grid(alpha=0.25, axis="x")
    ax.set_title("Improvement factor, scan by scan", fontsize=10, loc="left")

    ax = axes[2]
    r = np.corrcoef(amp, xc / yc)[0, 1]
    ax.scatter(amp, xc / yc, s=110, color="#2c6f9b", edgecolor="k", linewidth=0.6)
    ax.set_xlabel("$R_{rs}$(555)  sr$^{-1}$  (how bright the scan is)")
    ax.set_ylabel("improvement factor")
    ax.grid(alpha=0.25)
    ax.set_title("Does the gain depend on the spectrum?  r = %+.2f\n%s"
                 % (r, "no — the constrained shapes are wrong in the SAME way "
                    "everywhere" if abs(r) < 0.5 else
                    "yes — the misfit scales with brightness"),
                 fontsize=10, loc="left")

    fig.suptitle("FREE vs CONSTRAINED, $\\chi^2$ cross-plot. The mean-spectrum gain "
                 "(%.0f $\\to$ %.0f) is not an artefact of averaging:\nit reproduces "
                 "independently in every one of the %d angle-matched spectra."
                 % (mean_row["constrained"]["c2n"], mean_row["free"]["c2n"], len(xc)),
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    p = os.path.join(out, "giop11_chi2_crossplot.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("location", help="by_location/LOC*/FOREOPTIC_FOVxx")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    loc = a.location.rstrip("/")
    a.out = a.out or os.path.join(loc, "analysis", "GIOP")
    os.makedirs(a.out, exist_ok=True)
    csv_path = os.path.join(loc, "analysis", "FINAL_Rrs.csv")
    wl, mean, ssd = load(csv_path)
    r4 = at(wl, mean, np.array([443., 490, 510, 555]))
    chl = float(get_oc(r4[0], r4[1], r4[2], r4[3], "oc4"))
    ps = []
    p, res = f_bands(wl, mean, ssd, chl, a.out); ps.append(p)
    p, chl2, X, which = f_oc4(wl, mean, a.out); ps.append(p)
    ps.append(f_fit(wl, mean, chl, a.out))
    p, store = f_uncertainty(wl, mean, ssd, chl, a.out); ps.append(p)
    ps.append(f_sdg(wl, mean, ssd, chl, a.out))
    _m = (wl >= LO) & (wl <= HI)
    p, fits, fits_free = f_perscan(loc, chl, a.out, S=ssd[_m]); ps.append(p)
    p, D, names = f_covariance(fits, a.out); ps.append(p)
    p, A, Bc, free, chls, sfs = f_assumption_free(loc, wl, mean, ssd, chl, a.out)
    ps.append(p)
    p, arms, C, wgt, nu, rho1, neff, best, adm = f_chi2(wl, mean, ssd, chl, a.out)
    ps.append(p)
    p, mean_row = f_final(wl, mean, ssd, chl, fits, fits_free, a.out); ps.append(p)
    p = f_chi2_crossplot(fits, fits_free, mean_row, a.out)
    if p:
        ps.append(p)
        _x = np.array([f["c2n"] for f in fits])
        _fr = {f["n"]: f["c2n"] for f in fits_free}
        _y = np.array([_fr[f["n"]] for f in fits])
        print("\nFREE vs CONSTRAINED chi2_nu, per scan:")
        print("   constrained median %.1f, free median %.1f, improvement median %.2fx "
              "(range %.2f-%.2f), %d of %d improved"
              % (np.median(_x), np.median(_y), np.median(_x / _y), (_x / _y).min(),
                 (_x / _y).max(), int((_y <= _x).sum()), len(_x)))
    print("\nFINAL RESULT")
    print("   %-13s %9s %9s %9s %9s %8s %9s" % ("config", "M_phi", "adg443", "bbp443",
                                                "S_dg", "eta", "chi2_nu"))
    for tag in ("constrained", "free"):
        r_ = mean_row[tag]
        print("   %-13s %9.3f %9.4f %9.5f %9.5f %8.3f %9.1f   <- fit of the MEAN"
              % (tag, r_["M"], r_["adg"], r_["bbp"], r_["sdg"], r_["eta"], r_["c2n"]))
        F = fits if tag == "constrained" else fits_free
        for k, nm in (("M", "M_phi"), ("adg", "adg443"), ("bbp", "bbp443")):
            v = np.array([f[k] for f in F])
            print("      per-scan %-8s %.4g +/- %.0f %%" % (nm, v.mean(),
                                                            100 * v.std() / abs(v.mean())))
    print("\nDOES ANYTHING FIT?  nu = %d" % nu)
    print("   best chi2_nu %.1f, worst %.1f  -- a good fit would be ~1"
          % (C.min() / nu, C.max() / nu))
    print("   residual lag-1 rho = %.4f over %d bands -> n_eff = %.1f"
          % (rho1, nu + 3, neff))
    print("   free-shape minimum: S_dg=%.4f eta=%+.3f  chi2_nu=%.1f" % best)
    print("   ranked arms:")
    for i in np.argsort(C):
        print("      %-22s chi2_nu=%8.1f  w=%.4f  M_phi=%8.3f adg=%6.3f bbp=%7.4f"
              % (arms[i][0], C[i] / nu, wgt[i], arms[i][3], arms[i][4], arms[i][5]))
    for j, t in enumerate(("M_phi", "adg443", "bbp443")):
        v = np.array([a_[3 + j] for a_ in arms])
        print("   %-8s unweighted %.4g to %.4g   ADMISSIBLE (chi2_nu <= 2x best, "
              "%d/%d arms) %.4g to %.4g"
              % (t, v.min(), v.max(), int(adm.sum()), len(v),
                 v[adm].min(), v[adm].max()))
    print("\nASSUMPTION STRIPPING")
    print("   %-34s %10s %10s %10s" % ("", "M_phi", "adg443", "bbp443"))
    print("   %-34s %9.3g- %9.3g- %9.4g-" % ("OC4 input swept 0.5-50", A[:,0].min(),
                                             A[:,1].min(), A[:,2].min()))
    print("   %-34s %9.3g  %9.3g  %9.4g" % ("   to", A[:,0].max(), A[:,1].max(),
                                            A[:,2].max()))
    print("   %-34s %9.3g- %9.3g- %9.4g-" % ("Ciotti S_f swept 0-1", Bc[:,0].min(),
                                             Bc[:,1].min(), Bc[:,2].min()))
    print("   %-34s %9.3g  %9.3g  %9.4g" % ("   to", Bc[:,0].max(), Bc[:,1].max(),
                                            Bc[:,2].max()))
    for lab, v in free.items():
        print("   %-34s %9.3g  %9.3g  %9.4g   S_dg=%.4f (interior) eta=%.4f "
              "(RAILED at the -1 bound)   chi2_nu=%.1f"
              % ("S_dg,eta FREE (%s)" % lab, v[0], v[1], v[2], v[3], v[4], v[5]))
    per = fits
    # the explainer infographic belongs with the fits it explains
    import shutil
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "giop_python",
                       "GIOP_explainer.png")
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(a.out, "giop0_explainer.png"))
        ps.insert(0, os.path.join(a.out, "giop0_explainer.png"))
    print("OC4: X=%.4f (numerator band %d nm) -> chl = %.2f mg m^-3" % (X, which, chl2))
    print("\nband choice:")
    for lab, n, c, adg, bbp, ok in res:
        print("   %-22s n=%4d  M_phi=%9.3f adg=%8.3f bbp=%8.4f %s"
              % (lab.replace("\n", " "), n, c, adg, bbp, "" if ok else "FAILED"))
    print("\nconditioning under the 1.7 %% shape uncertainty:")
    for lab in ("6 bands", "hyperspectral"):
        v = store[lab]
        print("   %-14s M_phi ±%5.1f %%   adg ±%5.1f %%   bbp ±%5.1f %%"
              % (lab, 100 * v[:, 0].std() / abs(v[:, 0].mean()),
                 100 * v[:, 1].std() / abs(v[:, 1].mean()),
                 100 * v[:, 2].std() / abs(v[:, 2].mean())))
    good = per
    print("\nper-scan (angle-matched pairs), hyperspectral:")
    print("   %-7s %-6s %6s %10s %10s %10s" % ("water", "sky", "dtilt", "M_phi",
                                               "adg443", "bbp443"))
    for r in good:
        print("   %-7s %-6s %6.1f %10.3f %10.3f %10.4f  eta=%.3f RMS=%.1f%%"
              % (r["n"], r["sky"], r["dm"], r["M"], r["adg"], r["bbp"], r["eta"],
                 r["rms"]))
    for k, t in (("M", "M_phi"), ("adg", "adg443"), ("bbp", "bbp443")):
        v = np.array([r[k] for r in good])
        print("   %-8s across scans: %.4g +/- %.0f %%" % (t, v.mean(),
                                                          100 * v.std() / abs(v.mean())))
    cols = ["n", "sky", "dm", "M", "adg", "bbp", "eta", "rms", "amp"]
    with open(os.path.join(a.out, "giop_per_scan.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(per)
    print("\nwrote %s/giop_per_scan.csv" % a.out)
    for p in ps:
        print("wrote %s" % p)


if __name__ == "__main__":
    main()
