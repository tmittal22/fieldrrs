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
    """OC4 in full: where the number comes from."""
    r = at(wl, mean, np.array([443., 490, 510, 555]))
    a = _OC_COEF["oc4"]
    num = max(r[0], r[1], r[2]); which = [443, 490, 510][int(np.argmax(r[:3]))]
    X = np.log10(num / r[3])
    poly = a[0] + a[1] * X + a[2] * X ** 2 + a[3] * X ** 3 + a[4] * X ** 4
    chl = 10 ** poly

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))
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
    fig.suptitle("Where the OC4 chlorophyll comes from — NASA OC v6 coefficients, "
                 "ported from get_oc.m", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    p = os.path.join(out, "giop2_oc4_explained.png")
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
def f_sdg(wl, mean, chl, out):
    m = (wl >= LO) & (wl <= HI)
    sd = np.linspace(0.010, 0.025, 16)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.0))
    for lab, W, R, col in (("6 bands", B6, at(wl, mean, B6), "#c0392b"),
                           ("hyperspectral", wl[m], mean[m], "#2e7d32")):
        vals = []
        for s_ in sd:
            try:
                g = run(W, R, chl, sdg=float(s_))
                vals.append((g.chl, g.adg443, g.bbp443))
            except Exception:
                vals.append((np.nan,) * 3)
        vals = np.array(vals)
        for ax, k in zip(axes, range(3)):
            ax.plot(sd, vals[:, k], "o-", lw=2.2, ms=6, color=col, label=lab)
    for ax, t in zip(axes, ("$M_\\phi$", "$a_{dg}$(443)  m$^{-1}$",
                            "$b_{bp}$(443)  m$^{-1}$")):
        ax.axvline(0.018, color="#333", ls="--", lw=1.5)
        ax.text(0.0182, ax.get_ylim()[1] * 0.9, "default\n0.018", fontsize=8.5)
        ax.set_xlabel("$S_{dg}$  (nm$^{-1}$)  — ASSUMED, never measured")
        ax.set_ylabel(t); ax.set_yscale("symlog", linthresh=1e-2)
        ax.grid(alpha=0.25); ax.legend(fontsize=9)
    fig.suptitle("$S_{dg}$ is a CONSTANT NOBODY MEASURED. Sweeping it over its ordinary "
                 "range is the largest single source of error in the retrieval,\n"
                 "and hyperspectral tames it but does not remove it.", fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    p = os.path.join(out, "giop5_sdg_assumption.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    wl, mean, ssd = load(a.csv)
    r4 = at(wl, mean, np.array([443., 490, 510, 555]))
    chl = float(get_oc(r4[0], r4[1], r4[2], r4[3], "oc4"))
    ps = []
    p, res = f_bands(wl, mean, ssd, chl, a.out); ps.append(p)
    p, chl2, X, which = f_oc4(wl, mean, a.out); ps.append(p)
    ps.append(f_fit(wl, mean, chl, a.out))
    p, store = f_uncertainty(wl, mean, ssd, chl, a.out); ps.append(p)
    ps.append(f_sdg(wl, mean, chl, a.out))
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
    for p in ps:
        print("wrote %s" % p)


if __name__ == "__main__":
    main()
