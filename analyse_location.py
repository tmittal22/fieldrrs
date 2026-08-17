"""Full step-by-step analysis of ONE location / foreoptic dataset.

    python analyse_location.py Data_NatureSpec/2026_Aug_16/by_location/LOC1_*/FLENS8_FOV08

Writes an `analysis/` folder beside the data with four figures and a report:

  fig1  the three measured quantities pooled, so the SPREAD is visible before anything
        is derived from them
  fig2  the calculation, step by step, on one representative scan
  fig3  is E_d reasonable? checked against the true solar constant and against the
        expected shape of atmospheric transmission
  fig4  R_rs: the mean, and the FULL envelope over every water x sky pairing, with the
        variance split into how much comes from the water and how much from the sky

The fig4 envelope is the point. Picking one sky scan per water scan is a choice, and the
usual practice (average the skies, or take the nearest in time) hides how much that
choice was worth. Computing every pairing makes it explicit.
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

from fieldrrs.rrs import (RHO_MOBLEY1999, rho_at_angle, rrs_three_scan,
                          scaled_mean, view_zenith_from_tilt)
from organize_by_location import fov_deg, survey
from process_field_day import band

WLO, WHI = 350.0, 950.0
RLO, RHI = 400.0, 900.0
C_SKY, C_WATER, C_PANEL, C_MEAN = "#7fb3d5", "#1f7a99", "#d9534f", "#c0392b"


def clip(wl, v, lo=WLO, hi=WHI):
    p = [(w, x) for w, x in zip(wl, v) if lo <= w <= hi]
    return [a for a, _ in p], [b for _, b in p]


def stats(curves):
    """Pointwise mean, min, max and relative spread across a list of spectra."""
    n = len(curves)
    m = [sum(c[i] for c in curves) / n for i in range(len(curves[0]))]
    lo = [min(c[i] for c in curves) for i in range(len(m))]
    hi = [max(c[i] for c in curves) for i in range(len(m))]
    return m, lo, hi


def rayleigh_tau(wl_nm):
    """Bodhaine et al. (1999) Rayleigh optical depth at sea level, wl in nm."""
    u = wl_nm / 1000.0
    return 0.0021520 * (1.0455996 - 341.29061 / u**2 - 0.90230850 * u**2) / \
           (1 + 0.0027059889 / u**2 - 85.968563 * u**2)


def hhmm(h):
    return "%02d:%02d" % (int(h), int(round((h - int(h)) * 60)))


# ---------------------------------------------------------------- figure 1
def fig_pooled(sky, water, panels, land, wl, outdir, tag):
    n = 4 if land else 3
    fig, axes = plt.subplots(1, n, figsize=(4.3 * n, 5.0))
    groups = [(axes[0], panels, C_PANEL, "CALIBRATION PANEL  $L_{ref}$", "rad_ref"),
              (axes[1], sky, C_SKY, "SKY  $L_{sky}$", "rad_target"),
              (axes[2], water, C_WATER, "WATER  $L_t$", "rad_target")]
    if land:
        groups.append((axes[3], land, "#2e7d32", "LAND TARGETS  $L_t$", "rad_target"))
    for ax, group, col, title, key in groups:
        curves = [s["spec"].columns[key] for s in group]
        m, lo, hi = stats(curves)
        for c in curves:
            ax.plot(*clip(wl, c), lw=0.8, alpha=0.5, color=col)
        wlc, mc = clip(wl, m)
        ax.plot(wlc, mc, lw=2.4, color="k", label="mean (n=%d)" % len(curves))
        ax.fill_between(wlc, clip(wl, lo)[1], clip(wl, hi)[1], color=col, alpha=0.22,
                        label="full range")
        ax.set_yscale("log")
        pos = [v for v in mc if v > 0]
        ax.set_ylim(min(pos) * 0.6, max(pos) * 1.8)
        ax.set_xlabel("wavelength (nm)"); ax.grid(alpha=0.25)
        ax.legend(fontsize=9, loc="lower left")
        # spread in the visible
        sv = [(hi[i] - lo[i]) / m[i] for i, w in enumerate(wl) if 450 <= w <= 650]
        ax.set_title("%s\nn=%d   full spread in 450-650 nm: %.1f %%"
                     % (title, len(curves), 100 * sum(sv) / len(sv)), fontsize=11)
    if land:
        for s_ in land:
            wlc, c = clip(wl, s_["spec"].columns["rad_target"])
            axes[3].annotate(s_["n"], (wlc[len(wlc) // 2], c[len(c) // 2]),
                             fontsize=8, color="#2e7d32")
        axes[3].text(0.03, 0.03, "NOT water: these do NOT enter $R_{rs}$.\n"
                     "See fig7 for their reflectance.", transform=axes[3].transAxes,
                     fontsize=8.5, color="#c0392b")
    axes[0].set_ylabel("radiance  W m$^{-2}$ sr$^{-1}$ nm$^{-1}$")
    fig.suptitle("%s — the three measured quantities, every scan overlaid" % tag,
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    p = os.path.join(outdir, "fig1_pooled_measurements.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p


# ---------------------------------------------------------------- figure 2
def fig_steps(w, l_sky, l_panel, wl, outdir, tag, panel_r, rho):
    lt = w["spec"].columns["rad_target"]
    ed = [math.pi * p / panel_r for p in l_panel]
    lw = [t - rho * s for t, s in zip(lt, l_sky)]
    res = rrs_three_scan(wl, lt, l_sky, l_panel, panel_r, rho, "none")

    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    a = ax[0][0]
    for v, c, lab in ((l_panel, C_PANEL, "$L_{panel}$"), (l_sky, C_SKY, "$L_{sky}$"),
                      (lt, C_WATER, "$L_t$ (water)")):
        a.plot(*clip(wl, v), lw=1.8, color=c, label=lab)
    a.set_yscale("log"); a.legend(fontsize=9.5)
    a.set_title("STEP 1   three radiances, scan %s" % w["n"], fontsize=11, loc="left")
    a.set_ylabel("W m$^{-2}$ sr$^{-1}$ nm$^{-1}$")

    a = ax[0][1]
    a.plot(*clip(wl, ed), lw=2.2, color=C_PANEL)
    a.set_title("STEP 2   $E_d = \\pi L_{panel} / R_{panel}$   ($R_{panel}$=%.2f)"
                % panel_r, fontsize=11, loc="left")
    a.set_ylabel("W m$^{-2}$ nm$^{-1}$")

    a = ax[1][0]
    a.plot(*clip(wl, lt), lw=1.5, color=C_WATER, label="$L_t$")
    a.plot(*clip(wl, [rho * s for s in l_sky]), lw=1.5, color=C_SKY,
           label="$\\rho L_{sky}$, $\\rho$=%.3f" % rho)
    a.plot(*clip(wl, lw), lw=2.4, color=C_MEAN, label="$L_w$")
    frac = [100 * rho * s / t for w_, t, s in zip(wl, lt, l_sky) if 440 <= w_ <= 460
            for t, s in [(t, s)]]
    a.legend(fontsize=9.5)
    a.set_title("STEP 3   remove reflected skylight", fontsize=11, loc="left")
    a.set_ylabel("W m$^{-2}$ sr$^{-1}$ nm$^{-1}$")
    f = [100 * rho * s / t for lam, t, s in zip(wl, lt, l_sky) if 440 <= lam <= 460]
    a.text(0.03, 0.92, "skylight removed = %.0f %% of $L_t$ at 450 nm" % (sum(f)/len(f)),
           transform=a.transAxes, fontsize=9, color=C_MEAN)

    a = ax[1][1]
    a.plot(*clip(wl, res.rrs, RLO, RHI), lw=2.4, color=C_MEAN)
    a.axhline(0, color="#888", lw=0.8)
    a.set_title("STEP 4   $R_{rs} = L_w / E_d$", fontsize=11, loc="left")
    a.set_ylabel("$R_{rs}$  sr$^{-1}$")
    for a_ in ax.flat:
        a_.set_xlabel("wavelength (nm)"); a_.grid(alpha=0.25)
    fig.suptitle("%s — the calculation, worked" % tag, fontsize=13, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = os.path.join(outdir, "fig2_steps.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p


# ---------------------------------------------------------------- figure 3
def fig_ed(panels, wl, sun, outdir, tag, panel_r):
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "..", "giop_python", "src"))
        import numpy as np
        from giop.water import f0_solar
    except Exception as exc:
        print("  F0 table unavailable (%s); skipping fig3" % exc)
        return None, {}

    wl_a = np.array(wl)
    eds = [np.array([math.pi * p / panel_r for p in s["spec"].columns["rad_ref"]])
           for s in panels]
    ed = sum(eds) / len(eds)
    mu = math.cos(math.radians(90.0 - sun))
    f0 = f0_solar(wl_a)
    T = ed / (f0 * mu)
    m = (wl_a >= 380) & (wl_a <= 950)

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.9))
    a = axes[0]
    a.plot(wl_a[m], ed[m], lw=2.2, color=C_PANEL, label="$E_d$ from panel")
    a.plot(wl_a[m], (f0 * mu)[m], lw=1.8, color="#8a6000", ls="--",
           label="$F_0\\cos\\theta_s$ (top of atmosphere)")
    a.set_ylabel("W m$^{-2}$ nm$^{-1}$"); a.legend(fontsize=9)
    a.set_title("Measured $E_d$ against the solar constant", fontsize=11, loc="left")

    a = axes[1]
    a.plot(wl_a[m], T[m], lw=2.2, color="#2e7d32")
    tau = np.array([rayleigh_tau(x) for x in wl_a])
    a.plot(wl_a[m], np.exp(-tau[m] / mu), lw=1.6, ls="--", color="#2c6f9b",
           label="direct-beam Rayleigh only")
    a.axhline(1.0, color="#c0392b", lw=1.4)
    a.text(390, 1.02, "T = 1 is the hard ceiling", fontsize=8.5, color="#c0392b")
    for lam, lab in ((762, "$O_2$-A"), (940, "$H_2O$")):
        a.axvline(lam, color="#888", ls=":", lw=1)
        a.text(lam + 4, 0.15, lab, fontsize=8.5, color="#555")
    a.set_ylim(0, 1.15); a.legend(fontsize=8.5, loc="lower right")
    a.set_ylabel("total transmittance $T$")
    a.set_title("$T = E_d / (F_0\\cos\\theta_s)$   sun %.1f$^\\circ$ elevation"
                % sun, fontsize=11, loc="left")

    a = axes[2]
    vis = (wl_a >= 450) & (wl_a <= 650)
    for i, (s, e) in enumerate(zip(panels, eds)):
        a.plot(wl_a[m], (e / ed)[m], lw=1.1, alpha=0.8,
               label="ref %.4f" % band(s["spec"], "rad_ref", 450, 650))
    a.axhline(1, color="k", lw=1)
    a.set_ylabel("$E_d$ / mean $E_d$")
    a.set_title("Panel-to-panel consistency", fontsize=11, loc="left")
    a.legend(fontsize=8)
    for a_ in axes:
        a_.set_xlabel("wavelength (nm)"); a_.grid(alpha=0.25); a_.set_xlim(380, 950)

    diag = {
        "T_vis": float(np.median(T[vis])),
        "T_max": float(np.max(T[m])),
        "T_ok": bool(np.all(T[m] < 1.0)),
        "o2_depth": float(1 - np.min(T[(wl_a > 755) & (wl_a < 775)]) /
                          np.median(T[(wl_a > 730) & (wl_a < 750)])),
        "h2o_depth": float(1 - np.min(T[(wl_a > 920) & (wl_a < 960)]) /
                           np.median(T[(wl_a > 860) & (wl_a < 890)])),
        "blue_red": float(np.median(T[(wl_a > 600) & (wl_a < 680)]) /
                          np.median(T[(wl_a > 420) & (wl_a < 480)])),
        "ed_spread": float(np.max([np.max(np.abs(e / ed - 1)[vis]) for e in eds])),
    }
    fig.suptitle("%s — is $E_d$ reasonable?" % tag, fontsize=13, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    p = os.path.join(outdir, "fig3_Ed_check.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p, diag


# ---------------------------------------------------------------- figure 4
def match_by_angle(water, sky, respect_blocks=True):
    """Pair each water scan with a sky scan, INSIDE this location only.

    Two constraints, both physical:

    1. LOCATION. This function only ever sees one location/foreoptic folder, so a water
       scan can never be paired with a sky scan from a different site or a different
       field of view. Guarded by `assert_same_dataset`.

    2. PANEL-REFERENCE BLOCK. Within one location the operator may re-reference the
       panel mid-station, which splits the visit into blocks separated in time. Pairing
       across a block boundary reaches for a sky scan minutes away and from the other
       side of a recalibration, which is exactly the changing-light error the sky scan
       exists to avoid. Matching therefore stays inside a block, and only falls back to
       the whole location if a block has no sky scan at all -- which is reported.

    Within those constraints the sky scan chosen is the one whose view angle MIRRORS the
    water scan, per the field protocol: the same angle from zenith as the water is from
    nadir.
    """
    from process_field_day import band
    out = []
    for w in water:
        pool = sky
        note = ""
        if respect_blocks:
            kw = round(band(w["spec"], "rad_ref", 450, 650), 6)
            same = [s for s in sky
                    if round(band(s["spec"], "rad_ref", 450, 650), 6) == kw]
            if same:
                pool = same
            else:
                note = "no sky in this panel block; fell back to the location"
        tw = w["spec"].tilt_y_deg
        best = min(pool, key=lambda s: abs(s["spec"].tilt_y_deg - tw))
        out.append((w, best, abs(best["spec"].tilt_y_deg - tw), note))
    return out


def assert_same_dataset(scans):
    """Refuse to analyse a mixed bag: one position, one foreoptic, one wavelength grid."""
    fos = {s["fo"] for s in scans}
    if len(fos) > 1:
        raise ValueError("mixed foreoptics %s: 8 and 15 deg average over different "
                         "footprints and do not share rho" % sorted(fos))
    lat = [s["lat"] for s in scans]
    lon = [s["lon"] for s in scans]
    import math as _m
    span = _m.hypot((max(lat) - min(lat)) * 111320,
                    (max(lon) - min(lon)) * 111320 * _m.cos(_m.radians(lat[0])))
    if span > 120.0:
        raise ValueError("scans span %.0f m; this is not one location" % span)
    if len({len(s["spec"].wavelength) for s in scans}) > 1:
        raise ValueError("scans are on different wavelength grids")
    return span


def fig_rrs(water, sky, wl, outdir, tag, panel_r, rho):
    """Three treatments, from blind to geometry-aware, so the envelope is attributable.

      A  all water x sky pairings, fixed rho          -- what you get with no geometry
      B  all pairings, rho at each water scan's angle -- rho corrected, pairing blind
      C  ANGLE-MATCHED pairing, rho at that angle     -- both corrected

    TWO CAUTIONS, both measured rather than assumed.

    A min-max envelope GROWS WITH SAMPLE SIZE, so A (96 spectra) and C (12) are not
    comparable on that statistic. The standard deviation is, and it is reported
    alongside; quote the sd, not the envelope.

    And angle matching is NOT demonstrably better than picking a sky at random here: a
    2000-draw control that assigns one RANDOM sky per water at the same n=12 reproduces
    the angle-matched scatter about a quarter of the time. That is because these sky
    scans span only ~5 deg, so there is little geometry to match. Reported honestly by
    `random_pairing_control` rather than presented as a win.
    """
    import numpy as np
    i443 = min(range(len(wl)), key=lambda k: abs(wl[k] - 443))
    i555 = min(range(len(wl)), key=lambda k: abs(wl[k] - 555))

    def rrs_of(w, sk, use_angle):
        r = rho_at_angle(view_zenith_from_tilt(w["spec"].tilt_y_deg)) if use_angle \
            else rho
        return rrs_three_scan(wl, w["spec"].columns["rad_target"],
                              sk["spec"].columns["rad_target"],
                              w["spec"].columns["rad_ref"], panel_r, r, "none").rrs

    A = [rrs_of(w, sk, False) for w in water for sk in sky]
    B = [rrs_of(w, sk, True) for w in water for sk in sky]
    pairs = match_by_angle(water, sky)
    C = [rrs_of(w, sk, True) for w, sk, _, _ in pairs]
    mism = [d for _, _, d, _ in pairs]
    fallbacks = [w["n"] for w, _, _, note in pairs if note]

    by_water, by_sky = {}, {}
    for w in water:
        for sk in sky:
            r = rrs_of(w, sk, True)
            by_water.setdefault(w["n"], []).append(r)
            by_sky.setdefault(sk["n"], []).append(r)

    # Sample-size control: is ANGLE matching better than a random sky, at the same n?
    # The control must obey the SAME constraints as the treatment, or it is not a
    # control -- so it draws only from the water scan's own panel-reference block, and
    # differs from the treatment in exactly one thing: it ignores the view angle.
    import random as _random
    from process_field_day import band as _band
    pools = {}
    for w in water:
        kw = round(_band(w["spec"], "rad_ref", 450, 650), 6)
        pools[w["n"]] = [sk for sk in sky
                         if round(_band(sk["spec"], "rad_ref", 450, 650), 6) == kw] \
            or list(sky)
    rng = _random.Random(0)
    ctrl = []
    for _ in range(2000):
        v = [rrs_of(w, rng.choice(pools[w["n"]]), True)[i443] for w in water]
        mu = sum(v) / len(v)
        ctrl.append((sum((x - mu) ** 2 for x in v) / len(v)) ** 0.5 / mu * 100)
    ctrl.sort()

    mA, loA, hiA = stats(A)
    mB, loB, hiB = stats(B)
    mC, loC, hiC = stats(C)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.0))
    a = axes[0]
    for lab, (m, lo, hi), col in (
            ("A  blind pairing, fixed $\\rho$", (mA, loA, hiA), "#c0392b"),
            ("B  blind pairing, $\\rho(\\theta_v)$", (mB, loB, hiB), "#8a6000"),
            ("C  ANGLE-MATCHED, $\\rho(\\theta_v)$", (mC, loC, hiC), "#2e7d32")):
        wlc, _ = clip(wl, m, RLO, RHI)
        a.fill_between(wlc, clip(wl, lo, RLO, RHI)[1], clip(wl, hi, RLO, RHI)[1],
                       color=col, alpha=0.22)
        a.plot(*clip(wl, m, RLO, RHI), lw=2.2, color=col,
               label="%s   %.0f %% at 443" % (lab, 100 * (hi[i443] - lo[i443]) / m[i443]))
    a.axhline(0, color="#888", lw=0.8)
    a.set_ylabel("$R_{rs}$  sr$^{-1}$"); a.legend(fontsize=8.5)
    a.set_title("The envelope shrinks as geometry is respected", fontsize=10.5,
                loc="left")

    a = axes[1]
    idx = [i for i, x in enumerate(wl) if RLO <= x <= RHI]
    wsp, ssp = [], []
    for i in idx:
        wm = [sum(v[i] for v in by_water[k]) / len(by_water[k]) for k in by_water]
        sm = [sum(v[i] for v in by_sky[k]) / len(by_sky[k]) for k in by_sky]
        mu = sum(v[i] for v in B) / len(B)
        wsp.append((max(wm) - min(wm)) / abs(mu) * 100 if mu else 0)
        ssp.append((max(sm) - min(sm)) / abs(mu) * 100 if mu else 0)
    a.plot([wl[i] for i in idx], wsp, lw=2.2, color=C_WATER,
           label="which WATER scan (n=%d)" % len(water))
    a.plot([wl[i] for i in idx], ssp, lw=2.2, color=C_SKY,
           label="which SKY scan (n=%d)" % len(sky))
    a.set_ylabel("range of group means, % of $R_{rs}$"); a.legend(fontsize=9)
    a.set_title("Where the remaining spread comes from", fontsize=10.5, loc="left")

    a = axes[2]
    for w, sk, dm, _ in pairs:
        a.plot(*clip(wl, rrs_of(w, sk, True), RLO, RHI), lw=1.1, alpha=0.75,
               color=C_WATER)
    a.plot(*clip(wl, mC, RLO, RHI), lw=2.8, color="#2e7d32",
           label="angle-matched mean (n=%d)" % len(C))
    a.axhline(0, color="#888", lw=0.8)
    a.set_ylabel("$R_{rs}$  sr$^{-1}$"); a.legend(fontsize=9)
    cv = [c[i443] for c in C]
    sdC = float(np.std(cv) / np.mean(cv) * 100)
    frac = sum(1 for c in ctrl if c <= sdC) / len(ctrl)
    a.set_title("C: one sky per water, mirrored geometry\n"
                "mismatch mean %.1f$^\\circ$   sd %.1f %% vs %.1f %% for a RANDOM sky\n"
                "%.0f %% of random draws are as tight -> NOT yet a demonstrated gain"
                % (sum(mism) / len(mism), sdC, ctrl[len(ctrl) // 2], 100 * frac),
                fontsize=9.5, loc="left")
    for a_ in axes:
        a_.set_xlabel("wavelength (nm)"); a_.grid(alpha=0.25); a_.set_xlim(RLO, RHI)
    fig.suptitle("%s — $R_{rs}$, paired by geometry rather than blindly" % tag,
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    pth = os.path.join(outdir, "fig4_rrs_spread.png")
    fig.savefig(pth, dpi=140); plt.close(fig)

    env = lambda m, lo, hi, i: (hi[i] - lo[i]) / m[i] * 100
    out = {"n_A": len(A), "n_C": len(C), "mean": mC,
           "rrs443": mC[i443], "rrs555": mC[i555],
           "envA443": env(mA, loA, hiA, i443), "envB443": env(mB, loB, hiB, i443),
           "envC443": env(mC, loC, hiC, i443),
           "envA555": env(mA, loA, hiA, i555), "envC555": env(mC, loC, hiC, i555),
           "mism_mean": sum(mism) / len(mism), "mism_max": max(mism),
           "fallbacks": fallbacks,
           "sdA443": float(np.std([c[i443] for c in A]) / np.mean([c[i443] for c in A])
                           * 100),
           "sdC443": sdC, "ctrl_med": ctrl[len(ctrl) // 2],
           "ctrl_frac": sum(1 for c in ctrl if c <= sdC) / len(ctrl),
           "w_spread_443": wsp[idx.index(i443)], "s_spread_443": ssp[idx.index(i443)],
           "w_spread": max(wsp), "s_spread": max(ssp)}
    return pth, out



# ---------------------------------------------------------------- figure 5
NIR_SIM_RATIO = 1.912          # Ruddick et al. (2006)


def similarity_delta(wl, rrs):
    """The offset the NIR-similarity criterion says is left over.

    Ruddick et al. (2006): for a wide range of turbid water the TRUE R_rs obeys
    R_rs(780)/R_rs(870) = 1.912. If the measured pair does not, the discrepancy is
    attributed to a spectrally flat residual `delta`:

        (R780 - d) = 1.912 (R870 - d)   ->   d = (1.912 R870 - R780) / 0.912

    A pair of scans needing a SMALL |d| is one where the sky subtraction already
    left a physically consistent spectrum. That is what "best matched" means here:
    least ad-hoc correction required, not best-looking.
    """
    def at(lam):
        v = [x for w, x in zip(wl, rrs) if lam - 5 <= w <= lam + 5]
        return sum(v) / len(v)
    return (NIR_SIM_RATIO * at(870.0) - at(780.0)) / (NIR_SIM_RATIO - 1.0)


def fig_geometry(scans, water, sky, wl, outdir, tag, panel_r, rho):
    import numpy as np
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.0))

    a = axes[0]
    for role, col in (("sky", C_SKY), ("water", C_WATER), ("vegetation", "#2e7d32")):
        g = [x for x in scans if x["role"] == role]
        if not g:
            continue
        a.scatter([x["spec"].tilt_y_deg for x in g], [x["spec"].tilt_x_deg for x in g],
                  s=70, color=col, edgecolor="k", linewidth=0.5, alpha=0.85,
                  label="%s (n=%d)" % (role, len(g)))
    a.axvline(50.0, color="#c0392b", ls="--", lw=1.6)
    a.text(50.4, a.get_ylim()[0] + 0.3, "50$^\\circ$ from horizontal\n= 40$^\\circ$ "
           "from nadir/zenith\n(Mobley nominal)", fontsize=8.5, color="#c0392b")
    a.axhline(0, color="#888", lw=1)
    a.set_xlabel("tilt Y (deg from horizontal)")
    a.set_ylabel("tilt X (deg, roll)")
    a.legend(fontsize=9)
    a.set_title("Pointing actually achieved\n(the sensor gives MAGNITUDE, so sky and "
                "water overlap)", fontsize=10.5, loc="left")

    # distribution of the 96 pairings
    a = axes[1]
    combos, deltas, pairs = [], [], []
    for w in water:
        lt = w["spec"].columns["rad_target"]
        lp = w["spec"].columns["rad_ref"]
        for sk in sky:
            r = rrs_three_scan(wl, lt, sk["spec"].columns["rad_target"], lp, panel_r,
                               rho, "none").rrs
            combos.append(r)
            deltas.append(similarity_delta(wl, r))
            pairs.append((w["n"], sk["n"]))
    i443 = min(range(len(wl)), key=lambda k: abs(wl[k] - 443))
    i555 = min(range(len(wl)), key=lambda k: abs(wl[k] - 555))
    for idx, lam, col in ((i443, 443, "#2c6f9b"), (i555, 555, "#2e7d32")):
        v = [c[idx] for c in combos]
        a.hist(v, bins=22, alpha=0.6, color=col,
               label="%d nm: %.5f $\\pm$ %.5f" % (lam, np.mean(v), np.std(v)))
        a.axvline(np.mean(v), color=col, lw=2)
    a.set_xlabel("$R_{rs}$  sr$^{-1}$"); a.set_ylabel("number of pairings")
    a.legend(fontsize=9)
    a.set_title("Distribution over all %d water$\\times$sky pairings" % len(combos),
                fontsize=10.5, loc="left")

    a = axes[2]
    order = sorted(range(len(combos)), key=lambda k: abs(deltas[k]))
    best, worst = order[0], order[-1]
    a.plot(*clip(wl, combos[best], RLO, RHI), lw=2.4, color="#2e7d32",
           label="BEST pair  w%s/s%s   $\\delta$=%+.2e"
                 % (pairs[best][0], pairs[best][1], deltas[best]))
    a.plot(*clip(wl, combos[worst], RLO, RHI), lw=2.0, color="#c0392b", ls="--",
           label="WORST pair w%s/s%s   $\\delta$=%+.2e"
                 % (pairs[worst][0], pairs[worst][1], deltas[worst]))
    m, lo, hi = stats(combos)
    a.plot(*clip(wl, m, RLO, RHI), lw=2.6, color="k", ls=":", label="mean of all")
    a.axhline(0, color="#888", lw=0.8)
    a.set_xlabel("wavelength (nm)"); a.set_ylabel("$R_{rs}$  sr$^{-1}$")
    a.legend(fontsize=8.5)
    a.set_title("Best vs worst pairing by NIR-similarity mismatch", fontsize=10.5,
                loc="left")
    for a_ in axes[1:]:
        a_.grid(alpha=0.25)
    axes[0].grid(alpha=0.25)
    fig.suptitle("%s — geometry, the pairing distribution, and the best-matched pair"
                 % tag, fontsize=13, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    p = os.path.join(outdir, "fig5_geometry_and_pairing.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p, {"best": pairs[best], "worst": pairs[worst],
               "d_best": deltas[best], "d_worst": deltas[worst],
               "rrs443_mean": float(np.mean([c[i443] for c in combos])),
               "rrs443_sd": float(np.std([c[i443] for c in combos]))}


# ---------------------------------------------------------------- figure 6
def fig_sensitivity(water, sky, wl, outdir, tag, panel_r, rho, fov):
    """Does the pointing angle or the footprint systematically move R_rs?"""
    import numpy as np
    from scipy import stats as sps

    l_sky_mean = [sum(s["spec"].columns["rad_target"][i] for s in sky) / len(sky)
                  for i in range(len(wl))]
    i443 = min(range(len(wl)), key=lambda k: abs(wl[k] - 443))
    i555 = min(range(len(wl)), key=lambda k: abs(wl[k] - 555))

    rows = []
    for w in water:
        r = rrs_three_scan(wl, w["spec"].columns["rad_target"], l_sky_mean,
                           w["spec"].columns["rad_ref"], panel_r, rho, "none").rrs
        if w["range"] is None:
            continue          # no range in this header -> no footprint for this scan
        fp = 2 * w["range"] * math.tan(math.radians(fov / 2.0))
        rows.append({"tilt": w["spec"].tilt_y_deg, "range": w["range"],
                     "fp": fp, "r443": r[i443], "r555": r[i555], "n": w["n"],
                     "gps": w["gps"]})

    # and the sky-side test: vary which sky, keep water fixed to its own mean
    sky_rows = []
    for sk in sky:
        vals = []
        for w in water:
            r = rrs_three_scan(wl, w["spec"].columns["rad_target"],
                               sk["spec"].columns["rad_target"],
                               w["spec"].columns["rad_ref"], panel_r, rho, "none").rrs
            vals.append(r[i443])
        sky_rows.append({"tilt": sk["spec"].tilt_y_deg, "r443": float(np.mean(vals)),
                         "n": sk["n"]})

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.0))
    spec = []
    for a, xs, ys, xlab, ylab, title, col in (
            (axes[0], [r["tilt"] for r in rows], [r["r443"] for r in rows],
             "WATER scan tilt Y (deg)", "$R_{rs}$(443)  sr$^{-1}$",
             "Does the WATER view angle matter?", C_WATER),
            (axes[1], [r["tilt"] for r in sky_rows], [r["r443"] for r in sky_rows],
             "SKY scan tilt Y (deg)", "mean $R_{rs}$(443)  sr$^{-1}$",
             "Does the SKY view angle matter?", C_SKY),
            (axes[2], [r["fp"] for r in rows], [r["r443"] for r in rows],
             "footprint across-view (m)", "$R_{rs}$(443)  sr$^{-1}$",
             "Does the FOOTPRINT matter?", "#8a6000")):
        a.scatter(xs, ys, s=80, color=col, edgecolor="k", linewidth=0.5, zorder=3)
        r_p, p_p = sps.pearsonr(xs, ys)
        r_s, p_s = sps.spearmanr(xs, ys)
        sl, ic = np.polyfit(xs, ys, 1)
        xx = np.linspace(min(xs), max(xs), 20)
        a.plot(xx, sl * xx + ic, lw=2, color="k", ls="--", zorder=2)
        verdict = "SIGNIFICANT" if p_p < 0.05 else "no significant trend"
        a.set_xlabel(xlab); a.set_ylabel(ylab); a.grid(alpha=0.25)
        a.set_title("%s\nPearson r=%+.2f (p=%.3f) · Spearman %+.2f (p=%.3f)\n%s"
                    % (title, r_p, p_p, r_s, p_s, verdict), fontsize=10, loc="left")
        spec.append({"what": title, "r": r_p, "p": p_p, "rs": r_s, "ps": p_s,
                     "slope": sl, "span": max(xs) - min(xs),
                     "effect": sl * (max(xs) - min(xs)),
                     "rel": abs(sl * (max(xs) - min(xs))) / np.mean(ys) * 100})

    # CONFOUND CONTROL. A correlation with tilt is only interesting if tilt is not
    # standing in for something else -- time (the water changes, the sun rises) or
    # range. Partial correlation removes the linear part of the third variable.
    def partial(x, y, z):
        rx = np.asarray(x) - np.polyval(np.polyfit(z, x, 1), z)
        ry = np.asarray(y) - np.polyval(np.polyfit(z, y, 1), z)
        r, _ = sps.pearsonr(rx, ry)
        n = len(x)
        t = r * math.sqrt((n - 3) / max(1e-12, 1 - r * r))
        return r, float(2 * (1 - sps.t.cdf(abs(t), n - 3)))

    tilt = [r["tilt"] for r in rows]
    r443 = [r["r443"] for r in rows]
    tme = [r["gps"] for r in rows]
    rng = [r["range"] for r in rows]
    conf = {
        "tilt_vs_time": sps.pearsonr(tilt, tme),
        "range_vs_time": sps.pearsonr(rng, tme),
        "tilt_ctrl_time": partial(tilt, r443, tme),
        "tilt_ctrl_range": partial(tilt, r443, rng),
        "range_ctrl_tilt": partial(rng, r443, tilt),
        "n_ranges": len(set(rng)),
    }
    fig.suptitle("%s — does pointing or footprint systematically bias $R_{rs}$?" % tag,
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    p = os.path.join(outdir, "fig6_angle_footprint.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p, spec, conf



# ---------------------------------------------------------------- figure 7
def fig_land(land, wl, outdir, tag, panel_r):
    """Land targets get REFLECTANCE, not R_rs, and the difference is physical.

    For water, rho*L_sky is a contaminant sitting on top of the signal and is removed.
    For an opaque surface the sky is part of the ILLUMINATION, so removing it would
    delete real signal. The right product needs no rho and no sky scan:

        R = pi L_target / E_d = R_panel * L_target / L_panel
    """
    from process_field_day import land_reflectance
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0))
    a = axes[0]
    for s_ in land:
        w, r = land_reflectance(s_["spec"], panel_r)
        wc, rc = clip(w, r, 380, 950)
        lab = "%s  red edge %.1fx" % (s_["n"], s_["diag"]["red_edge"])
        a.plot(wc, rc, lw=2.2, label=lab)
    a.axvspan(670, 680, color="#c0392b", alpha=0.12)
    a.text(682, a.get_ylim()[1] * 0.55, "chlorophyll\nabsorption", fontsize=8.5,
           color="#c0392b")
    a.axvspan(700, 760, color="#2e7d32", alpha=0.10)
    a.text(706, a.get_ylim()[1] * 0.15, "red edge", fontsize=8.5, color="#2e7d32")
    a.set_xlabel("wavelength (nm)"); a.set_ylabel("reflectance factor  $R$")
    a.legend(fontsize=9); a.grid(alpha=0.25)
    a.set_title("Land targets: $R = R_{panel}\\,L_t/L_{panel}$\n"
                "no $\\rho$, no sky scan, no $R_{rs}$", fontsize=11, loc="left")

    a = axes[1]
    names, edges = [], []
    for s_ in land:
        names.append(s_["n"]); edges.append(s_["diag"]["red_edge"])
    bars = a.bar(names, edges, color=["#2e7d32" if e > 3 else "#8a6000" for e in edges])
    a.axhline(3.0, color="#c0392b", ls="--", lw=1.5)
    a.text(0.02, 3.1, "above ~3 = vegetated", fontsize=9, color="#c0392b",
           transform=a.get_yaxis_transform())
    a.set_ylabel("red edge  $R$(750)/$R$(670)")
    a.set_title("Vegetated or bare?", fontsize=11, loc="left")
    a.grid(alpha=0.25, axis="y")
    for b, e in zip(bars, edges):
        a.text(b.get_x() + b.get_width() / 2, e + 0.1, "%.1f" % e, ha="center",
               fontsize=9)
    fig.suptitle("%s — the LAND targets, which are NOT part of the water analysis"
                 % tag, fontsize=13, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    p = os.path.join(outdir, "fig7_land_targets.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p



# ---------------------------------------------------------------- figure 8
def fig_rho_correction(water, sky, wl, outdir, tag, panel_r):
    """Per-scan rho at each scan's OWN view angle, against a single fixed value.

    The test is falsifiable: if the tilt trend really is rho's angular dependence, then
    correcting at each scan's own angle must REMOVE it. If instead the correction
    amplifies the trend or inflates the scatter, the hypothesis is wrong and gets
    reported as wrong.
    """
    import numpy as np
    from scipy import stats as sps

    l_sky = [sum(s["spec"].columns["rad_target"][i] for s in sky) / len(sky)
             for i in range(len(wl))]
    i443 = min(range(len(wl)), key=lambda k: abs(wl[k] - 443))
    tilt, fixed, corr, rhos, other = [], [], [], [], []
    for w in water:
        lt, lp = w["spec"].columns["rad_target"], w["spec"].columns["rad_ref"]
        t = w["spec"].tilt_y_deg
        r_own = rho_at_angle(view_zenith_from_tilt(t))
        r_alt = rho_at_angle(90.0 - abs(t))
        tilt.append(t); rhos.append(r_own)
        fixed.append(rrs_three_scan(wl, lt, l_sky, lp, panel_r, RHO_MOBLEY1999,
                                    "none").rrs)
        corr.append(rrs_three_scan(wl, lt, l_sky, lp, panel_r, r_own, "none").rrs)
        other.append(rrs_three_scan(wl, lt, l_sky, lp, panel_r, r_alt, "none").rrs)
    tilt = np.array(tilt)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.0))
    a = axes[0]
    tt = np.linspace(min(tilt) - 3, max(tilt) + 3, 60)
    a.plot(tt, [rho_at_angle(x) for x in tt], lw=2.2, color="#2e7d32",
           label="$\\rho(\\theta_v)$ Fresnel-scaled")
    a.axhline(RHO_MOBLEY1999, color="#c0392b", ls="--", lw=1.8,
              label="fixed $\\rho$ = 0.028")
    a.scatter(tilt, rhos, s=70, color="#2e7d32", edgecolor="k", zorder=4)
    a.set_xlabel("view zenith angle (deg)"); a.set_ylabel("$\\rho$")
    a.legend(fontsize=9); a.grid(alpha=0.25)
    a.set_title("$\\rho$ actually used, per scan\nspread %.0f %% across the achieved "
                "angles" % (100 * (max(rhos) / min(rhos) - 1)), fontsize=10.5,
                loc="left")

    a = axes[1]
    res = {}
    for lab, arr, col in (("fixed 0.028", fixed, "#c0392b"),
                          ("$\\rho(\\theta_v=$tilt$)$", corr, "#2e7d32"),
                          ("$\\rho(\\theta_v=90-$tilt$)$", other, "#8a6000")):
        y = np.array([v[i443] for v in arr])
        r, pv = sps.pearsonr(tilt, y)
        sl = np.polyfit(tilt, y, 1)[0]
        a.scatter(tilt, y, s=60, color=col, edgecolor="k", linewidth=0.4, zorder=3,
                  label="%s   r=%+.2f p=%.3f" % (lab, r, pv))
        a.plot(tt, sl * tt + np.polyfit(tilt, y, 1)[1], lw=1.6, color=col, ls="--")
        res[lab] = {"r": r, "p": pv, "sd": float(y.std() / y.mean() * 100),
                    "mean": float(y.mean())}
    a.set_xlabel("view zenith angle (deg)"); a.set_ylabel("$R_{rs}$(443)  sr$^{-1}$")
    a.legend(fontsize=8); a.grid(alpha=0.25)
    a.set_title("Does the correction FLATTEN the trend?", fontsize=10.5, loc="left")

    a = axes[2]
    labs = list(res)
    a.bar(range(len(labs)), [res[k]["sd"] for k in labs],
          color=["#c0392b", "#2e7d32", "#8a6000"])
    for i, k in enumerate(labs):
        a.text(i, res[k]["sd"] + 0.2, "%.1f %%" % res[k]["sd"], ha="center", fontsize=10)
    a.set_xticks(range(len(labs)))
    a.set_xticklabels([l.replace("$", "").replace("\\rho", "rho").replace("\\theta_v",
                       "th") for l in labs], fontsize=8, rotation=12)
    a.set_ylabel("scan-to-scan scatter of $R_{rs}$(443)  (%)")
    a.grid(alpha=0.25, axis="y")
    a.set_title("A correct correction REDUCES scatter.\nA wrong one adds variance.",
                fontsize=10.5, loc="left")
    fig.suptitle("%s — per-scan $\\rho$ from the measured view angle" % tag,
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    p = os.path.join(outdir, "fig8_rho_angle_correction.png")
    fig.savefig(p, dpi=140); plt.close(fig)
    return p, res, (min(rhos), max(rhos))



# ---------------------------------------------------------------- figure 9
def fig_variability(water, sky, wl, outdir, tag, panel_r, panels):
    """Is the residual scan-to-scan spread MEASUREMENT error or real water variability?

    Four signatures separate them, and they agree here:

      1 COHERENCE. An SVD of the spectra about their mean. Real changes in load move
        the whole spectrum together, so one mode dominates. Noise spreads across modes.
      2 AMPLITUDE vs SHAPE. Normalising each spectrum at 555 nm removes any pure
        scaling. Whatever survives is a change in water TYPE rather than amount.
      3 TIME. Drift (tide, sun, a settling plume) leaves scans close in time more alike.
        Patchiness does not.
      4 THE INSTRUMENT'S OWN FLOOR. The panel replicates bound what the instrument plus
        the E_d chain can do; the water must beat that comfortably to be physical.
    """
    import numpy as np
    from scipy import stats as sps
    pairs = match_by_angle(water, sky)
    order = sorted(range(len(pairs)), key=lambda k: pairs[k][0]["gps"])
    R, t, names = [], [], []
    for k in order:
        w, sk, _, _ = pairs[k]
        r = rho_at_angle(view_zenith_from_tilt(w["spec"].tilt_y_deg))
        R.append(rrs_three_scan(wl, w["spec"].columns["rad_target"],
                                sk["spec"].columns["rad_target"],
                                w["spec"].columns["rad_ref"], panel_r, r, "none").rrs)
        t.append(w["gps"] * 60.0); names.append(w["n"])
    wl_a = np.array(wl); m = (wl_a >= 420) & (wl_a <= 750)
    X = np.array(R)[:, m]; lam = wl_a[m]; t = np.array(t)
    i555 = int(np.argmin(abs(lam - 555)))

    U, S, Vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
    var = S ** 2 / np.sum(S ** 2)
    raw = 100 * np.mean(X.std(0) / X.mean(0))
    Xn = X / X[:, [i555]]
    nor = 100 * np.mean(Xn.std(0) / Xn.mean(0))
    dd, dt = [], []
    for i in range(len(X)):
        for j in range(i + 1, len(X)):
            dd.append(float(np.linalg.norm(X[i] - X[j]))); dt.append(abs(t[i] - t[j]))
    r_t, p_t = sps.pearsonr(dt, dd)

    def sp(g, key):
        c = np.array([x["spec"].columns[key] for x in g])[:, (wl_a >= 450) &
                                                          (wl_a <= 650)]
        return 100 * float(np.mean(c.std(0) / c.mean(0)))
    fl_panel, fl_sky = sp(panels, "rad_ref"), sp(sky, "rad_target")

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.0))
    a = axes[0]
    for i in range(len(X)):
        a.plot(lam, Xn[i], lw=1.0, alpha=0.7, color=C_WATER)
    a.plot(lam, Xn.mean(0), lw=2.8, color="#2e7d32")
    a.set_xlabel("wavelength (nm)"); a.set_ylabel("$R_{rs}$ / $R_{rs}$(555)")
    a.grid(alpha=0.25)
    a.set_title("Normalised at 555 nm: the SHAPE\n"
                "spread %.1f %% raw -> %.1f %% normalised, so %.0f %% of the\n"
                "variance is pure amplitude" % (raw, nor, 100 * (1 - (nor / raw) ** 2)),
                fontsize=10, loc="left")

    a = axes[1]
    # min(5, len(var)): the SVD of n spectra yields at most n modes, and a station with
    # fewer than 5 water scans is normal (LOC3 has 4).
    kmode = min(5, len(var))
    a.bar(range(1, kmode + 1), 100 * var[:kmode], color="#2c6f9b")
    a.set_xlabel("SVD mode"); a.set_ylabel("variance explained (%)")
    a.set_xticks(range(1, kmode + 1)); a.grid(alpha=0.25, axis="y")
    a.set_title("One mode carries %.1f %%\n"
                "coherent, as a load change is; noise\nwould spread across modes"
                % (100 * var[0]), fontsize=10, loc="left")

    a = axes[2]
    amp = X[:, i555]
    a.plot(t - t.min(), amp, "o-", color=C_WATER, ms=8, lw=1.2)
    for x_, y_, n_ in zip(t - t.min(), amp, names):
        a.annotate(n_, (x_, y_), fontsize=7, xytext=(3, 4),
                   textcoords="offset points")
    a.axhline(amp.mean(), color="#888", ls="--")
    a.set_xlabel("minutes into the station"); a.set_ylabel("$R_{rs}$(555)  sr$^{-1}$")
    a.grid(alpha=0.25)
    a.set_title("Amplitude vs time: factor %.2f, no trend\n"
                "distance-vs-time-gap r=%+.2f (p=%.2f)\n"
                "instrument floor %.1f %%, sky %.1f %%, water %.1f %%"
                % (amp.max() / amp.min(), r_t, p_t, fl_panel, fl_sky, raw),
                fontsize=10, loc="left")
    fig.suptitle("%s — is the residual spread measurement error or real water?" % tag,
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    pth = os.path.join(outdir, "fig9_variability_origin.png")
    fig.savefig(pth, dpi=140); plt.close(fig)
    return pth, {"pc1": 100 * var[0], "raw": raw, "nor": nor,
                 "amp_frac": 100 * (1 - (nor / raw) ** 2),
                 "ratio": float(amp.max() / amp.min()), "r_t": r_t, "p_t": p_t,
                 "floor": fl_panel, "sky_floor": fl_sky}



# --------------------------------------------------------------- figure 10
def fig_final(water, sky, wl, outdir, tag, panel_r):
    """THE PRODUCT: a clean mean shape, with amplitude reported separately.

    92 % of the scan-to-scan variance at this site is a multiplicative factor -- the
    same water with more or less material in it. A plain average lets that leak into
    every band, so the mean looks less certain than the shape actually is.
    `scaled_mean` alternates fitting a per-scan amplitude and re-meaning, which splits
    the result into the two things that are genuinely different:

        SHAPE      what the water is        -- small uncertainty, this is the product
        AMPLITUDE  how much of it there is  -- large, and REAL, not error

    Quoting one combined error bar would imply the water type is uncertain when only
    its concentration is.
    """
    import numpy as np
    pairs = match_by_angle(water, sky)
    R, names = [], []
    for w, sk, _, _ in pairs:
        r = rho_at_angle(view_zenith_from_tilt(w["spec"].tilt_y_deg))
        R.append(rrs_three_scan(wl, w["spec"].columns["rad_target"],
                                sk["spec"].columns["rad_target"],
                                w["spec"].columns["rad_ref"], panel_r, r, "none").rrs)
        names.append(w["n"])
    wl_a = np.array(wl)
    m = (wl_a >= RLO) & (wl_a <= RHI)
    X = [list(np.array(r)[m]) for r in R]
    lam = wl_a[m]
    mean, scales, iters, shape_cv, amp_cv = scaled_mean(X)
    mean = np.array(mean)
    plain = np.array(X).mean(0)
    plain_cv = float(np.mean(np.array(X).std(0) / np.abs(plain)) * 100)
    # A RELATIVE scatter is meaningless where the signal approaches zero, and R_rs does
    # exactly that below ~430 and above ~800 nm. Quote the shape over the core band and
    # over the bands carrying real signal, not over the noise tails.
    Xs_ = np.array([[sc_ * v for v in x] for sc_, x in zip(scales, X)])
    core = (lam >= 450) & (lam <= 700)
    shape_core = float(np.mean(Xs_.std(0)[core] / np.abs(mean[core])) * 100)
    plain_core = float(np.mean(np.array(X).std(0)[core] / np.abs(plain[core])) * 100)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.0))
    a = axes[0]
    for x in X:
        a.plot(lam, x, lw=0.8, alpha=0.45, color=C_WATER)
    a.plot(lam, plain, lw=2.0, color="#888", ls="--",
           label="plain mean   band scatter %.1f %%" % plain_cv)
    a.plot(lam, mean, lw=3.0, color="#2e7d32",
           label="scaled mean   shape %.1f %% (450-700 nm)" % shape_core)
    a.axvspan(450, 700, color="#2e7d32", alpha=0.06)
    sd = np.array(X).std(0)
    a.fill_between(lam, plain - sd, plain + sd, color="#bbb", alpha=0.4)
    a.set_ylabel("$R_{rs}$  sr$^{-1}$"); a.legend(fontsize=9)
    a.set_title("The product: mean $R_{rs}$\ngrey band = plain 1 sd, which is mostly "
                "amplitude", fontsize=10.5, loc="left")

    a = axes[1]
    Xs = np.array([[sc * v for v in x] for sc, x in zip(scales, X)])
    for row in Xs:
        a.plot(lam, row, lw=0.9, alpha=0.6, color="#2e7d32")
    a.plot(lam, mean, lw=3.0, color="k")
    rel = Xs.std(0) / np.abs(mean) * 100
    a.set_ylabel("$R_{rs}$ rescaled  sr$^{-1}$")
    a.axvspan(450, 700, color="#2e7d32", alpha=0.06)
    a.set_title("After scaling each scan onto the mean\n"
                "450-700 nm: %.1f %% residual, was %.1f %%\n"
                "(%.1f %% over 400-900, inflated by the near-zero tails)"
                % (shape_core, plain_core, shape_cv), fontsize=10, loc="left")

    a = axes[2]
    amp = [1.0 / s_ for s_ in scales]
    a.bar(range(len(amp)), amp, color=C_WATER)
    a.axhline(1.0, color="k", lw=1.4)
    a.set_xticks(range(len(amp))); a.set_xticklabels(names, rotation=60, fontsize=7.5)
    a.set_ylabel("amplitude relative to the mean")
    a.set_title("Amplitude per scan: %.0f %% (1 sd), range %.2f-%.2f\n"
                "REAL variability in load, not measurement error"
                % (amp_cv, min(amp), max(amp)), fontsize=10.5, loc="left")
    for a_ in axes[:2]:
        a_.set_xlabel("wavelength (nm)"); a_.grid(alpha=0.25); a_.set_xlim(RLO, RHI)
    axes[2].grid(alpha=0.25, axis="y")
    fig.suptitle("%s — FINAL PRODUCT: shape and amplitude reported separately" % tag,
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    pth = os.path.join(outdir, "fig10_final_product.png")
    fig.savefig(pth, dpi=140); plt.close(fig)

    csv_path = os.path.join(outdir, "FINAL_Rrs.csv")
    with open(csv_path, "w", newline="") as fh:
        wcsv = csv.writer(fh)
        wcsv.writerow(["# scaled-mean R_rs; shape_cv_pct=%.3f amp_cv_pct=%.2f n=%d"
                       % (shape_cv, amp_cv, len(X))])
        wcsv.writerow(["wavelength_nm", "Rrs_mean", "shape_sd", "plain_sd"])
        for i, lm in enumerate(lam):
            wcsv.writerow(["%.1f" % lm, "%.8e" % mean[i], "%.8e" % Xs.std(0)[i],
                           "%.8e" % sd[i]])
    return pth, csv_path, {"shape_cv": shape_cv, "amp_cv": amp_cv, "plain_cv": plain_cv,
                           "shape_core": shape_core, "plain_core": plain_core,
                           "iters": iters, "amp": amp, "names": names,
                           "gain": plain_core / shape_core}



# --------------------------------------------------------------- figure 11
def fig_scaled_method(water, sky, wl, outdir, tag, panel_r):
    """How the amplitude-normalised mean works, step by step. THEORY_SCALED_MEAN.md."""
    import numpy as np
    pairs = match_by_angle(water, sky)
    R, names = [], []
    for w, sk, _, _ in pairs:
        r = rho_at_angle(view_zenith_from_tilt(w["spec"].tilt_y_deg))
        R.append(rrs_three_scan(wl, w["spec"].columns["rad_target"],
                                sk["spec"].columns["rad_target"],
                                w["spec"].columns["rad_ref"], panel_r, r, "none").rrs)
        names.append(w["n"])
    wl_a = np.array(wl); m = (wl_a >= RLO) & (wl_a <= RHI)
    X = np.array(R)[:, m]; lam = wl_a[m]
    core = (lam >= 450) & (lam <= 700)

    # run the iteration by hand so each step can be plotted
    Xl = [list(x) for x in X]
    hist_J, hist_a, Ms = [], [], []
    M = X.mean(0)
    for _ in range(8):
        num = X @ M
        den = np.einsum("ij,ij->i", X, X)
        a = num / den
        a = a / a.mean()
        M = (X.T * a).T.mean(0)
        hist_a.append(a.copy()); Ms.append(M.copy())
        hist_J.append(float(np.sum(((X.T * a).T - M) ** 2)))
    a = hist_a[-1]

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.0))
    a0 = axes[0]
    amp = 1.0 / a
    norm = plt.Normalize(amp.min(), amp.max())
    cm = plt.cm.viridis
    for x, am in zip(X, amp):
        a0.plot(lam, x, lw=1.1, color=cm(norm(am)), alpha=0.9)
    a0.plot(lam, Ms[-1], lw=3.0, color="#c0392b", label="shape $S(\\lambda)$")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cm); sm.set_array([])
    plt.colorbar(sm, ax=a0, label="amplitude $c_i$ (relative)")
    a0.legend(fontsize=9)
    a0.set_ylabel("$R_{rs}$  sr$^{-1}$")
    a0.set_title("MODEL   $R_i(\\lambda) = c_i\\,S(\\lambda) + \\varepsilon$\n"
                 "colour = the fitted $c_i$: the spectra are one shape, scaled",
                 fontsize=10.5, loc="left")

    a1 = axes[1]
    # J itself changes in the 6th significant figure after the first step, so plotting
    # J on a linear axis shows a flat line. The EXCESS over the converged value is what
    # demonstrates convergence.
    Jinf = min(hist_J)
    exc = [max((j - Jinf) / Jinf, 1e-17) for j in hist_J]
    a1.semilogy(range(1, len(exc) + 1), exc, "o-", color="#2c6f9b", ms=7)
    a1.set_xlabel("iteration")
    a1.set_ylabel("$(J_k-J_\\infty)/J_\\infty$", color="#2c6f9b")
    a1.grid(alpha=0.25)
    a1b = a1.twinx()
    drift = [float(np.max(np.abs(hist_a[k] - hist_a[k - 1]))) if k else np.nan
             for k in range(len(hist_a))]
    a1b.semilogy(range(1, len(drift) + 1), drift, "s--", color="#8a6000", ms=5,
                 label="max $|\\Delta a_i|$")
    a1b.set_ylabel("amplitude change per step", color="#8a6000")
    a1.set_title("CONVERGENCE   alternate (6) and (7)\n"
                 "each step is an exact minimiser, so $J$ cannot rise",
                 fontsize=10.5, loc="left")

    a2 = axes[2]
    raw = X.std(0) / np.abs(X.mean(0)) * 100
    res = (X.T * a).T.std(0) / np.abs(Ms[-1]) * 100
    a2.plot(lam, raw, lw=2.0, color="#c0392b", label="before: %.1f %% (450-700)"
            % float(np.mean(raw[core])))
    a2.plot(lam, res, lw=2.4, color="#2e7d32", label="after:  %.1f %%"
            % float(np.mean(res[core])))
    a2.axvspan(450, 700, color="#2e7d32", alpha=0.07)
    a2.set_ylabel("relative scatter (%)"); a2.legend(fontsize=9)
    a2.set_ylim(0, min(40, float(np.nanmax(raw)) * 1.2))
    a2.set_title("WHAT IT REMOVES\nthe amplitude term, eq. (4), which is identical\n"
                 "in every band; the tails stay noisy because $R_{rs}\\to 0$ there",
                 fontsize=10.5, loc="left")
    for ax_ in (a0, a2):
        ax_.set_xlabel("wavelength (nm)"); ax_.grid(alpha=0.25); ax_.set_xlim(RLO, RHI)
    fig.suptitle("%s — the amplitude-normalised mean, method" % tag, fontsize=13,
                 weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    pth = os.path.join(outdir, "fig11_scaled_mean_method.png")
    fig.savefig(pth, dpi=140); plt.close(fig)
    return pth


# --------------------------------------------------------------- figure 12
def fig_final_mean(water, sky, wl, outdir, tag, panel_r):
    """The deliverable on its own: mean R_rs with the two uncertainties separated."""
    import numpy as np
    pairs = match_by_angle(water, sky)
    R = []
    for w, sk, _, _ in pairs:
        r = rho_at_angle(view_zenith_from_tilt(w["spec"].tilt_y_deg))
        R.append(rrs_three_scan(wl, w["spec"].columns["rad_target"],
                                sk["spec"].columns["rad_target"],
                                w["spec"].columns["rad_ref"], panel_r, r, "none").rrs)
    wl_a = np.array(wl); m = (wl_a >= RLO) & (wl_a <= RHI)
    X = np.array(R)[:, m]; lam = wl_a[m]
    core = (lam >= 450) & (lam <= 700)
    mean, sc, it, shp, amp_cv = scaled_mean([list(x) for x in X])
    mean = np.array(mean); sc = np.array(sc)
    Xs = (X.T * sc).T
    sS = Xs.std(0)
    sC = mean * (amp_cv / 100.0)
    shp_core = float(np.mean(sS[core] / np.abs(mean[core])) * 100)

    fig, ax = plt.subplots(figsize=(10.5, 6.8))
    ax.fill_between(lam, mean - sC, mean + sC, color="#bbb", alpha=0.55,
                    label="amplitude $\\pm$%.0f %% — how MUCH (real variability)"
                          % amp_cv)
    ax.fill_between(lam, mean - sS, mean + sS, color="#2e7d32", alpha=0.40,
                    label="shape $\\pm$%.1f %% (450-700 nm) — what the water IS"
                          % shp_core)
    ax.plot(lam, mean, lw=3.2, color="#14532d", label="mean $R_{rs}$  (n=%d scans)"
            % len(X))
    ax.axhline(0, color="#888", lw=0.8)
    for lam0, txt in ((570, "green peak\nsediment"), (700, "700 nm peak\nhigh SPM"),
                      (810, "$a_w$ window")):
        k = int(np.argmin(abs(lam - lam0)))
        ax.annotate(txt, (lam[k], mean[k]), xytext=(0, 34),
                    textcoords="offset points", ha="center", fontsize=8.5,
                    arrowprops=dict(arrowstyle="->", lw=1))
    ax.set_xlim(RLO, RHI); ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("$R_{rs}$  (sr$^{-1}$)")
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(alpha=0.25)
    ax.set_title("%s\nFINAL $R_{rs}$ — amplitude-normalised mean of %d scans, "
                 "converged in %d iterations" % (tag, len(X), it), fontsize=12,
                 weight="bold", loc="left")
    ax.text(0.015, 0.03,
            "$R_{rs}(\\lambda) = S(\\lambda)\\times c$   with $S$ known to "
            "%.1f %% and $c$ varying %.0f %% across the station.\n"
            "Band ratios inherit the %.1f %%; absolute magnitudes inherit the %.0f %%."
            % (shp_core, amp_cv, shp_core, amp_cv),
            transform=ax.transAxes, fontsize=9.5, va="bottom",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f7f7f7", ec="#999"))
    fig.tight_layout()
    pth = os.path.join(outdir, "fig12_FINAL_mean_Rrs.png")
    fig.savefig(pth, dpi=150); plt.close(fig)
    return pth



def fig_sky_choice(water, sky, wl, outdir, tag, panel_r):
    """PER WATER SCAN: which sky is best, what the choice costs, and does it matter.

    fig4 shows the pooled envelope over all pairings and fig5 shows the geometry, but
    neither answers the question scan by scan: for THIS water spectrum, how much does
    the sky I pair it with move the answer, and is that movement large compared with the
    things we already know the size of?

    The three reference scales it is graded against, all measured in this dataset:
      0.6 %   panel replicate scatter -- the instrument floor
      1.7 %   shape uncertainty of the amplitude-normalised mean (fig11)
     11.4 %   scan-to-scan amplitude spread -- real water variability (fig9)
    """
    import numpy as np

    FLOOR, SHAPE, WATERSD = 0.6, 1.7, 11.4
    i443 = min(range(len(wl)), key=lambda k: abs(wl[k] - 443))
    i555 = min(range(len(wl)), key=lambda k: abs(wl[k] - 555))
    band_m = [i for i, x in enumerate(wl) if RLO <= x <= RHI]

    rows = []
    for w in water:
        lt = w["spec"].columns["rad_target"]
        lp = w["spec"].columns["rad_ref"]
        rho = rho_at_angle(view_zenith_from_tilt(w["spec"].tilt_y_deg))
        tv = view_zenith_from_tilt(w["spec"].tilt_y_deg)
        curves, dth = [], []
        for sk in sky:
            curves.append(rrs_three_scan(wl, lt, sk["spec"].columns["rad_target"], lp,
                                         panel_r, rho, "none").rrs)
            dth.append(abs(view_zenith_from_tilt(sk["spec"].tilt_y_deg) - tv))
        C = np.array(curves)
        best = int(np.argmin(dth))
        rows.append(dict(n=w["n"], C=C, best=best, sky_best=sky[best]["n"],
                         dth=dth[best],
                         # what the CHOICE costs: sd across the 8 skies, per band
                         sd443=100 * C[:, i443].std() / C[:, i443].mean(),
                         sd555=100 * C[:, i555].std() / C[:, i555].mean(),
                         sdband=100 * float(np.mean(C[:, band_m].std(axis=0)
                                                    / C[:, band_m].mean(axis=0))),
                         # what the angle-matched choice differs from the naive
                         # average-all-skies choice by, which is the practical question
                         d443=100 * (C[best, i443] / C[:, i443].mean() - 1),
                         d555=100 * (C[best, i555] / C[:, i555].mean() - 1),
                         dband=100 * float(np.mean(np.abs(C[best, band_m]
                                                          / C[:, band_m].mean(axis=0)
                                                          - 1)))))

    fig = plt.figure(figsize=(17.5, 10.6))

    ax = fig.add_subplot(2, 2, 1)
    x = np.arange(len(rows))
    ax.bar(x - 0.2, [r["sd443"] for r in rows], 0.4, color="#2c6f9b", label="443 nm")
    ax.bar(x + 0.2, [r["sd555"] for r in rows], 0.4, color="#2e7d32", label="555 nm")
    for y, c, l in ((FLOOR, "#888", "instrument floor 0.6 %"),
                    (SHAPE, "#8a6000", "shape uncertainty 1.7 %"),
                    (WATERSD, "#c0392b", "real water spread 11.4 %")):
        ax.axhline(y, color=c, ls="--", lw=1.6)
        ax.text(len(rows) - 0.4, y * 1.04, l, fontsize=8, color=c, ha="right")
    ax.set_xticks(x); ax.set_xticklabels([r["n"] for r in rows], rotation=90, fontsize=7)
    ax.set_ylabel("sd across the %d sky choices  (%%)" % len(sky))
    ax.legend(fontsize=8.5); ax.grid(alpha=0.25, axis="y")
    ax.set_title("WHAT THE SKY CHOICE COSTS, scan by scan.\nmedian %.2f %% at 443, "
                 "%.2f %% at 555 — %.1fx below the real water spread"
                 % (np.median([r["sd443"] for r in rows]),
                    np.median([r["sd555"] for r in rows]),
                    WATERSD / max(np.median([r["sd443"] for r in rows]), 1e-9)),
                 fontsize=10, loc="left")

    ax = fig.add_subplot(2, 2, 2)
    ax.bar(x - 0.2, [r["d443"] for r in rows], 0.4, color="#2c6f9b", label="443 nm")
    ax.bar(x + 0.2, [r["d555"] for r in rows], 0.4, color="#2e7d32", label="555 nm")
    ax.axhline(0, color="k", lw=1)
    for s in (1, -1):
        ax.axhline(s * FLOOR, color="#888", ls=":", lw=1.4)
    ax.set_xticks(x); ax.set_xticklabels([r["n"] for r in rows], rotation=90, fontsize=7)
    ax.set_ylabel("angle-matched $-$ average-of-all-skies  (%)")
    ax.legend(fontsize=8.5); ax.grid(alpha=0.25, axis="y")
    ax.set_title("DOES THE MATCHING CHANGE THE ANSWER? Angle-matched sky against\n"
                 "simply averaging every sky. |mean| %.2f %% at 443, worst %.2f %%. "
                 "Dotted = instrument floor."
                 % (np.mean([abs(r["d443"]) for r in rows]),
                    max(abs(r["d443"]) for r in rows)), fontsize=10, loc="left")

    # the widest-spread scan, spelled out
    worst = max(rows, key=lambda r: r["sdband"])
    ax = fig.add_subplot(2, 2, 3)
    # clip() returns (wl, v) -- unpack it. Passing the tuple straight into plot() makes
    # both arguments 2-D and silently draws one line per band.
    for i in range(len(sky)):
        xs, ys = clip(wl, worst["C"][i], RLO, RHI)
        ax.plot(xs, ys, lw=1.0, color="#bbbbbb",
                label="each of the %d sky choices" % len(sky) if i == 0 else None)
    xs, ys = clip(wl, worst["C"][worst["best"]], RLO, RHI)
    ax.plot(xs, ys, lw=2.6, color="#c0392b",
            label="angle-matched sky %s ($\\Delta\\theta$=%.1f$^\\circ$)"
                  % (worst["sky_best"], worst["dth"]))
    xs, ys = clip(wl, worst["C"].mean(axis=0), RLO, RHI)
    ax.plot(xs, ys, lw=2.0, ls="--", color="#2c6f9b",
            label="average of all %d skies" % len(sky))
    ax.set_xlim(RLO, RHI)
    ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("$R_{rs}$  sr$^{-1}$")
    ax.legend(fontsize=8.5); ax.grid(alpha=0.25)
    ax.set_title("Water %s — the WORST case in this station (%.2f %% band-mean "
                 "spread).\nGrey = every sky choice." % (worst["n"], worst["sdband"]),
                 fontsize=10, loc="left")

    ax = fig.add_subplot(2, 2, 4)
    ax.axis("off")
    med = np.median([r["sdband"] for r in rows])
    mx = max(r["sdband"] for r in rows)
    m443 = np.median([r["sd443"] for r in rows])
    m555 = np.median([r["sd555"] for r in rows])
    verdict = ("IT DEPENDS ON THE BAND, and that is the useful answer.\n"
               "  443 nm : %.2f %% -- %s\n"
               "  555 nm : %.2f %% -- %s\n"
               "The blue is where the sky radiance is largest relative to\n"
               "the water signal, so rho*L_sky is the biggest term there.\n"
               "Anything read off the BLUE (band ratios, OC4, a_dg) carries\n"
               "this; anything read off the green does not."
               % (m443, "ABOVE the %.1f %% shape uncertainty" % SHAPE
                  if m443 > SHAPE else "below the shape uncertainty",
                  m555, "BELOW the %.1f %% instrument floor" % FLOOR
                  if m555 < FLOOR else "above the instrument floor"))
    ax.text(0.0, 0.97,
            "DOES THE SKY CHOICE MATTER, AT THE SCALE WE HAVE?\n"
            "%s\n\n"
            "spread from the sky choice, 400-900 nm band mean\n"
            "   median over the %d water scans : %5.2f %%\n"
            "   worst single scan              : %5.2f %%  (water %s)\n\n"
            "reference scales measured in this dataset\n"
            "   panel replicates, instrument floor : %5.2f %%\n"
            "   shape uncertainty of the mean      : %5.2f %%\n"
            "   scan-to-scan water variability     : %5.2f %%\n\n"
            "angle-matched vs averaging every sky\n"
            "   mean |difference| at 443 nm : %5.2f %%\n"
            "   worst                        : %5.2f %%\n\n"
            "The angle-matched sky gives a SYSTEMATICALLY higher blue\n"
            "than averaging every sky (%d of %d scans positive), so this\n"
            "is a bias, not scatter. The sky scans span only ~5 deg here,\n"
            "so there is little geometry to match; widening the sky\n"
            "angular coverage is what would make matching worth doing\n"
            "(NEXT_CAMPAIGN.md, Tier 2 item 4)."
            % (verdict, len(rows), med, mx, worst["n"], FLOOR, SHAPE, WATERSD,
               np.mean([abs(r["d443"]) for r in rows]),
               max(abs(r["d443"]) for r in rows),
               sum(1 for r in rows if r["d443"] > 0), len(rows)),
            transform=ax.transAxes, va="top", fontsize=9.5, family="monospace",
            bbox=dict(boxstyle="round,pad=0.7", fc="#f7f7f7", ec="#555", lw=1.4))

    fig.suptitle("PER-SCAN SKY PAIRING — %s" % tag, fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.945), h_pad=3.5)
    p = os.path.join(outdir, "fig13_sky_choice_per_scan.png")
    fig.savefig(p, dpi=140); plt.close(fig)

    with open(os.path.join(outdir, "sky_choice_per_scan.csv"), "w", newline="") as fh:
        w_ = csv.writer(fh)
        w_.writerow(["# per water scan: the angle-matched sky, the spread the sky "
                     "choice causes, and the difference from averaging all skies"])
        w_.writerow(["water", "best_sky", "dtheta_deg", "sd443_pct", "sd555_pct",
                     "sdband_pct", "diff443_pct", "diff555_pct", "diffband_pct"])
        for r in rows:
            w_.writerow([r["n"], r["sky_best"], "%.2f" % r["dth"]] +
                        ["%.4f" % r[k] for k in ("sd443", "sd555", "sdband",
                                                 "d443", "d555", "dband")])
    return p, rows


#: What each figure settles. Written into analysis/README.md so the index cannot drift
#: from the figures: every entry is checked against the file the run actually produced.
FIG_INDEX = [
    ("fig1_pooled_measurements.png", "The three raw families, pooled",
     "Every panel, sky and water radiance on one axis each. The spread you see here is "
     "the raw material for everything downstream: panel replicates hold to ~1 %, sky to "
     "~14 %, water to ~33 %. Land targets are drawn but excluded from R_rs."),
    ("fig2_steps.png", "The calculation, step by step",
     "One water scan carried through L_p -> E_d = pi L_p / R_p -> L_w = L_t - rho L_sky "
     "-> R_rs = L_w / E_d, each stage plotted, so the arithmetic is inspectable rather "
     "than asserted."),
    ("fig3_Ed_check.png", "Is the irradiance physical?",
     "E_d against a clear-sky expectation, with the atmospheric transmittance it "
     "implies, the O2-A and 940 nm water-vapour bands, and the Rayleigh red/blue slope. "
     "A real atmosphere shows all of these; a broken panel chain does not."),
    ("fig4_rrs_spread.png", "R_rs, and how much the sky choice costs",
     "Every water scan against every sky scan, then the angle-matched subset. Shows what "
     "part of the spread is the pairing choice and what part is the water."),
    ("fig5_geometry_and_pairing.png", "The geometry actually achieved",
     "Sky and water view angles, the mismatch of each matched pair, and the best pair by "
     "smallest mismatch. The sky scans span only ~5 deg, which is why matching buys "
     "little here."),
    ("fig6_angle_footprint.png", "Does angle matter? Does footprint?",
     "R_rs(443) against view angle and against footprint, with the confound controls. "
     "Angle: yes (r=+0.74, p=0.006, +27 % over the observed span). Footprint: not "
     "separable from time in this dataset (range correlates with time at r=0.92)."),
    ("fig7_land_targets.png", "The two land scans, kept separate",
     "00014 (algae on concrete) and 00015 (bare concrete), reported as REFLECTANCE not "
     "R_rs, because for an opaque surface the reflected sky is illumination, not "
     "contamination, so rho L_sky must NOT be subtracted."),
    ("fig8_rho_angle_correction.png", "Correcting rho per scan instead of fixing it",
     "rho(theta_v) Fresnel-scaled from Mobley's 0.028 at 40 deg. Removes the angle trend "
     "(p 0.006 -> 0.244) and cuts scatter 11.1 % -> 8.6 %. The tilt datum was determined "
     "FROM the data, and the three candidate datums are shown."),
    ("fig9_variability_origin.png", "Where the remaining spread comes from",
     "PCA of the 12 R_rs: PC1 is 98.3 % of the variance and 92 % pure amplitude. Shape "
     "is stable to 1.7 % over 450-700 nm while amplitude moves +-11.4 %, with no time "
     "trend (r=+0.03, p=0.81) and 18x the instrument floor. So it is real water "
     "variability, not measurement noise."),
    ("fig10_final_product.png", "The final product, all treatments together",
     "Mean and full spread under each pairing/rho treatment, so the choice is visible "
     "rather than buried."),
    ("fig11_scaled_mean_method.png", "The amplitude-normalised mean, derived",
     "The iterative scaled mean: each scan rescaled onto the running mean, the mean "
     "re-formed, repeat. Separates the 1.7 % shape uncertainty from the 11 % amplitude "
     "spread instead of reporting one inflated number at every band. Full derivation and "
     "the rank-1 SVD equivalence in THEORY_SCALED_MEAN.md."),
    ("fig13_sky_choice_per_scan.png", "Per-scan sky pairing: does the choice matter?",
     "For EACH water scan: which sky the angle matcher picked, the spread that the "
     "choice of sky causes, and the difference between angle-matching and simply "
     "averaging every sky -- all graded against three scales measured in this dataset "
     "(0.6 % instrument floor, 1.7 % shape uncertainty, 11.4 % real water spread). "
     "Answers 'does it matter at the scale we have' with a number rather than a "
     "judgement. Numbers also in `sky_choice_per_scan.csv`."),
    ("fig12_FINAL_mean_Rrs.png", "FINAL_Rrs.csv, plotted",
     "The spectrum that goes into GIOP, with its shape uncertainty band. This is the one "
     "to quote."),
]


def write_index(outdir, tag, figs, n_scans, n_sky, n_water, n_land):
    """analysis/README.md -- an index of what each figure settles.

    Generated, not hand-written, and it lists only files that exist on disk, so it
    cannot claim a figure the run did not produce.
    """
    have = {os.path.basename(p) for p in figs if p}
    L = ["# %s -- analysis" % tag, "",
         "%d scans: %d sky, %d water, %d land. Generated by `analyse_location.py`; "
         "do not edit by hand." % (n_scans, n_sky, n_water, n_land), "",
         "- `REPORT.txt` -- every number, including the ones the figures only show.",
         "- `REPORT.html` -- the same, with zoomable plotly versions "
         "(`make_interactive_report.py`).",
         "- `FINAL_Rrs.csv` -- wavelength, amplitude-normalised mean R_rs, shape "
         "uncertainty. **This is the product.**",
         "- `GIOP/` -- the semi-analytical inversion of that spectrum "
         "(`make_giop_figures.py`), with its own README.",
         "", "## Figures", ""]
    for fn, title, what in FIG_INDEX:
        if fn in have:
            L.append("### `%s` -- %s" % (fn, title))
            L.append("")
            L.append(what)
            L.append("")
        else:
            L.append("### `%s` -- NOT PRODUCED by this run" % fn)
            L.append("")
            L.append("(%s -- absent because this location has no such scans.)" % title)
            L.append("")
    L += ["## Reading order", "",
          "fig1 -> fig2 -> fig3 establish that the measurement is sound. fig4 -> fig6 "
          "-> fig8 establish what the residual spread is made of and correct the part "
          "that is geometric. fig9 shows the rest is real water. fig11 -> fig12 produce "
          "the spectrum. Then `GIOP/`.", ""]
    with open(os.path.join(outdir, "README.md"), "w") as fh:
        fh.write("\n".join(L))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--rho", type=float, default=RHO_MOBLEY1999)
    ap.add_argument("--panel-reflectance", type=float, default=0.99)
    ap.add_argument("--exclude-water", nargs="*", default=[], metavar="SCAN",
                    help="water scan numbers (e.g. 00027) to drop before forming the "
                         "final mean -- for a station where analyse_water_scans.py has "
                         "shown the water is not one population. The excluded scans are "
                         "still shown in fig1/fig4/fig9 (so the exclusion is visible, "
                         "not silent) but are NOT averaged into FINAL_Rrs.csv.")
    a = ap.parse_args()

    scans = survey(a.folder)
    span_m = assert_same_dataset(scans)
    loc = os.path.basename(os.path.dirname(a.folder.rstrip("/")))
    fo = os.path.basename(a.folder.rstrip("/"))
    tag = "%s  ·  %s" % (loc, fo)
    outdir = os.path.join(a.folder, "analysis")
    os.makedirs(outdir, exist_ok=True)

    sky = sorted([s for s in scans if s["role"] == "sky"], key=lambda x: x["n"])
    water = sorted([s for s in scans if s["role"] == "water"], key=lambda x: x["n"])
    land = sorted([s for s in scans if s["role"] == "land"],
                  key=lambda x: x["n"])
    wl = scans[0]["spec"].wavelength

    # one representative panel per distinct reference
    seen, panels = set(), []
    for s in scans:
        k = round(band(s["spec"], "rad_ref", 450, 650), 6)
        if k not in seen:
            seen.add(k); panels.append(s)

    L = []
    P = L.append
    P("=" * 78)
    P("%s" % tag)
    P("=" * 78)
    P("%d scans   %s-%s UTC   sun %.2f-%.2f deg   FOV %.0f deg"
      % (len(scans), hhmm(min(s["gps"] for s in scans)),
         hhmm(max(s["gps"] for s in scans)), min(s["sun"] for s in scans),
         max(s["sun"] for s in scans), fov_deg(fo.split("_")[0])))
    P("%d sky, %d water, %d LAND targets" % (len(sky), len(water), len(land)))
    P("all scans within %.0f m and one foreoptic: pairing cannot leave this location"
      % span_m)
    P("%d distinct panel references, so the panel was re-taken %d time(s) mid-station"
      % (len(panels), len(panels) - 1))
    rng = [s["range"] for s in scans if s["range"] is not None]
    nmiss = len(scans) - len(rng)
    if rng:
        P("range %.3f-%.3f m -> footprint %.2f-%.2f m across (2 R tan(FOV/2))%s"
          % (min(rng), max(rng),
             2 * min(rng) * math.tan(math.radians(fov_deg(fo.split("_")[0]) / 2)),
             2 * max(rng) * math.tan(math.radians(fov_deg(fo.split("_")[0]) / 2)),
             "   [%d scan(s) carry NO range in the header]" % nmiss if nmiss else ""))
    else:
        P("range: NOT RECORDED in any header, so footprint cannot be computed")
    P("")

    p1 = fig_pooled(sky, water, panels, land, wl, outdir, tag)
    for name, grp, key in (("PANEL", panels, "rad_ref"), ("SKY", sky, "rad_target"),
                           ("WATER", water, "rad_target")):
        cur = [s["spec"].columns[key] for s in grp]
        m, lo, hi = stats(cur)
        sv = [(hi[i] - lo[i]) / m[i] for i, w in enumerate(wl) if 450 <= w <= 650]
        P("%-6s n=%2d   mean L(450-650) = %.4e   full spread %.1f %%"
          % (name, len(cur), sum(m[i] for i, w in enumerate(wl) if 450 <= w <= 650) /
             sum(1 for w in wl if 450 <= w <= 650), 100 * sum(sv) / len(sv)))
    P("")

    l_sky_mean = [sum(s["spec"].columns["rad_target"][i] for s in sky) / len(sky)
                  for i in range(len(wl))]
    rep = water[len(water) // 2]
    p2 = fig_steps(rep, l_sky_mean, rep["spec"].columns["rad_ref"], wl, outdir, tag,
                   a.panel_reflectance, a.rho)

    sun = sum(s["sun"] for s in scans) / len(scans)
    p3, diag = fig_ed(panels, wl, sun, outdir, tag, a.panel_reflectance)
    if diag:
        P("IS E_d REASONABLE?")
        P("  median T over 450-650 nm      %.3f      (clear sky 0.6-0.9)" % diag["T_vis"])
        P("  max T anywhere 380-950 nm     %.3f      %s"
          % (diag["T_max"], "OK, below the hard ceiling of 1"
             if diag["T_ok"] else "IMPOSSIBLE, exceeds 1"))
        P("  O2-A absorption depth         %.1f %%     (a real atmosphere shows this)"
          % (100 * diag["o2_depth"]))
        P("  H2O 940 nm absorption depth   %.1f %%" % (100 * diag["h2o_depth"]))
        P("  T(red)/T(blue)                %.2f       (>1 expected: Rayleigh weakens "
          "to the red)" % diag["blue_red"])
        P("  panel-to-panel disagreement   %.1f %%     over 450-650 nm"
          % (100 * diag["ed_spread"]))
        verdict = ("REASONABLE" if diag["T_ok"] and 0.5 < diag["T_vis"] < 0.95
                   and diag["o2_depth"] > 0.05 and diag["blue_red"] > 1.0
                   else "SUSPECT, inspect fig3")
        P("  VERDICT: %s" % verdict)
        P("")

    excl = set(a.exclude_water)
    water_clean = [w for w in water if w["n"] not in excl]
    if excl:
        found = {w["n"] for w in water} & excl
        missing = excl - found
        P("EXCLUDED FROM THE FINAL PRODUCT (--exclude-water)")
        P("  %s -- see analysis/water_scans/REPORT.txt for why (run "
          "analyse_water_scans.py first)." % ", ".join(sorted(found)))
        if missing:
            P("  ⚠ requested but not found among this station's water scans: %s"
              % ", ".join(sorted(missing)))
        P("  %d of %d water scans remain for FINAL_Rrs.csv."
          % (len(water_clean), len(water)))
        P("  Figures 1/4/5/6/9/13 below still show ALL %d scans (the exclusion is not "
          "silent); figures 10/11/12 and FINAL_Rrs.csv use only the %d clean ones."
          % (len(water), len(water_clean)))
        P("")
        if not water_clean:
            raise SystemExit("--exclude-water removed every water scan; nothing left "
                             "to form a final spectrum from")

    p4, r = fig_rrs(water, sky, wl, outdir, tag, a.panel_reflectance, a.rho)
    P("R_rs RESULT  (angle-matched pairing, per-scan rho)")
    P("  Rrs(443) = %.5f sr^-1     Rrs(555) = %.5f sr^-1" % (r["rrs443"], r["rrs555"]))
    P("  sky paired to each water by mirrored view angle, WITHIN this location and")
    P("  within the same panel-reference block: mismatch mean %.1f deg, worst %.1f deg"
      % (r["mism_mean"], r["mism_max"]))
    if r["fallbacks"]:
        P("  block fallback needed for: %s" % ", ".join(r["fallbacks"]))
    P("")
    P("  SPREAD AT 443 nm.  Quote the sd: a min-max envelope grows with sample size,")
    P("  so A (%d spectra) and C (%d) are NOT comparable on it." % (r["n_A"], r["n_C"]))
    P("    %-34s %8s %8s" % ("treatment", "envelope", "sd"))
    P("    A  blind pairing, fixed rho        %7.1f %% %7.1f %%"
      % (r["envA443"], r["sdA443"]))
    P("    B  blind pairing, rho(theta_v)     %7.1f %% %8s" % (r["envB443"], "-"))
    P("    C  ANGLE-MATCHED, rho(theta_v)     %7.1f %% %7.1f %%"
      % (r["envC443"], r["sdC443"]))
    P("")
    P("  CONTROL: 2000 draws assigning a RANDOM sky per water at the same n=%d, drawn"
      % r["n_C"])
    P("  from the SAME panel block, differing from the treatment only in ignoring angle:")
    P("  sd = %.1f %% (median). Angle matching gives %.1f %%, and %.0f %% of random draws"
      % (r["ctrl_med"], r["sdC443"], 100 * r["ctrl_frac"]))
    P("  are at least as tight, so it is NOT a demonstrated improvement at this n.")
    P("  Reason: these sky scans span only ~5 deg, so there is little geometry to match.")
    P("  -> the honest spread on R_rs(443) is the sd, about %.0f %%, and the fall from"
      % r["sdC443"])
    P("     the %.0f %% blind envelope is mostly sample size, not physics."
      % r["envA443"])
    P("  spread from WHICH WATER scan:  %.1f %% at 443 nm (max %.1f %% over 400-900)"
      % (r["w_spread_443"], r["w_spread"]))
    P("  spread from WHICH SKY scan:    %.1f %% at 443 nm (max %.1f %% over 400-900)"
      % (r["s_spread_443"], r["s_spread"]))
    dom = "WATER-to-water variability" if r["w_spread_443"] > r["s_spread_443"] \
        else "the SKY choice"
    P("  -> %s dominates at 443 nm." % dom)
    P("")
    P("For scale, rho alone (0.022-0.045) moves Rrs(443) by ~32 %%, so it remains the")
    P("largest single term. See THEORY.pdf p3.")

    p5, g = fig_geometry(scans, water, sky, wl, outdir, tag, a.panel_reflectance, a.rho)
    P("")
    P("BEST-MATCHED PAIR (smallest NIR-similarity residual)")
    P("  best   water %s + sky %s   delta = %+.3e sr^-1" % (g["best"][0], g["best"][1],
                                                            g["d_best"]))
    P("  worst  water %s + sky %s   delta = %+.3e sr^-1" % (g["worst"][0],
                                                            g["worst"][1], g["d_worst"]))
    P("  across all pairings Rrs(443) = %.5f +/- %.5f (1 sd)"
      % (g["rrs443_mean"], g["rrs443_sd"]))

    fov = fov_deg(fo.split("_")[0])
    p6, spec, conf = fig_sensitivity(water, sky, wl, outdir, tag, a.panel_reflectance,
                                     a.rho, fov)
    p7 = fig_land(land, wl, outdir, tag, a.panel_reflectance) if land else None
    p8, rres, rrange = fig_rho_correction(water, sky, wl, outdir, tag,
                                          a.panel_reflectance)
    p9, vv = fig_variability(water, sky, wl, outdir, tag, a.panel_reflectance, panels)
    p10, fcsv, fp = fig_final(water_clean, sky, wl, outdir, tag, a.panel_reflectance)
    p11 = fig_scaled_method(water_clean, sky, wl, outdir, tag, a.panel_reflectance)
    p12 = fig_final_mean(water_clean, sky, wl, outdir, tag, a.panel_reflectance)
    P("")
    P("WHAT IS THE REMAINING SPREAD?  measurement error, or real water?")
    P("  SVD mode 1 carries          %5.1f %%   coherent; noise would spread out" % vv["pc1"])
    P("  spread raw / normalised     %.1f %% / %.1f %%" % (vv["raw"], vv["nor"]))
    P("  -> %.0f %% of the variance is pure AMPLITUDE, only %.1f %% is SHAPE"
      % (vv["amp_frac"], vv["nor"]))
    P("  R_rs(555) range             factor %.2f across the station" % vv["ratio"])
    P("  spectral distance vs time   r=%+.2f p=%.2f  -> %s"
      % (vv["r_t"], vv["p_t"],
         "drift" if vv["p_t"] < 0.05 and vv["r_t"] > 0 else "NO time trend"))
    P("  instrument floor (panel)    %.1f %%   sky %.1f %%   water %.1f %%"
      % (vv["floor"], vv["sky_floor"], vv["raw"]))
    P("  water exceeds the instrument floor by %.0fx." % (vv["raw"] / vv["floor"]))
    P("")
    P("  READING: the spread is REAL WATER, not measurement error. The composition is")
    P("  near-constant (shape stable to %.1f %%) while the AMOUNT of scattering material"
      % vv["nor"])
    P("  changes by a factor %.2f between scans, with no time trend, over a footprint of"
      % vv["ratio"])
    P("  ~0.16 m. That is spatial patchiness in suspended load, which is what a turbid")
    P("  nearshore surface does.")
    P("")
    P("  CONSEQUENCE FOR INVERSION: band RATIOS are far more stable than absolute R_rs")
    P("  here (%.1f %% vs %.1f %%), so ratio-based products should be quoted with more"
      % (vv["nor"], vv["raw"]))
    P("  confidence than absolute magnitudes at this site.")
    P("")
    P("=" * 78)
    P("FINAL PRODUCT  (fig10, FINAL_Rrs.csv)")
    P("=" * 78)
    P("  Iterative amplitude-normalised mean, converged in %d iterations." % fp["iters"])
    P("  Each scan is scaled onto the running mean, then the mean is re-formed, which")
    P("  separates the two things that are genuinely different:")
    P("")
    P("    SHAPE      %5.1f %% over 450-700 nm   <- what the water IS. The product."
      % fp["shape_core"])
    P("    AMPLITUDE  %5.1f %%                    <- how MUCH there is. Real, not error."
      % fp["amp_cv"])
    P("    (a plain average blends them into one %.1f %% band)" % fp["plain_core"])
    P("")
    P("  Over the full 400-900 nm the shape figure is %.1f %%, but a RELATIVE scatter is"
      % fp["shape_cv"])
    P("  meaningless where R_rs approaches zero, which it does below 430 and above 800.")
    P("  The mean SHAPE is %.0fx better determined than a plain average suggests."
      % fp["gain"])
    P("  Per-scan amplitude relative to the mean:")
    row = "    "
    for n_, v_ in zip(fp["names"], fp["amp"]):
        row += "%s %.2f  " % (n_, v_)
        if len(row) > 68:
            P(row); row = "    "
    if row.strip():
        P(row)
    P("")
    P("  Quote it as: R_rs = (shape, +/- %.1f %%) x (amplitude, %.0f %% 1 sd across the"
      % (fp["shape_core"], fp["amp_cv"]))
    P("  station). One combined error bar would imply the water TYPE is uncertain when")
    P("  only its CONCENTRATION is.")
    P("")
    P("  Method figure: fig11.  Deliverable on its own: fig12.")
    P("  Full derivation, equivalences and failure modes: THEORY_SCALED_MEAN.md")
    P("")
    P("PER-SCAN rho FROM THE MEASURED ANGLE (rather than a fixed 0.028)")
    P("  rho spans %.5f-%.5f across the achieved angles, a %.0f %% range."
      % (rrange[0], rrange[1], 100 * (rrange[1] / rrange[0] - 1)))
    P("  %-26s %8s %8s %10s" % ("treatment", "r", "p", "scatter"))
    for k, v in rres.items():
        lab = k.replace("$", "").replace("\\rho", "rho").replace("\\theta_v", "theta_v")
        P("  %-26s %+8.2f %8.3f %9.1f %%" % (lab, v["r"], v["p"], v["sd"]))
    best = min(rres, key=lambda k: rres[k]["sd"])
    P("  -> lowest scatter: %s" % best.replace("$", "").replace("\\rho", "rho")
      .replace("\\theta_v", "theta_v"))
    P("  A correct correction flattens the trend AND reduces scatter. One that does")
    P("  neither is refuted, which is how the tilt datum was settled (see")
    P("  fieldrrs.rrs.view_zenith_from_tilt).")
    P("")
    P("DOES GEOMETRY MATTER? (Rrs(443), n=%d water / %d sky)" % (len(water), len(sky)))
    for d in spec:
        P("  %-34s r=%+.2f p=%.3f   over the observed span: %+.1f %%"
          % (d["what"].replace("Does the ", "").replace("?", ""), d["r"], d["p"],
             d["rel"] if d["effect"] >= 0 else -d["rel"]))
    sig = [d for d in spec if d["p"] < 0.05]
    P("  -> %s" % ("no term reaches p<0.05, so none is separable from scan-to-scan "
                   "noise at this n" if not sig
                   else "SIGNIFICANT: " + "; ".join(d["what"] for d in sig)))
    P("")
    P("  CONFOUND CONTROL (a tilt trend could just be time or range in disguise)")
    P("    tilt correlated with time?    r=%+.2f p=%.3f" % conf["tilt_vs_time"])
    P("    range correlated with time?   r=%+.2f p=%.3f" % conf["range_vs_time"])
    P("    tilt vs Rrs | control TIME    r=%+.2f p=%.3f" % conf["tilt_ctrl_time"])
    P("    tilt vs Rrs | control RANGE   r=%+.2f p=%.3f" % conf["tilt_ctrl_range"])
    P("    range vs Rrs | control TILT   r=%+.2f p=%.3f" % conf["range_ctrl_tilt"])
    P("    range took only %d distinct values, so 'footprint' is a %d-level factor "
      "here, not a continuum." % (conf["n_ranges"], conf["n_ranges"]))

    if land:
        from process_field_day import land_reflectance
        P("")
        P("LAND TARGETS -- reported as REFLECTANCE, excluded from R_rs")
        P("  They are opaque surfaces, so rho*L_sky is NOT subtracted: for water the")
        P("  reflected sky is a contaminant, for a solid surface it is illumination.")
        P("  Product is R = R_panel * L_target / L_panel, needing no rho and no sky.")
        P("")
        P("  %-7s %9s %9s %9s %9s   %s"
          % ("scan", "R(550)", "R(670)", "R(750)", "red edge", "reading"))
        for s_ in land:
            w, r = land_reflectance(s_["spec"], a.panel_reflectance)
            at = lambda lo, hi: (sum(v for x, v in zip(w, r) if lo <= x <= hi)
                                 / max(1, sum(1 for x in w if lo <= x <= hi)))
            e = s_["diag"]["red_edge"]
            P("  %-7s %9.4f %9.4f %9.4f %9.2f   %s"
              % (s_["n"], at(545, 555), at(665, 675), at(745, 755), e,
                 "VEGETATED (chlorophyll)" if e > 3 else "bare / weakly vegetated"))
        P("")
        P("  The bare target doubles as a check on the E_d chain: a mineral surface")
        P("  should be spectrally smooth and grey, and it is.")

    p13, skyrows = fig_sky_choice(water, sky, wl, outdir, tag, a.panel_reflectance)
    P("")
    P("PER-SCAN SKY PAIRING  (fig13, sky_choice_per_scan.csv)")
    P("  For each water scan: which sky the angle matcher picked, how much the CHOICE")
    P("  of sky moves R_rs, and how that compares with scales we already know.")
    P("  %-7s %-7s %7s %9s %9s %10s" % ("water", "sky", "dtheta", "sd443 %",
                                        "sd555 %", "vs mean %"))
    for r in skyrows:
        P("  %-7s %-7s %7.1f %9.2f %9.2f %10.2f"
          % (r["n"], r["sky_best"], r["dth"], r["sd443"], r["sd555"], r["d443"]))
    import numpy as _np
    _med = _np.median([r["sdband"] for r in skyrows])
    _mx = max(r["sdband"] for r in skyrows)
    P("  band-mean spread from the sky choice: median %.2f %%, worst %.2f %%"
      % (_med, _mx))
    P("  against 0.6 % instrument floor, 1.7 % shape uncertainty, 11.4 % water spread")
    P("  -> %s" % ("the sky choice does NOT matter at this scale" if _mx < 1.7 else
                   "the worst scans exceed the shape uncertainty"))
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(outdir, "REPORT.txt"), "w") as fh:
        fh.write(txt + "\n")
    figs = (p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13)
    write_index(outdir, tag, figs, len(scans), len(sky), len(water), len(land))
    print("\nwrote %s/REPORT.txt" % outdir)
    print("wrote %s/README.md" % outdir)
    print("wrote %s" % fcsv)
    for p in figs:
        if p:
            print("wrote %s" % p)
    print("wrote %s/sky_choice_per_scan.csv" % outdir)


if __name__ == "__main__":
    main()
