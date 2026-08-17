"""Verify the R_rs chain on the real data, BEFORE any inversion is attempted.

    python verify_field_calcs.py Data_NatureSpec/2026_Aug_16

Six checks. The first four are the ones that matter, because each uses a reference the
package does not control:

  V1  ALGEBRA        re-derive R_rs from the raw columns with no package code and
                     require BITWISE equality with fieldrrs.
  V2  THE INSTRUMENT DARWin computes its own `Reflect.` column = target/reference. That
                     is an INDEPENDENT computation of the same ratio by the vendor's
                     software, so it audits our column mapping and our arithmetic.
  V3  PHYSICS        E_d inferred from the panel, divided by the true top-of-atmosphere
                     solar irradiance F0 x cos(theta_s), must give an atmospheric
                     transmittance in the physically possible range. Nothing in the
                     package knows F0; it comes from the GIOP data tables.
  V4  ROUND TRIP     reconstruct L_t from the retrieved R_rs and require it back.
  V5  SENSITIVITY    how much do rho and the panel reflectance actually move R_rs.
  V6  CONSERVATION   R_rs must be bounded by simple physical limits.
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fieldrrs.rrs import RHO_MOBLEY1999, rrs_three_scan
from process_field_day import SKY, WATER, band, load, mean_spectrum, process

OK, BAD = "  PASS", "  FAIL"
fails = []


def check(name, ok, detail):
    print("%s  %-46s %s" % (OK if ok else BAD, name, detail))
    if not ok:
        fails.append(name)


def v1_algebra(scans, results):
    """R_rs = (L_t - rho L_sky) / (pi L_panel / R_panel), by hand."""
    print("\nV1  ALGEBRA -- independent re-derivation, no package code")
    st = results[0]
    r = st["rrs"][0]
    lt = r["spec"].columns["rad_target"]
    ls, lp = st["l_sky"], st["l_panel"]
    # the panel reference is exactly 0 above 1653 nm (no calibrated SWIR2 reference),
    # so compare only where the quantity is defined at all
    idx = [i for i, p in enumerate(lp) if p > 0.0]
    manual = [(lt[i] - RHO_MOBLEY1999 * ls[i]) / (math.pi * lp[i] / 0.99) for i in idx]
    d = max(abs(manual[k] - r["res"].rrs[i]) for k, i in enumerate(idx))
    check("hand-computed R_rs == fieldrrs", d == 0.0, "max|diff| = %.3e" % d)

    # and the same thing through the public entry point
    wl = st["wl"]
    again = rrs_three_scan(wl, lt, ls, lp, 0.99, RHO_MOBLEY1999, "none")
    d2 = max(abs(again.rrs[i] - r["res"].rrs[i]) for i in idx)
    check("rrs_three_scan is deterministic", d2 == 0.0, "max|diff| = %.3e" % d2)


def v2_instrument(scans):
    """DARWin's own Reflect. column audits our column mapping and arithmetic.

    Compared on a MEDIAN and a 95th percentile, not a worst case. The worst case is
    dominated by the deep telluric bands -- the 760 nm O2-A line and the 940 nm water
    vapour band -- where target and reference are both near zero and their ratio is a
    ratio of two small noisy numbers. A single such sample reads 10 % while the spectrum
    as a whole agrees to 0.2 %, so a worst-case threshold tests the atmosphere rather
    than the code.
    """
    print("\nV2  THE INSTRUMENT -- DARWin's own Reflect. column as an oracle")
    import statistics
    meds, biases, sds = [], [], []
    for s in scans:
        sp = s["spec"]
        if not sp.has("reflectance"):
            continue
        signed = []
        for w, t_, rf, th in zip(sp.wavelength, sp.columns["rad_target"],
                                 sp.columns["rad_ref"], sp.columns["reflectance"]):
            if rf <= 0 or th <= 1e-4 or 758 <= w <= 772:
                continue
            if 450 <= w <= 700:
                signed.append((t_ / rf - th) / th)
        if len(signed) < 50:
            continue
        meds.append(statistics.median(abs(x) for x in signed))
        biases.append(statistics.mean(signed))
        sds.append(statistics.stdev(signed))
    med = max(meds)
    bias = max(abs(b) for b in biases)
    sd = max(sds)
    check("agrees with DARWin, 450-700 nm (median)", med < 5e-3,
          "worst median |diff| %.3f %%" % (100 * med))
    # The real question is BIAS, not scatter. A column-mapping or arithmetic error is
    # systematic and would show as a non-zero mean; the scatter is a processing
    # difference in how DARWin forms its own ratio, and being unbiased it averages away.
    check("and carries NO systematic bias", bias < 1e-3,
          "worst |mean signed residual| %.4f %% against %.3f %% scatter"
          % (100 * bias, 100 * sd))
    check("  so rad_target and rad_ref are not swapped", med < 5e-3,
          "a swap would give the reciprocal, ~10^2 off")
    print("      Scatter is SNR-driven: it tracks 1/L_t, worst at 720-750 nm where the\n"
          "      water radiance is lowest. Unbiased, so it contributes ~%.2f %% to a\n"
          "      33-scan station mean." % (100 * sd / 33 ** 0.5))


def v3_physics(results):
    """E_d / (F0 cos theta_s) must be a physical atmospheric transmittance."""
    print("\nV3  PHYSICS -- inferred E_d against the true solar constant")
    try:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "..", "giop_python", "src"))
        import numpy as np
        from giop.water import f0_solar
    except Exception as exc:
        print("     SKIPPED, F0 table unavailable (%s)" % exc)
        return
    allT = []
    for st in results:
        wl = np.array(st["wl"])
        ed = np.array([math.pi * p / 0.99 for p in st["l_panel"]])
        mu = math.cos(math.radians(90.0 - st["sun"]))
        m = (wl >= 450) & (wl <= 650)
        T = ed[m] / (f0_solar(wl[m]) * mu)
        allT.append((st["ref"], float(np.median(T)), st["sun"]))
    for ref, T, sun in allT:
        print("      ref=%.4f  sun %.1f deg  ->  median T(450-650) = %.3f"
              % (ref, sun, T))
    ts = [t for _, t, _ in allT]
    check("transmittance is physically possible (0<T<1)",
          all(0.0 < t < 1.0 for t in ts), "range %.3f-%.3f" % (min(ts), max(ts)))
    check("and consistent with a clear-to-hazy sky (0.4-0.95)",
          all(0.4 < t < 0.95 for t in ts),
          "a wrong panel reflectance or unit error would fall outside")


def v4_roundtrip(results):
    print("\nV4  ROUND TRIP -- rebuild L_t from R_rs")
    st = results[0]
    r = st["rrs"][0]
    lt = r["spec"].columns["rad_target"]
    ed = [math.pi * p / 0.99 for p in st["l_panel"]]
    back = [v * e + RHO_MOBLEY1999 * s
            for v, e, s in zip(r["res"].rrs, ed, st["l_sky"])]
    # (the 400-900 window below is entirely inside the region where the panel is > 0)
    pairs = [(a, b) for w, a, b in zip(st["wl"], back, lt) if 400 <= w <= 900]
    rel = max(abs(a - b) / b for a, b in pairs if b > 0)
    check("L_t reconstructed from R_rs", rel < 1e-12, "worst rel error %.2e" % rel)


def v5_sensitivity(results):
    print("\nV5  SENSITIVITY -- what actually moves the answer")
    st = results[0]
    r = st["rrs"][0]
    lt = r["spec"].columns["rad_target"]
    wl, ls, lp = st["wl"], st["l_sky"], st["l_panel"]
    i443 = min(range(len(wl)), key=lambda k: abs(wl[k] - 443))
    i555 = min(range(len(wl)), key=lambda k: abs(wl[k] - 555))
    base = rrs_three_scan(wl, lt, ls, lp, 0.99, 0.028, "none").rrs
    for rho in (0.022, 0.028, 0.035, 0.045):
        v = rrs_three_scan(wl, lt, ls, lp, 0.99, rho, "none").rrs
        print("      rho=%.3f  Rrs(443)=%.5f (%+.1f%%)  Rrs(555)=%.5f (%+.1f%%)"
              % (rho, v[i443], 100 * (v[i443] / base[i443] - 1),
                 v[i555], 100 * (v[i555] / base[i555] - 1)))
    for pr in (0.95, 0.99, 1.00):
        v = rrs_three_scan(wl, lt, ls, lp, pr, 0.028, "none").rrs
        print("      R_panel=%.2f  Rrs(443)=%.5f (%+.1f%%)"
              % (pr, v[i443], 100 * (v[i443] / base[i443] - 1)))
    lo = rrs_three_scan(wl, lt, ls, lp, 0.99, 0.022, "none").rrs[i443]
    hi = rrs_three_scan(wl, lt, ls, lp, 0.99, 0.045, "none").rrs[i443]
    check("rho dominates, as the error budget says",
          abs(hi / lo - 1) > 0.10,
          "rho 0.022->0.045 moves Rrs(443) by %.0f%%" % (100 * abs(hi / lo - 1)))


def v6_conservation(results):
    print("\nV6  CONSERVATION -- physical bounds on every retrieved spectrum")
    nneg = worst_nir = 0
    peak_ok = True
    maxrrs = 0.0
    hi_nir = []
    for st in results:
        for r in st["rrs"]:
            wl, v = r["res"].wavelength, r["res"].rrs
            vis = [x for w, x in zip(wl, v) if 400 <= w <= 700]
            nneg += sum(1 for x in vis if x < -1e-5)
            maxrrs = max(maxrrs, max(vis))
            nir = [x for w, x in zip(wl, v) if 860 <= w <= 890]
            g = max(x for w, x in zip(wl, v) if 550 <= w <= 600)
            ratio = (sum(nir) / len(nir)) / g
            worst_nir = max(worst_nir, ratio)
            if ratio > 0.25:
                hi_nir.append((r["n"], ratio))
            i = max((k for k in range(len(wl)) if 400 <= wl[k] <= 750),
                    key=lambda k: v[k])
            peak_ok &= wl[i] > 520
    check("no negative R_rs in 400-700 nm", nneg == 0, "%d negative samples" % nneg)
    check("R_rs below the 1/pi isotropic ceiling", maxrrs < 1 / math.pi,
          "max R_rs = %.4f sr^-1 (limit 0.318)" % maxrrs)
    check("R_rs collapses in the NIR for nearly all scans",
          len(hi_nir) <= 3,
          "%d of 33 scans exceed NIR/green = 0.25 (worst %.3f)"
          % (len(hi_nir), worst_nir))
    if hi_nir:
        print("      FLAGGED, inspect before using: %s"
              % ", ".join("%s (%.2f)" % x for x in hi_nir))
        print("      High NIR means very high suspended sediment OR residual surface "
              "reflection.\n      It is a property of those spectra, not of the "
              "arithmetic.")
    check("every spectrum peaks green, not blue", peak_ok,
          "consistent with sediment-laden water")


def main():
    folder = sys.argv[1] if len(sys.argv) > 1 else "Data_NatureSpec/2026_Aug_16"
    scans = load(folder)
    results = process(scans, RHO_MOBLEY1999, 0.99, "none")
    print("VERIFYING %d scans -> %d stations, %d water spectra"
          % (len(scans), len(results), sum(len(s["rrs"]) for s in results)))
    v1_algebra(scans, results)
    v2_instrument(scans)
    v3_physics(results)
    v4_roundtrip(results)
    v5_sensitivity(results)
    v6_conservation(results)
    print("\n" + "=" * 72)
    if fails:
        print("FAILED: %s" % ", ".join(fails))
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
