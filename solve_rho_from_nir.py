"""Solve for rho per scan from the NIR similarity constraint, instead of assuming it.

    python solve_rho_from_nir.py Data_NatureSpec/2026_Aug_16

WHY. `residual_correction(method='nir_similarity')` already assumes Ruddick et al.
(2006): for turbid water R_rs(780)/R_rs(870) = 1.912. It then removes the discrepancy
as a SPECTRALLY FLAT offset. But the dominant thing that discrepancy is made of, a
mis-specified rho, is NOT flat: its residual is rho_err * L_sky(lambda)/E_d(lambda),
and L_sky is strongly blue-weighted by Rayleigh scattering. Correcting an L_sky-shaped
error with a flat offset necessarily under-corrects the blue, which is exactly where
the surface term carries 31-46 % of the raw signal at these stations.

So: use the SAME published assumption, but let it constrain the parameter that has the
right spectral shape. With  R_rs(l) = [L_t(l) - rho*L_sky(l)] / E_d(l)  and
R_rs(780) = alpha*R_rs(870), solve the (linear in rho) equation

    rho = [ L_t(780)/E_d(780) - alpha*L_t(870)/E_d(870) ]
          -----------------------------------------------
          [ L_sky(780)/E_d(780) - alpha*L_sky(870)/E_d(870) ]

This is a MEASUREMENT of rho under the wind and sky that actually occurred, from bands
with real signal-to-noise (unlike the 1.6-2.2 um windows, where L_sky collapses and the
ratio is contaminated by any additive floor or direct-beam glint).

It buys nothing for free: it trades "rho is 0.028" for "the Ruddick ratio is 1.912",
and the second is an assumption too. What it does is put the assumption where the data
constrains it and where the error has a known spectral shape. Reports the blast radius
at 443 nm, the band with the most surface leverage and the least tolerance for error.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from analyse_location import match_by_angle
from fieldrrs.rrs import (RHO_MOBLEY1999, SIMILARITY_780_870, rho_at_angle,
                          rrs_three_scan, view_zenith_from_tilt)
from organize_by_location import survey

STATIONS = [
    ("LOC1_66.89718N_162.60290W/FLENS8_FOV08", "LOC1", 4.73),
    ("LOC2a_66.89677N_162.57953W/FLENS8_FOV08", "LOC2a", 4.54),
    ("LOC2b_66.89677N_162.57953W/FLENS8_FOV08", "LOC2b", 4.54),
    ("LOC3_66.89235N_162.59149W/FIBR15_FOV15", "LOC3-FIBR15", 5.62),
    ("LOC3_66.89235N_162.59149W/FIBR15_FOV15_murky", "LOC3-FIBR15murky", 5.62),
    ("LOC3_66.89235N_162.59149W/FLENS8_FOV08", "LOC3-FLENS8", 5.32),
]


def _at(wl, y, target, half=4.0):
    m = (wl >= target - half) & (wl <= target + half)
    return float(np.nanmedian(np.asarray(y, float)[m]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("day_folder")
    ap.add_argument("--panel-reflectance", type=float, default=0.99)
    a = ap.parse_args()
    byloc = os.path.join(a.day_folder.rstrip("/"), "by_location")
    alpha = SIMILARITY_780_870

    print("=" * 104)
    print("rho SOLVED FROM THE NIR SIMILARITY CONSTRAINT (Ruddick et al. 2006, "
          "alpha = R_rs(780)/R_rs(870) = %.3f)" % alpha)
    print("=" * 104)
    print("%-18s %6s  %14s  %22s  %26s" %
          ("station", "wind", "rho assumed", "rho SOLVED (16-50-84%)", "R_rs(443) shift vs assumed"))
    print("-" * 104)

    for rel, label, wind in STATIONS:
        folder = os.path.join(byloc, rel)
        if not os.path.isdir(folder):
            continue
        scans = survey(folder)
        wl = np.array(scans[0]["spec"].wavelength, float)
        water = sorted([s for s in scans if s["role"] == "water"], key=lambda x: x["n"])
        sky = [s for s in scans if s["role"] == "sky"]
        if not water or not sky:
            continue

        rhos, rho_assumed, d443_flat, d443_rho = [], [], [], []
        for w, sk, _dm, _ in match_by_angle(water, sky):
            lt = w["spec"].columns["rad_target"]
            ls = sk["spec"].columns["rad_target"]
            lp = w["spec"].columns["rad_ref"]
            ed = np.pi * np.asarray(lp, float) / a.panel_reflectance

            t780, t870 = _at(wl, lt, 780), _at(wl, lt, 870)
            s780, s870 = _at(wl, ls, 780), _at(wl, ls, 870)
            e780, e870 = _at(wl, ed, 780), _at(wl, ed, 870)
            num = t780 / e780 - alpha * t870 / e870
            den = s780 / e780 - alpha * s870 / e870
            if not np.isfinite(num) or not np.isfinite(den) or abs(den) < 1e-9:
                continue
            rho_fit = num / den
            rhos.append(rho_fit)

            rho0 = rho_at_angle(view_zenith_from_tilt(w["spec"].tilt_y_deg))
            rho_assumed.append(rho0)

            # R_rs(443) three ways: assumed rho, assumed rho + flat NIR offset, solved rho
            base = rrs_three_scan(wl, lt, ls, lp, a.panel_reflectance, rho0, "none")
            flat = rrs_three_scan(wl, lt, ls, lp, a.panel_reflectance, rho0, "nir_similarity")
            solv = rrs_three_scan(wl, lt, ls, lp, a.panel_reflectance, rho_fit, "none")
            i443 = int(np.argmin(abs(wl - 443)))
            b = base.rrs[i443]
            if b and np.isfinite(b) and b != 0:
                d443_flat.append(100.0 * (flat.rrs[i443] - b) / b)
                d443_rho.append(100.0 * (solv.rrs[i443] - b) / b)

        if not rhos:
            print("%-18s  (no usable pairs)" % label)
            continue
        r = np.array(rhos, float)
        q16, q50, q84 = np.nanpercentile(r, [16, 50, 84])
        print("%-18s %5.1f  %8.4f      %6.4f  %6.4f  %6.4f   flat %+6.1f %%   rho-shaped %+6.1f %%"
              % (label, wind, float(np.nanmedian(rho_assumed)), q16, q50, q84,
                 float(np.nanmedian(d443_flat)), float(np.nanmedian(d443_rho))))

    print("-" * 104)
    print("'flat' = what the shipped nir_similarity correction does to R_rs(443).")
    print("'rho-shaped' = what solving for rho instead does to R_rs(443).")
    print("The GAP between those two columns is the part of the blue the flat")
    print("correction cannot reach, because a rho error is not spectrally flat.")
    print("=" * 104)


if __name__ == "__main__":
    main()
