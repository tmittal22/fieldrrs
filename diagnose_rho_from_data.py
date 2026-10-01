"""Can rho be MEASURED from these scans, instead of assumed, under the wind we had?

    python diagnose_rho_from_data.py Data_NatureSpec/2026_Aug_16

Three diagnostics, all computed from data already on disk, no new field work:

1. SKY STRUCTURE INDEX,  pi*L_sky/E_d  =  R_panel * L_sky / L_panel.
   The whole reason wind perturbs rho is that a roughened surface samples a WIDER
   PART OF THE SKY DOME. Under a perfectly uniform (overcast) sky that resampling
   returns the same radiance whatever the slope distribution does, so rho is nearly
   wind-independent; under a clear structured sky it is not. A uniform sky has
   L_sky = E_d/pi in every direction, so this index is 1.0 and SPECTRALLY FLAT. A
   clear sky puts most of E_d in the direct beam and leaves the 40-deg sky radiance
   well below the hemispheric average, so the index is << 1 and falls steeply with
   wavelength (Rayleigh keeps the blue sky bright, the red/NIR sky is dark). The
   spectral SLOPE is the discriminator, not the level.

2. SURFACE FRACTION,  rho*L_sky / L_t  in the visible. The leverage term. If the
   surface term is 5 % of what the sensor sees, a 20 % error in rho is a 1 % error
   in R_rs and this whole question is academic. If it is 40 %, the same rho error is
   an 8 % error and it is the dominant term in the budget.

3. rho MEASURED DIRECTLY IN THE SWIR. Past ~1000 nm pure-water absorption is
   10-1000 m^-1, penetration depth is millimetres to centimetres, and the
   water-leaving radiance is zero for ANY water, turbid included. So there

       L_t(SWIR) = rho * L_sky(SWIR)   exactly,   rho = L_t/L_sky.

   No Mobley table, no wind parameterisation, no assumption about sky state: this is
   rho under the exact conditions of that scan. The catch is SNR (both terms are
   small in the SWIR) so the scatter across bands and across scans is reported
   alongside, and the strong water-vapour bands are excluded rather than averaged
   through.

Windows avoid the 1340-1460 and 1790-1960 nm H2O absorption bands, where the
atmosphere is opaque and both radiances collapse into noise.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from analyse_location import match_by_angle
from fieldrrs.rrs import RHO_MOBLEY1999, rho_at_angle, view_zenith_from_tilt
from organize_by_location import survey

#: (path relative to by_location, label, ASOS wind m/s over the occupation)
#: wind from RHO_METHODOLOGY_REVIEW.md (Iowa Environmental Mesonet, station PAOT).
STATIONS = [
    ("LOC1_66.89718N_162.60290W/FLENS8_FOV08", "LOC1", "4.73 (3.55-5.32)"),
    ("LOC2a_66.89677N_162.57953W/FLENS8_FOV08", "LOC2a", "4.54 (4.14-4.73)"),
    ("LOC2b_66.89677N_162.57953W/FLENS8_FOV08", "LOC2b", "4.54 (4.14-4.73)"),
    ("LOC3_66.89235N_162.59149W/FIBR15_FOV15", "LOC3-FIBR15", "5.62 (5.32-5.92)"),
    ("LOC3_66.89235N_162.59149W/FIBR15_FOV15_murky", "LOC3-FIBR15murky", "5.62 (5.32-5.92)"),
    ("LOC3_66.89235N_162.59149W/FLENS8_FOV08", "LOC3-FLENS8", "5.32"),
]

SKY_BANDS = [443, 555, 665, 750, 865]
#: (lo, hi, label) -- SWIR windows clear of the strong H2O bands
SWIR_WINDOWS = [(1000, 1060, "1.03um"), (1240, 1300, "1.27um"),
                (1600, 1700, "1.65um"), (2150, 2250, "2.20um")]


def _bm(wl, y, lo, hi):
    m = (wl >= lo) & (wl <= hi)
    if not m.any():
        return np.nan
    v = np.asarray(y, float)[m]
    return float(np.nanmedian(v))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("day_folder")
    ap.add_argument("--panel-reflectance", type=float, default=0.99)
    a = ap.parse_args()
    byloc = os.path.join(a.day_folder.rstrip("/"), "by_location")

    print("=" * 100)
    print("CAN rho BE MEASURED RATHER THAN ASSUMED?  (all numbers from scans on disk)")
    print("=" * 100)

    for rel, label, wind in STATIONS:
        folder = os.path.join(byloc, rel)
        if not os.path.isdir(folder):
            print("  ! missing, skipped: %s" % folder)
            continue
        scans = survey(folder)
        wl = np.array(scans[0]["spec"].wavelength, float)
        water = sorted([s for s in scans if s["role"] == "water"], key=lambda x: x["n"])
        sky = [s for s in scans if s["role"] == "sky"]
        if not water or not sky:
            print("  ! no water/sky at %s" % label)
            continue

        print("\n" + "-" * 100)
        print("%s   n_water=%d  n_sky=%d   ASOS wind %s m/s" % (label, len(water), len(sky), wind))
        print("-" * 100)

        # ---- 1. sky structure index, pi*L_sky/E_d  =  R_p * L_sky / L_panel
        idx = {b: [] for b in SKY_BANDS}
        for w, sk, _dm, _ in match_by_angle(water, sky):
            lp = w["spec"].columns["rad_ref"]          # panel radiance in the water file
            ls = sk["spec"].columns["rad_target"]
            for b in SKY_BANDS:
                num = _bm(wl, ls, b - 5, b + 5)
                den = _bm(wl, lp, b - 5, b + 5)
                if den and np.isfinite(den) and den > 0:
                    idx[b].append(a.panel_reflectance * num / den)
        vals = [np.nanmedian(idx[b]) if idx[b] else np.nan for b in SKY_BANDS]
        print("  SKY STRUCTURE INDEX  pi*L_sky/E_d   (uniform/overcast sky = 1.0 and FLAT)")
        print("     " + "  ".join("%4d nm %6.3f" % (b, v) for b, v in zip(SKY_BANDS, vals)))
        if np.isfinite(vals[0]) and np.isfinite(vals[-1]) and vals[-1] > 0:
            ratio = vals[0] / vals[-1]
            verdict = ("UNIFORM  -> rho nearly wind-INSENSITIVE" if ratio < 1.6 else
                      "STRUCTURED -> rho IS wind-sensitive here")
            print("     blue/NIR index ratio = %.2f   %s" % (ratio, verdict))

        # ---- 2. surface fraction rho*L_sky/L_t in the visible, at Mobley rho
        print("  SURFACE FRACTION  rho*L_sky/L_t at rho=%.3f  (the leverage on R_rs)"
              % RHO_MOBLEY1999)
        frac = {b: [] for b in (443, 555, 665)}
        for w, sk, _dm, _ in match_by_angle(water, sky):
            rho = rho_at_angle(view_zenith_from_tilt(w["spec"].tilt_y_deg))
            lt = w["spec"].columns["rad_target"]
            ls = sk["spec"].columns["rad_target"]
            for b in frac:
                t = _bm(wl, lt, b - 5, b + 5)
                s = _bm(wl, ls, b - 5, b + 5)
                if t and np.isfinite(t) and t > 0:
                    frac[b].append(rho * s / t)
        print("     " + "  ".join("%4d nm %5.1f %%" % (b, 100 * np.nanmedian(frac[b]))
                                 for b in (443, 555, 665)))

        # ---- 3. rho measured directly in the SWIR, where L_w == 0
        print("  rho MEASURED IN THE SWIR   (L_w = 0 there, so rho = L_t/L_sky)")
        for lo, hi, name in SWIR_WINDOWS:
            rr, lts, lss = [], [], []
            for w, sk, _dm, _ in match_by_angle(water, sky):
                t = _bm(wl, w["spec"].columns["rad_target"], lo, hi)
                s = _bm(wl, sk["spec"].columns["rad_target"], lo, hi)
                lts.append(t); lss.append(s)
                if s and np.isfinite(s) and s > 0 and np.isfinite(t):
                    rr.append(t / s)
            if not rr:
                print("     %-7s no usable bands" % name); continue
            rr = np.array(rr, float)
            q16, q50, q84 = np.nanpercentile(rr, [16, 50, 84])
            print("     %-7s rho = %6.3f  [16-84%%: %6.3f - %6.3f]   "
                  "L_t=%.3g  L_sky=%.3g  (n=%d)"
                  % (name, q50, q16, q84, np.nanmedian(lts), np.nanmedian(lss), len(rr)))

    print("\n" + "=" * 100)
    print("Mobley (1999) rho = 0.028 is the 40deg/135deg, wind<5 m/s, clear-sky reference.")
    print("A SWIR rho far above it means the real surface term is larger than assumed;")
    print("far below, or negative, means the SWIR is noise and this route is unavailable.")
    print("=" * 100)


if __name__ == "__main__":
    main()
