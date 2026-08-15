"""Spectral resampling and export. Pure standard library.

Resample R_rs, never a retrieved IOP: the ocean-colour forward model is nonlinear in the
IOPs, so band-averaging does not commute with an inversion.
"""

from __future__ import annotations

import csv
import math

__all__ = ["interp_linear", "bin_spectrum", "gaussian_resample",
           "SATELLITE_BANDS", "write_rrs_csv", "write_batch_csv"]

SATELLITE_BANDS = {
    "PACE OCI (GIOP set)": [412.0, 443, 490, 510, 555, 670],
    "MODIS-Aqua": [412.0, 443, 488, 531, 547, 667],
    "SeaWiFS": [412.0, 443, 490, 510, 555, 670],
    "VIIRS": [410.0, 443, 486, 551, 671],
    "Sentinel-3 OLCI": [400.0, 412.5, 442.5, 490, 510, 560, 620, 665, 674, 681.25, 709],
    "1 nm, 400-700": [float(w) for w in range(400, 701)],
    "5 nm, 400-900": [float(w) for w in range(400, 901, 5)],
}


def interp_linear(x, y, xq):
    """Linear interpolation, NaN outside the data range (MATLAB/numpy convention)."""
    out = []
    n = len(x)
    for q in xq:
        if q < x[0] or q > x[-1]:
            out.append(float("nan"))
            continue
        lo, hi = 0, n - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if x[mid] <= q:
                lo = mid
            else:
                hi = mid
        if x[hi] == x[lo]:
            out.append(y[lo])
        else:
            t = (q - x[lo]) / (x[hi] - x[lo])
            out.append(y[lo] * (1 - t) + y[hi] * t)
    return out


def bin_spectrum(wl, values, centres, width=10.0, min_bands=1):
    """Boxcar average onto `centres` over a full width of `width` nm.

    Bins with too few contributing channels come back NaN rather than being quietly
    extrapolated. Returns (centres, binned, n_per_bin).
    """
    out, counts = [], []
    half = width / 2.0
    for c in centres:
        vals = [v for w, v in zip(wl, values)
                if c - half <= w <= c + half and v == v]
        counts.append(len(vals))
        out.append(sum(vals) / len(vals) if len(vals) >= min_bands and vals
                   else float("nan"))
    return list(centres), out, counts


def gaussian_resample(wl, values, centres, fwhm=10.0):
    """Convolve with a Gaussian of the given FWHM at each centre.

    A centre whose response is more than half outside the measured range returns NaN,
    measured against the full Gaussian rather than against the truncated sum.
    """
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    lo, hi = wl[0], wl[-1]
    root2 = math.sqrt(2.0)
    out = []
    for c in centres:
        coverage = 0.5 * (math.erf((hi - c) / (sigma * root2))
                          - math.erf((lo - c) / (sigma * root2)))
        if coverage <= 0.5:
            out.append(float("nan"))
            continue
        num = den = 0.0
        for w, v in zip(wl, values):
            if v != v:
                continue
            g = math.exp(-0.5 * ((w - c) / sigma) ** 2)
            num += g * v
            den += g
        out.append(num / den if den > 0 else float("nan"))
    return list(centres), out


def write_rrs_csv(path, result, extra_columns=None):
    """One spectrum per file: wavelength, R_rs, and any diagnostics available."""
    cols = [("wavelength_nm", result.wavelength), ("Rrs_sr-1", result.rrs)]
    if getattr(result, "sd", None) is not None:
        cols.append(("Rrs_sd_sr-1", result.sd))
    if result.ed is not None:
        cols.append(("Ed_W_m-2_nm-1", result.ed))
    if result.lw is not None:
        cols.append(("Lw_W_m-2_sr-1_nm-1", result.lw))
    for name, vals in (extra_columns or []):
        cols.append((name, vals))

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# fieldrrs above-water Rrs"])
        w.writerow(["# rho", result.rho, "residual_method", result.method,
                    "offset_sr-1", result.offset])
        for k, v in sorted(result.meta.items()):
            w.writerow(["#", k, v])
        for note in result.notes:
            w.writerow(["# WARNING", note])
        w.writerow([c[0] for c in cols])
        for i in range(len(result.wavelength)):
            w.writerow(["%.6g" % c[1][i] for c in cols])
    return path


def write_batch_csv(path, named_results, centres=None, width=10.0):
    """Many stations in one wide table: one row per station, one column per band.

    This is the file to hand to an inversion (GIOP, titanspec) or to a spreadsheet.
    """
    named_results = list(named_results)
    if not named_results:
        raise ValueError("nothing to write")
    if centres is None:
        centres = named_results[0][1].wavelength

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["station"] + ["%.1f" % c for c in centres])
        for name, res in named_results:
            _, binned, _ = bin_spectrum(res.wavelength, res.rrs, centres, width)
            w.writerow([name] + ["%.6g" % v for v in binned])
    return path
