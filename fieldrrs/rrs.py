"""Above-water remote-sensing reflectance from field radiometer scans.

Pure standard library. The physics is three equations:

    E_d(l)   = pi * L_panel(l) / R_panel(l)                              (1)
    L_w(l)   = L_target(l) - rho * L_sky(l)                              (2)
    R_rs(l)  = L_w(l) / E_d(l)  -  delta                                 (3)

References
----------
Mobley, C.D. (1999), Applied Optics 38(36), 7442-7455, doi:10.1364/AO.38.007442
    the 40 deg / 135 deg geometry and rho = 0.028.
Ruddick, K.G., et al. (2006), L&O 51(2), 1167-1179, doi:10.4319/lo.2006.51.2.1167
    the NIR residual corrections.
"""

from __future__ import annotations

import math

__all__ = [
    "RHO_MOBLEY1999", "SIMILARITY_780_870", "RrsResult",
    "rrs_three_scan", "rrs_from_sed", "residual_correction", "rho_advice",
]

#: Mobley (1999) effective sea-surface radiance reflectance for the recommended
#: geometry: 40 deg from nadir, 135 deg azimuth from the sun, wind below ~5 m/s,
#: clear sky. NOT the Fresnel coefficient, and NOT valid at higher wind.
RHO_MOBLEY1999 = 0.028

#: Ruddick et al. (2006) NIR similarity ratio R_rs(780)/R_rs(870).
SIMILARITY_780_870 = 1.912

#: Typical Spectralon plaque reflectance. Use YOUR panel's calibration if you have it.
DEFAULT_PANEL_REFLECTANCE = 0.99

#: A band counts as materially over-subtracted when L_w is negative by more than this
#: fraction of L_target. Below it, the excursion is indistinguishable from a genuinely
#: zero water-leaving signal in the NIR.
OVER_SUBTRACTION_FRAC = 0.02


class RrsResult(object):
    def __init__(self, wavelength, rrs, ed, lw, rho, offset, method, notes, meta=None):
        self.wavelength = wavelength
        self.rrs = rrs
        self.ed = ed
        self.lw = lw
        self.rho = rho
        self.offset = offset
        self.method = method
        self.notes = notes
        self.meta = meta or {}

    def clipped(self, lo, hi):
        keep = [i for i, w in enumerate(self.wavelength) if lo <= w <= hi]
        return ([self.wavelength[i] for i in keep], [self.rrs[i] for i in keep])

    def value_at(self, target):
        i = min(range(len(self.wavelength)),
                key=lambda j: abs(self.wavelength[j] - target))
        return self.wavelength[i], self.rrs[i]


def rho_advice(wind_ms):
    """What to do about rho at a given wind speed. Returns (rho_or_None, message).

    There is no wind-dependent rho table bundled here. Mobley (2015) published one and
    it is not redistributable, so above the validity range of the 0.028 value this
    function refuses to invent a number and says so. Using 0.028 in a stiff breeze
    biases blue R_rs high by tens of percent, and that error is invisible in the output.
    """
    if wind_ms is None:
        return RHO_MOBLEY1999, ("No wind speed recorded. Using rho = 0.028, which "
                                "assumes wind below ~5 m/s and clear sky. Record the "
                                "wind, it is the largest error in this measurement.")
    if wind_ms <= 5.0:
        return RHO_MOBLEY1999, ("rho = 0.028 (Mobley 1999) is appropriate for "
                                "%.1f m/s at 40 deg / 135 deg, clear sky." % wind_ms)
    return None, (
        "Wind is %.1f m/s, above the ~5 m/s validity limit of rho = 0.028. The correct "
        "value is larger and comes from the Mobley (2015) wind- and geometry-dependent "
        "lookup table, which is NOT bundled with this software. Options: (a) proceed "
        "with 0.028 and record that blue R_rs is biased HIGH, (b) supply a rho from the "
        "published table, (c) re-measure in calmer conditions." % wind_ms)


def rrs_three_scan(wavelength, l_target, l_sky, l_panel,
                   panel_reflectance=DEFAULT_PANEL_REFLECTANCE,
                   rho=RHO_MOBLEY1999, residual="none",
                   residual_window=(750.0, 800.0)):
    """Equations (1)-(3). All radiance lists must share one wavelength grid."""
    n = len(wavelength)
    for name, arr in (("l_target", l_target), ("l_sky", l_sky), ("l_panel", l_panel)):
        if len(arr) != n:
            raise ValueError(
                "%s has %d points but the wavelength grid has %d. All three scans must "
                "come from the same instrument on the same grid." % (name, len(arr), n))

    if isinstance(panel_reflectance, (int, float)):
        rp = [float(panel_reflectance)] * n
    else:
        rp = list(panel_reflectance)
        if len(rp) != n:
            raise ValueError("panel_reflectance spectrum length does not match the grid")

    notes = []
    ed, lw, rrs = [], [], []
    n_bad = 0
    n_over = 0
    for i in range(n):
        e = math.pi * l_panel[i] / rp[i] if rp[i] > 0 else float("nan")
        w = l_target[i] - rho * l_sky[i]
        ed.append(e)
        lw.append(w)
        if e and e == e and e != 0.0:
            rrs.append(w / e)
        else:
            rrs.append(float("nan"))
            n_bad += 1
        # MATERIAL over-subtraction only. In the far NIR the true water-leaving signal is
        # essentially zero, so L_w hovers at +/- rounding and a bare `w < 0` test fires on
        # a perfectly good measurement. A warning that cries wolf on correct data is worse
        # than no warning, so this requires the negative excursion to be a real fraction of
        # what was measured.
        if w < 0 and abs(w) > OVER_SUBTRACTION_FRAC * abs(l_target[i]):
            n_over += 1

    if n_bad:
        notes.append("%d bands have zero or invalid panel radiance and are NaN." % n_bad)
    if n_over:
        frac = 100.0 * n_over / n
        notes.append(
            "%d bands (%.0f%%) have L_target below rho*L_sky by more than %.0f%% of the "
            "measured radiance: the glint subtraction is removing materially more than "
            "was measured. Either rho is too large for this sea state, or the sky scan "
            "does not mirror the water view (check the bearing)."
            % (n_over, frac, 100 * OVER_SUBTRACTION_FRAC))

    offset, method, more = residual_correction(wavelength, rrs, residual, residual_window)
    notes.extend(more)
    rrs = [v - offset for v in rrs]

    neg_vis = sum(1 for w, v in zip(wavelength, rrs)
                  if 400 <= w <= 700 and v == v and v < 0)
    if neg_vis:
        notes.append(
            "%d visible bands are NEGATIVE after correction. R_rs cannot be negative; "
            "this indicates over-subtraction (rho too high, or the NIR correction is "
            "invalid for this water)." % neg_vis)

    return RrsResult(list(wavelength), rrs, ed, lw, float(rho), float(offset),
                     method, notes)


def residual_correction(wavelength, rrs, method="none", window=(750.0, 800.0)):
    """Residual glint / offset correction. Returns (offset, method, notes).

    'none'           subtract nothing. The default, deliberately.
    'nir_zero'       assume R_rs = 0 over `window` and subtract its mean.
                     VALID in clear oceanic water. INVALID in turbid water, where it
                     deletes real signal.
    'nir_similarity' Ruddick et al. (2006), uses the fixed R_rs(780)/R_rs(870) ratio.
                     Usable in turbid water.
    """
    notes = []
    if method in (None, "none", ""):
        return 0.0, "none", notes

    if method == "nir_zero":
        vals = [v for w, v in zip(wavelength, rrs)
                if window[0] <= w <= window[1] and v == v]
        if not vals:
            raise ValueError(
                "nir_zero needs bands in %.0f-%.0f nm; this spectrum spans %.0f-%.0f nm"
                % (window[0], window[1], wavelength[0], wavelength[-1]))
        off = sum(vals) / len(vals)
        notes.append(
            "nir_zero subtracted %.3e sr^-1 (mean of %.0f-%.0f nm, %d bands). This "
            "assumes zero water-leaving signal in the NIR: INVALID in turbid or "
            "sediment-laden water." % (off, window[0], window[1], len(vals)))
        return off, "nir_zero", notes

    if method == "nir_similarity":
        i780 = min(range(len(wavelength)), key=lambda j: abs(wavelength[j] - 780.0))
        i870 = min(range(len(wavelength)), key=lambda j: abs(wavelength[j] - 870.0))
        if abs(wavelength[i780] - 780) > 10 or abs(wavelength[i870] - 870) > 10:
            raise ValueError(
                "nir_similarity needs bands near 780 and 870 nm; nearest are %.1f and "
                "%.1f nm" % (wavelength[i780], wavelength[i870]))
        a = SIMILARITY_780_870
        off = (a * rrs[i870] - rrs[i780]) / (a - 1.0)
        notes.append(
            "nir_similarity subtracted %.3e sr^-1 using the Ruddick et al. (2006) "
            "ratio %.3f at %.0f/%.0f nm." % (off, a, wavelength[i780], wavelength[i870]))
        return off, "nir_similarity", notes

    raise ValueError("unknown residual method %r" % (method,))


def rrs_from_sed(water, sky, panel=None,
                 panel_reflectance=DEFAULT_PANEL_REFLECTANCE,
                 rho=RHO_MOBLEY1999, residual="none", source="radiance"):
    """Build R_rs from DARWin .sed scans.

    TWO WORKFLOWS, both supported:

    * **Two files (the DARWin-native way).** Scan the panel as the REFERENCE, then take
      the sky and the water as TARGETS. Both files then carry the panel in their
      ``Rad. (Ref.)`` column, so ``panel=None`` and the panel radiance is read from the
      water file's own reference column. Fewer files, and the panel is guaranteed
      contemporaneous with the target.
    * **Three files.** A separate panel scan passed as ``panel``.

    ``source='radiance'`` uses the calibrated radiance columns and is preferred.
    ``source='reflectance'`` falls back to the ``Reflect.`` ratio columns, valid only
    when every scan shares one reference panel, since the panel radiance then cancels.
    """
    wl = water.wavelength
    if list(sky.wavelength) != list(wl):
        raise ValueError(
            "water (%s) and sky (%s) are on different wavelength grids. They must come "
            "from the same instrument configuration." % (water.name, sky.name))

    meta = {"water_file": water.name, "sky_file": sky.name,
            "panel_file": panel.name if panel is not None else
                          "%s [Rad. (Ref.) column]" % water.name,
            "water_time": water.when, "comment": water.comment,
            "lat": water.latitude, "lon": water.longitude}

    if source == "radiance":
        lt = water.radiance_target
        ls = sky.radiance_target
        if panel is not None:
            if list(panel.wavelength) != list(wl):
                raise ValueError("panel scan is on a different wavelength grid")
            lp = panel.radiance_target
        else:
            lp = water.radiance_reference
        res = rrs_three_scan(wl, lt, ls, lp, panel_reflectance, rho, residual)
        res.meta = meta
        return res

    if source == "reflectance":
        rt = water.reflectance
        rs = sky.reflectance
        rp = (panel_reflectance if isinstance(panel_reflectance, (int, float))
              else panel_reflectance)
        scale = (rp if isinstance(rp, (int, float)) else 1.0) / math.pi
        rrs = [(rt[i] - rho * rs[i]) * scale for i in range(len(wl))]
        off, method, notes = residual_correction(wl, rrs, residual)
        rrs = [v - off for v in rrs]
        notes.insert(0, "reflectance path: assumes water, sky and panel scans all used "
                        "the SAME reference panel, so the panel radiance cancels.")
        res = RrsResult(list(wl), rrs, None, None, float(rho), float(off), method, notes)
        res.meta = meta
        return res

    raise ValueError("source must be 'radiance' or 'reflectance'")


def average_results(results):
    """Mean and standard deviation across replicate R_rs spectra.

    Replicates are the only uncertainty estimate available in the field, and they cover
    only the random part: they say nothing about rho or the panel calibration, which are
    systematic and shared by every replicate.
    """
    if not results:
        raise ValueError("no results to average")
    wl = results[0].wavelength
    for r in results[1:]:
        if list(r.wavelength) != list(wl):
            raise ValueError("replicates are on different wavelength grids")
    n = len(results)
    mean, sd = [], []
    for i in range(len(wl)):
        vals = [r.rrs[i] for r in results if r.rrs[i] == r.rrs[i]]
        if not vals:
            mean.append(float("nan")); sd.append(float("nan")); continue
        m = sum(vals) / len(vals)
        mean.append(m)
        sd.append(math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1))
                  if len(vals) > 1 else 0.0)
    out = RrsResult(list(wl), mean, None, None, results[0].rho, 0.0,
                    "mean of %d replicates" % n,
                    ["standard deviation is RANDOM error only; rho and panel "
                     "calibration are systematic and are not captured here."])
    out.sd = sd
    return out
