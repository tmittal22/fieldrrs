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
    "scaled_mean",
    "refractive_index_water", "fresnel_reflectance", "rho_at_angle",
    "view_zenith_from_tilt",
    "RHO_MOBLEY1999", "SIMILARITY_780_870", "RrsResult",
    "rrs_three_scan", "rrs_from_sed", "residual_correction", "rho_advice",
    "overcast_notes", "par_from_ed", "integrated_irradiance", "ed_stability",
    "cross_calibration_factor", "rrs_from_separate_ed",
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



# ---------------------------------------------------------------------------
# Angle-dependent rho
# ---------------------------------------------------------------------------

def refractive_index_water(wavelength_nm, temperature_c=10.0, salinity=30.0):
    """Real refractive index of seawater. Quan & Fry (1995) doi:10.1364/AO.34.003477.

    Valid 400-700 nm, 0-30 degC, 0-35 PSU. Defaults are Arctic coastal values.
    """
    wl = float(wavelength_nm)
    t, sal = float(temperature_c), float(salinity)
    n0, n1, n2, n3, n4 = 1.31405, 1.779e-4, -1.05e-6, 1.6e-8, -2.02e-6
    n5, n6, n7, n8, n9 = 15.868, 0.01155, -0.00423, -4382.0, 1.1455e6
    return (n0 + (n1 + n2 * t + n3 * t * t) * sal + n4 * t * t
            + (n5 + n6 * sal + n7 * t) / wl + n8 / (wl * wl) + n9 / (wl ** 3))


def fresnel_reflectance(theta_deg, n=1.34):
    """Unpolarised Fresnel reflectance of a FLAT air-water interface.

    theta_deg is measured from the surface normal; for the sky-glint ray it equals the
    view zenith angle.
    """
    ti = math.radians(float(theta_deg))
    st = math.sin(ti) / n
    if st >= 1.0:
        return 1.0
    tt = math.asin(st)
    ci, ct = math.cos(ti), math.cos(tt)
    rs = ((ci - n * ct) / (ci + n * ct)) ** 2
    rp = ((n * ci - ct) / (n * ci + ct)) ** 2
    return 0.5 * (rs + rp)


def rho_at_angle(view_zenith_deg, rho_ref=RHO_MOBLEY1999, theta_ref=40.0,
                 wavelength_nm=550.0, temperature_c=10.0, salinity=30.0):
    """rho scaled from its reference angle to the angle you actually achieved:

        rho(theta) = rho_ref * R_Fresnel(theta) / R_Fresnel(theta_ref)

    WHY A RATIO. rho = 0.028 is Mobley's (1999) EFFECTIVE value at 40 deg from nadir,
    135 deg azimuth, wind under ~5 m/s. Its magnitude comes from averaging over wave
    facets and is not Fresnel. Its ANGULAR DEPENDENCE does follow the flat-surface
    Fresnel curve closely for modest departures, because tilting the view changes the
    mean incidence angle on the facets much as it does on a flat surface. So the ratio
    is the defensible part; the absolute level is inherited from Mobley.

    It is worth doing. Over 42-50 deg -- an ordinary hand-aimed spread -- Fresnel
    reflectance rises about 32 %, and a fixed rho puts all of that into R_rs as a
    systematic tied to your pointing. Measured on LOC1 (n=12): correcting each scan at
    its own angle removed a p=0.006 trend and cut scan-to-scan scatter from 11.1 % to
    8.6 %.

    LIMITS. First order: it captures the angular trend only, not the wind or sky-state
    dependence, and should not be pushed far from theta_ref. Beyond roughly +/-15 deg,
    or above ~5 m/s wind, the real Mobley (2015) tables are needed and those are not
    redistributable. Assumes the reference azimuth is maintained.
    """
    n = refractive_index_water(wavelength_nm, temperature_c, salinity)
    return rho_ref * fresnel_reflectance(view_zenith_deg, n) / \
        fresnel_reflectance(theta_ref, n)


def view_zenith_from_tilt(tilt_y_deg):
    """View zenith angle from the NaturaSpec's reported |tilt|.

    The instrument reports a magnitude with no stated datum, and the two candidate
    readings differ by 90 deg, so this was settled AGAINST THE DATA. On LOC1 (12 water
    scans, 39.5-49.6 deg), re-deriving each scan with :func:`rho_at_angle` at its own
    angle gives, for the tilt-vs-R_rs(443) trend:

        fixed rho = 0.028      r=+0.74 (p=0.006), scatter 11.1 %
        theta_v = 90 - tilt    r=+0.90 (p<0.001), scatter 18.9 %   WORSE
        theta_v = tilt         r=+0.36 (p=0.244), scatter  8.6 %   trend gone

    So `tilt_y` IS the view zenith angle. That reading removes the systematic AND
    tightens the scatter by 23 %; a wrong correction cannot do both, it adds variance,
    as the other reading does. The values also bracket the nominal 40 deg from nadir,
    which is what an operator aiming for 40 deg produces.

    Evidence, not proof: n=12 at one location. Re-check on LOC2/LOC3. Under this reading
    the sky scans sit 44-50 deg from ZENITH.
    """
    return abs(float(tilt_y_deg))


def scaled_mean(spectra, iterations=25, tol=1e-12, weight=None):
    """Iterative amplitude-normalised mean: a clean SHAPE plus per-scan amplitudes.

    When replicate spectra differ mostly by a multiplicative factor -- which is what a
    changing amount of the same scattering material does -- a plain average is noisier
    than it needs to be, because the amplitude scatter leaks into every band. This
    alternates two steps to convergence:

        a_i  = <R_i, M> / <R_i, R_i>       best scale taking R_i onto the current mean
        M    = mean_i(a_i R_i)             re-mean the rescaled spectra

    then rescales M so the amplitudes have mean 1, which keeps the result in the
    original physical units instead of drifting.

    Returns (mean, scales, n_iter, shape_cv, amp_cv):
      mean      the shape at the mean amplitude, in sr^-1
      scales    a_i per input; 1/a_i is that scan's amplitude relative to the mean
      shape_cv  residual scatter AFTER rescaling, in per cent -- the shape uncertainty
      amp_cv    scatter of the amplitudes, in per cent -- REAL variability, not error

    Separating those two is the point. Reporting one combined error bar implies the
    water type is uncertain when only its concentration is.

    ONLY valid when the variation really is multiplicative. Check `shape_cv` against the
    unscaled spread: if it barely falls, the spectra differ in shape and this estimator
    is hiding that rather than helping. Pure standard library.
    """
    X = [list(map(float, r)) for r in spectra]
    n, m = len(X), len(X[0])
    if n == 0:
        raise ValueError("no spectra")
    w = [1.0] * m if weight is None else list(map(float, weight))
    if len(w) != m:
        raise ValueError("weight must match the spectrum length")

    def wdot(u, v):
        return sum(wi * ui * vi for wi, ui, vi in zip(w, u, v))

    mean = [sum(x[j] for x in X) / n for j in range(m)]
    scales = [1.0] * n
    used = 0
    for used in range(1, iterations + 1):
        num = [wdot(x, mean) for x in X]
        den = [wdot(x, x) for x in X]
        scales = [(nu / de if de > 0 else 1.0) for nu, de in zip(num, den)]
        g = sum(scales) / n
        if g <= 0:
            break
        scales = [a / g for a in scales]          # keep the physical scale
        new = [sum(a * x[j] for a, x in zip(scales, X)) / n for j in range(m)]
        shift = max(abs(a - b) for a, b in zip(new, mean)) / max(
            1e-30, max(abs(v) for v in mean))
        mean = new
        if shift < tol:
            break

    def cv(vals):
        mu = sum(vals) / len(vals)
        if mu == 0:
            return float("nan")
        var = sum((v - mu) ** 2 for v in vals) / len(vals)
        return 100.0 * math.sqrt(var) / abs(mu)

    res = []
    for j in range(m):
        col = [a * x[j] for a, x in zip(scales, X)]
        if mean[j] != 0 and w[j] > 0:
            res.append(cv(col))
    shape_cv = sum(res) / len(res) if res else float("nan")
    amp_cv = cv([1.0 / a if a else float("nan") for a in scales])
    return mean, scales, used, shape_cv, amp_cv


def rho_advice(wind_ms, sky="clear"):
    """What to do about rho at a given wind speed and sky. Returns (rho_or_None, msg).

    ``sky='overcast'`` changes the answer, and it changes it in the helpful direction.
    The wind sensitivity of rho comes from wave facets sampling DIFFERENT PARTS OF A
    NON-UNIFORM SKY: tilt a facet under a clear sky and it swings between bright
    horizon and dark zenith, so the effective reflectance depends on the facet
    distribution and therefore on wind. Under a uniformly overcast sky there is little
    to sample between, so facet orientation stops mattering and rho becomes close to the
    flat-surface Fresnel value and largely wind-insensitive.

    There is still no bundled wind-dependent table (Mobley 2015 is not redistributable),
    so under a CLEAR sky above ~5 m/s this refuses to invent a number.
    """
    if sky == "overcast":
        return RHO_MOBLEY1999, (
            "Uniform overcast: rho is close to the flat-surface Fresnel value and only "
            "weakly wind-dependent, because a uniform sky gives wave facets little to "
            "sample between. rho = 0.028 is a reasonable choice at %s and the wind "
            "caveat is much weaker than under clear sky. THE REAL RISK IS BROKEN CLOUD, "
            "not thick cloud: if the sky is patchy or the cover is changing, E_d moves "
            "between scans and neither rho nor the panel reference is stable."
            % ("%.1f m/s" % wind_ms if wind_ms is not None else "any wind"))
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


def overcast_notes():
    """What changes when the sky is fully overcast. Returned as operator-facing text."""
    return [
        "NO SUN GLINT. The dominant error in above-water radiometry is the specular "
        "sun beam, and under full cloud there is no beam. This is the one respect in "
        "which overcast is BETTER than clear sky.",
        "The 135 deg relative azimuth stops meaning anything: there is no sun direction "
        "to point away from. Keep the 40 deg view angle, and choose the bearing purely "
        "to avoid the platform, its shadow and its wake.",
        "rho is more stable and only weakly wind-dependent, because a uniform sky gives "
        "wave facets little to sample between.",
        "E_d is much lower, so increase integration time and take more replicates. The "
        "signal is real but the signal-to-noise is not what it was.",
        "DO NOT apply the BRDF / f-over-Q normalisation. The Morel tables describe a "
        "clear-sky light field with a direct beam; under full overcast the field is "
        "entirely diffuse and those tables do not describe it.",
        "No satellite match-up is possible: an optical sensor cannot see the water "
        "through the cloud either.",
        "BROKEN CLOUD IS THE WORST CASE, worse than either clear or fully overcast, "
        "because E_d changes between the panel scan and the target scan. If you have a "
        "cosine-collector irradiance channel, use source='irradiance': measuring E_d "
        "simultaneously removes that time lag entirely.",
    ]


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

    ``source='irradiance'`` uses a MEASURED E_d from a cosine collector instead of
    inferring it from a panel. Preferred whenever the light is changing, which is most
    of what makes cloudy-day work hard.

    ``source='radiance'`` uses the calibrated radiance columns and is preferred.
    ``source='reflectance'`` falls back to the ``Reflect.`` ratio columns, valid only
    when every scan shares one reference panel, since the panel radiance then cancels.
    """
    # The sky scan is the ONE mandatory extra input, so it gets an explicit refusal.
    # Without it, rho*L_sky is never subtracted and R_rs comes out far too high in the
    # blue -- a wrong answer that looks like a plausible spectrum. The GUI blocks this
    # earlier with its own message; this is for anyone driving the library directly, who
    # otherwise got AttributeError on NoneType and no idea why.
    if water is None:
        raise ValueError("rrs_from_sed needs a WATER scan.")
    if sky is None:
        raise ValueError(
            "rrs_from_sed needs a SKY scan, and it is not optional. Without it the "
            "reflected skylight rho*L_sky stays in the signal and R_rs comes out far "
            "too high in the blue. It is the 40-deg-from-ZENITH scan on the SAME "
            "compass bearing as the water scan. (The PANEL, by contrast, IS optional: "
            "omit it and the panel radiance is read from the water file's own "
            "'Rad. (Ref.)' column.)")

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
        elif water.has("rad_ref"):
            lp = water.radiance_reference
        else:
            # Water + sky with no panel anywhere. This is the commonest field mistake
            # and it is NOT recoverable: the sky scan removes reflected skylight, it
            # does not tell you how bright the light was. Raise something that names
            # all three cures instead of a KeyError on 'rad_ref'.
            raise ValueError(
                "No irradiance reference. R_rs = L_water / E_d, and nothing here "
                "supplies E_d: no panel scan was passed, and '%s' has no "
                "'Rad. (Ref.)' column (it has: %s). The SKY scan does not supply it. "
                "Supply ONE of: (1) a panel scan as `panel=`; (2) files re-exported "
                "from DARWin in REFLECTANCE mode, which write the panel into every "
                "file's 'Rad. (Ref.)' column; (3) a measured irradiance spectrum via "
                "`rrs_from_separate_ed(...)`."
                % (water.name, ", ".join(sorted(water.raw_columns))))
        res = rrs_three_scan(wl, lt, ls, lp, panel_reflectance, rho, residual)
        res.meta = meta
        return res

    if source == "irradiance":
        # E_d measured DIRECTLY by a cosine collector, instead of inferred from a panel.
        # Strictly better whenever the light is changing, and that is most of what makes
        # cloudy-day radiometry hard: the panel method assumes E_d is the same at the
        # moment of the panel scan and the moment of the target scan, and under moving
        # cloud it is not. A simultaneous irradiance channel removes the time lag, and
        # with it the panel reflectance (A11) and the panel-levelness error as well.
        lt = water.radiance_target
        ls = sky.radiance_target
        ed = None
        for src in (water, panel):
            if src is not None and src.has("irr_target"):
                ed = src.columns["irr_target"]
                break
            if src is not None and src.has("irr_ref"):
                ed = src.columns["irr_ref"]
                break
        if ed is None:
            raise KeyError(
                "source='irradiance' needs an 'Irr. (Target)' or 'Irr. (Ref.)' column "
                "from a cosine-collector scan. %s has: %s. Re-export from DARWin with "
                "the irradiance channel enabled, or take the E_d scan with the cosine "
                "diffuser fitted." % (water.name, sorted(water.raw_columns)))

        notes = []
        rrs = []
        n_bad = 0
        for i in range(len(wl)):
            w = lt[i] - rho * ls[i]
            if ed[i] > 0:
                rrs.append(w / ed[i])
            else:
                rrs.append(float("nan"))
                n_bad += 1
        if n_bad:
            notes.append("%d bands have zero or negative measured E_d." % n_bad)
        notes.append("E_d measured directly by a cosine collector, not inferred from a "
                     "panel. The panel reflectance and its levelness do not enter, and "
                     "there is no panel-to-target time lag.")
        off, method, more = residual_correction(wl, rrs, residual)
        notes.extend(more)
        rrs = [v - off for v in rrs]
        res = RrsResult(list(wl), rrs, list(ed), None, float(rho), float(off), method,
                        notes)
        res.meta = dict(meta, ed_source="measured irradiance (cosine collector)")
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

    raise ValueError("source must be 'radiance', 'irradiance' or 'reflectance'; "
                     "got %r" % (source,))


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


# ----------------------------------------------------------------------------------
# Products that need an ABSOLUTE radiometer, not just a reflectance instrument.

#: Energy of one mole of photons at wavelength lambda [nm] is 1.196e8 / lambda J/mol,
#: from N_A h c = 0.1196 J m /mol. Dividing an irradiance in W m^-2 nm^-1 by that and
#: scaling to micromoles gives the constant below.
_PHOTON_CONST = 119.6


def par_from_ed(wavelength, ed, lo=400.0, hi=700.0):
    """Photosynthetically Available Radiation from a measured E_d spectrum.

    PAR is a PHOTON flux, not an energy flux, because photosynthesis counts photons:

        PAR [umol photons m^-2 s^-1] = integral(400..700) E_d(l) * l / 119.6  dl

    with E_d in W m^-2 nm^-1 and l in nm. The constant is N_A h c = 0.1196 J m /mol
    expressed for nanometres and micromoles; nothing here is fitted or tabulated.

    This is a product a reflectance-only instrument cannot give you. It needs an
    absolutely calibrated irradiance channel, which is what a spectroradiometer with a
    cosine collector is. Full midday sunlight is around 2000 umol m^-2 s^-1; heavy
    overcast is one to two orders of magnitude lower, and measuring that drop is a
    direct record of how much light the water column actually received.

    Returns (PAR, n_bands_used). Trapezoidal integration on the instrument's own grid.
    """
    pts = [(w, e) for w, e in zip(wavelength, ed)
           if lo <= w <= hi and e == e and e >= 0.0]
    if len(pts) < 2:
        raise ValueError(
            "PAR needs at least two finite E_d samples in %.0f-%.0f nm; got %d. Was the "
            "irradiance channel recorded?" % (lo, hi, len(pts)))
    total = 0.0
    for (w0, e0), (w1, e1) in zip(pts[:-1], pts[1:]):
        f0 = e0 * w0 / _PHOTON_CONST
        f1 = e1 * w1 / _PHOTON_CONST
        total += 0.5 * (f0 + f1) * (w1 - w0)
    return total, len(pts)


def integrated_irradiance(wavelength, ed, lo=400.0, hi=700.0):
    """Broadband E_d over a window, W m^-2. The energy counterpart of PAR."""
    pts = [(w, e) for w, e in zip(wavelength, ed)
           if lo <= w <= hi and e == e and e >= 0.0]
    if len(pts) < 2:
        raise ValueError("need at least two finite E_d samples in %.0f-%.0f nm"
                         % (lo, hi))
    return sum(0.5 * (e0 + e1) * (w1 - w0)
               for (w0, e0), (w1, e1) in zip(pts[:-1], pts[1:]))


def ed_stability(ed_before, ed_after, wavelength=None, lo=400.0, hi=700.0,
                 tolerance=0.02):
    """Did the light change during the station? The QC a panel alone cannot give you.

    Take an E_d scan before the target sequence and another after. If they differ, the
    station was measured under changing illumination and every R_rs in it is suspect,
    because the whole method assumes E_d is the same for the reference and the target.

    With a panel you can only test this by re-scanning the panel, which costs a scan and
    still leaves the moment between them unsampled. With a cosine collector it is two
    cheap scans that bracket the station.

    Returns a dict with the mean and worst fractional change and a verdict. ``tolerance``
    is the fraction above which the station is flagged; 2 % is a reasonable field bar.
    """
    if len(ed_before) != len(ed_after):
        raise ValueError("the two E_d scans must be on the same wavelength grid")
    idx = range(len(ed_before))
    if wavelength is not None:
        if len(wavelength) != len(ed_before):
            raise ValueError("wavelength length does not match the E_d scans")
        idx = [i for i in idx if lo <= wavelength[i] <= hi]
    pairs = [(ed_before[i], ed_after[i]) for i in idx
             if ed_before[i] == ed_before[i] and ed_before[i] > 0
             and ed_after[i] == ed_after[i]]
    if not pairs:
        raise ValueError("no usable bands in %.0f-%.0f nm" % (lo, hi))

    frac = [abs(b - a) / a for a, b in pairs]
    mean_change = sum(frac) / len(frac)
    worst = max(frac)
    ratio = sum(b for _, b in pairs) / sum(a for a, _ in pairs)
    stable = worst <= tolerance
    return {
        "mean_change": mean_change,
        "worst_change": worst,
        "after_over_before": ratio,
        "n_bands": len(pairs),
        "stable": stable,
        "verdict": (
            "E_d stable to %.1f %% across the station." % (100 * worst) if stable else
            "E_d CHANGED by up to %.1f %% (mean %.1f %%, net %+.1f %%) during the "
            "station. The method assumes E_d is the same for reference and target, so "
            "every R_rs here is suspect. Re-measure, or flag these data."
            % (100 * worst, 100 * mean_change, 100 * (ratio - 1))),
    }


# ----------------------------------------------------------------------------------
# TWO-INSTRUMENT SETUP: a separate hemispherical irradiance sensor alongside the
# narrow-FOV radiance spectroradiometer.
#
# This is the classic above-water configuration (the arrangement a HyperSAS-style rig
# uses) and it is better than a single instrument in the one way that matters most:
# E_d is logged SIMULTANEOUSLY with the radiance, so the assumption that E_d is
# unchanged between reference and target is no longer an assumption.
#
# It buys that at the cost of a new error that a single-instrument setup does not have.
# R_rs = L_w / E_d now divides a radiance measured by instrument A by an irradiance
# measured by instrument B, so ANY offset between their absolute calibrations is a
# direct multiplicative bias on every R_rs. Two separately calibrated instruments agree
# to their combined calibration uncertainty, typically several percent, and that does
# not cancel anywhere.
#
# The fix is that you already own the transfer standard: the panel. Point the radiance
# sensor at the panel while the irradiance sensor sees the same sky, and the ratio of
# the two E_d estimates is the inter-instrument factor.


def cross_calibration_factor(wavelength, l_panel, ed_measured, ed_wavelength=None,
                             panel_reflectance=DEFAULT_PANEL_REFLECTANCE):
    """Tie a separate irradiance sensor to the radiance sensor, using the panel.

    Two independent estimates of the same downwelling irradiance:

        E_d from the panel  = pi * L_panel / R_panel     (radiance instrument)
        E_d measured        = the irradiance sensor's own reading

    Their ratio is the inter-instrument calibration factor C(lambda):

        C(lambda) = [pi * L_panel(lambda) / R_panel(lambda)] / E_d_measured(lambda)

    **Use it as a diagnostic first.** If C is near 1 and spectrally flat, the two
    instruments agree and the measured E_d can be trusted directly. If C is far from 1,
    or has structure in wavelength, something is wrong: a stale calibration, a tilted
    collector, a shaded panel, or a mismatch in how the two were calibrated.

    You *can* multiply by C to force agreement, but understand what that costs: it makes
    the result depend on the panel reflectance again, which is exactly the error the
    separate irradiance sensor was removing. Correct only if you have reason to trust the
    panel more than the irradiance calibration.

    ``ed_wavelength`` lets the irradiance sensor live on its own grid; E_d is
    interpolated onto the radiance grid, which is what a two-instrument setup needs.
    """
    from .resample import interp_linear

    n = len(wavelength)
    if len(l_panel) != n:
        raise ValueError("l_panel does not match the radiance wavelength grid")

    if ed_wavelength is not None:
        ed = interp_linear(list(ed_wavelength), list(ed_measured), list(wavelength))
    else:
        if len(ed_measured) != n:
            raise ValueError(
                "ed_measured has %d points against %d radiance bands. Pass "
                "ed_wavelength so the two grids can be matched: separate instruments "
                "rarely share a grid." % (len(ed_measured), n))
        ed = list(ed_measured)

    rp = ([float(panel_reflectance)] * n
          if isinstance(panel_reflectance, (int, float)) else list(panel_reflectance))

    c, used = [], []
    for i in range(n):
        e_panel = math.pi * l_panel[i] / rp[i] if rp[i] > 0 else float("nan")
        if ed[i] and ed[i] > 0 and e_panel == e_panel:
            c.append(e_panel / ed[i])
            used.append(wavelength[i])
        else:
            c.append(float("nan"))

    vis = [c[i] for i in range(n)
           if 400.0 <= wavelength[i] <= 700.0 and c[i] == c[i]]
    if not vis:
        raise ValueError("no usable bands in 400-700 nm to compare the two instruments")
    mean = sum(vis) / len(vis)
    spread = (max(vis) - min(vis)) / mean

    if abs(mean - 1.0) <= 0.03 and spread <= 0.03:
        verdict = ("The two instruments agree to %.1f %% (mean C = %.3f, spectral "
                   "spread %.1f %%). Use the measured E_d directly."
                   % (100 * abs(mean - 1), mean, 100 * spread))
    elif spread > 0.05:
        verdict = ("C varies by %.1f %% across 400-700 nm (mean %.3f). A spectrally "
                   "STRUCTURED disagreement is not a simple gain offset: suspect a "
                   "stale wavelength or radiometric calibration on one instrument, a "
                   "tilted collector, or a partly shaded panel. Do not just multiply "
                   "it out." % (100 * spread, mean))
    else:
        verdict = ("The two instruments differ by %.1f %% (mean C = %.3f) but the "
                   "disagreement is spectrally flat, so it behaves like a gain offset. "
                   "Decide which absolute scale you trust before correcting; applying C "
                   "re-introduces the panel reflectance you were trying to avoid."
                   % (100 * (mean - 1), mean))

    return {"factor": c, "mean": mean, "spread": spread, "n_bands": len(vis),
            "agree": abs(mean - 1.0) <= 0.03 and spread <= 0.03, "verdict": verdict}


def rrs_from_separate_ed(wavelength, l_target, l_sky, ed_measured, ed_wavelength=None,
                         rho=RHO_MOBLEY1999, residual="none", calibration_factor=None):
    """R_rs using E_d from a SEPARATE irradiance instrument.

        R_rs = (L_t - rho * L_sky) / (C * E_d_measured)

    ``ed_wavelength`` lets the irradiance sensor be on its own grid; it is interpolated
    onto the radiance grid. ``calibration_factor`` is optional and is what
    :func:`cross_calibration_factor` returns as ``factor``; leave it None to trust the
    irradiance sensor's own absolute scale, which is the point of having it.

    No panel is involved anywhere in this path. The panel reflectance, the panel
    levelness and the panel-to-target time lag all drop out, and because the two
    instruments log at the same moment there is no changing-light assumption left to
    violate.
    """
    from .resample import interp_linear

    n = len(wavelength)
    for name, arr in (("l_target", l_target), ("l_sky", l_sky)):
        if len(arr) != n:
            raise ValueError("%s does not match the radiance wavelength grid" % name)

    if ed_wavelength is not None:
        ed = interp_linear(list(ed_wavelength), list(ed_measured), list(wavelength))
    else:
        if len(ed_measured) != n:
            raise ValueError(
                "ed_measured has %d points against %d radiance bands; pass "
                "ed_wavelength." % (len(ed_measured), n))
        ed = list(ed_measured)

    if calibration_factor is not None:
        if len(calibration_factor) != n:
            raise ValueError("calibration_factor does not match the radiance grid")
        ed = [ed[i] * calibration_factor[i] for i in range(n)]

    notes = ["E_d from a separate irradiance instrument. No panel enters this "
             "calculation: panel reflectance, panel levelness and the panel-to-target "
             "time lag are all absent."]
    if calibration_factor is None:
        notes.append("No inter-instrument calibration factor applied. R_rs therefore "
                     "carries the COMBINED absolute calibration uncertainty of the two "
                     "instruments as a multiplicative bias. Run "
                     "cross_calibration_factor() against a panel scan to measure it.")

    rrs, n_bad, n_out = [], 0, 0
    for i in range(n):
        w = l_target[i] - rho * l_sky[i]
        if ed[i] and ed[i] > 0 and ed[i] == ed[i]:
            rrs.append(w / ed[i])
        else:
            rrs.append(float("nan"))
            n_bad += 1
        if ed_wavelength is not None and not (
                min(ed_wavelength) <= wavelength[i] <= max(ed_wavelength)):
            n_out += 1
    if n_bad:
        notes.append("%d bands have zero or invalid E_d." % n_bad)
    if n_out:
        notes.append("%d radiance bands fall OUTSIDE the irradiance sensor's "
                     "wavelength range and are extrapolated to NaN. The two "
                     "instruments do not cover the same spectrum." % n_out)

    off, method, more = residual_correction(wavelength, rrs, residual)
    notes.extend(more)
    rrs = [v - off for v in rrs]
    res = RrsResult(list(wavelength), rrs, list(ed), None, float(rho), float(off),
                    method, notes)
    res.meta = {"ed_source": "separate irradiance instrument",
                "cross_calibrated": calibration_factor is not None}
    return res
