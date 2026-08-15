"""Solar position, and the pointing geometry that follows from it. Pure stdlib.

Implements the NOAA Solar Calculator algorithm (after Meeus, *Astronomical Algorithms*,
2nd ed.), which is accurate to well under a degree over 1900-2100. That is far better
than you can aim a hand-held foreoptic, which is the point: the limiting error in the
field is your arm, not the ephemeris.

Everything here is deterministic. Given the same latitude, longitude and UTC timestamp
it returns the same answer, so it is tested against the analytic solstice and equinox
extremes rather than against another implementation.
"""

from __future__ import annotations

import math

__all__ = [
    "solar_position", "pointing", "SolarPosition", "Pointing",
    "DEFAULT_VIEW_ZENITH", "DEFAULT_RELATIVE_AZIMUTH",
    "SOLAR_ZENITH_MIN", "SOLAR_ZENITH_MAX",
]

#: Mobley (1999) recommended viewing geometry, endorsed by the IOCCG Protocol Series
#: (2019), "Protocols for Satellite Ocean Colour Data Validation: In Situ Optical
#: Radiometry" v3.0, Ch. 5: "a viewing angle theta of 40 deg and a relative azimuth phi
#: of 135 deg are the most appropriate to minimize sun-glint perturbations".
DEFAULT_VIEW_ZENITH = 40.0        # degrees from NADIR (straight down)
DEFAULT_RELATIVE_AZIMUTH = 135.0  # degrees in azimuth away from the sun

#: The same IOCCG chapter immediately qualifies 135 deg: it "may easily become the
#: source of perturbations ... because the radiometer necessarily looks at the sea close
#: to the deployment structure or at its shadow", worsening at large solar zenith, "which
#: would suggest that phi = 90 deg is a better solution". Their Fig. 5.1 shows 40/90 as
#: the geometry "commonly applied". Use 90 from a boat or pier; 135 is better from a
#: small or shore-based setup where nothing is in the way.
ALT_RELATIVE_AZIMUTH = 90.0

#: Above ~20 deg full-angle the sky-glint contribution varies too much across the field
#: of view (IOCCG v3.0 Ch. 5, "Field-of-view").
MAX_RECOMMENDED_FOV_DEG = 20.0

#: Usable solar-zenith window for above-water radiometry. Outside it, either the sun is
#: high enough that glint is hard to avoid, or low enough that the signal is weak and
#: the atmospheric path long.
SOLAR_ZENITH_MIN = 20.0
SOLAR_ZENITH_MAX = 60.0

_COMPASS = ("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
            "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")


class SolarPosition(object):
    def __init__(self, zenith, azimuth, elevation, declination, eq_time_min):
        self.zenith = zenith            # deg from vertical
        self.azimuth = azimuth          # deg clockwise from true North
        self.elevation = elevation      # deg above horizon
        self.declination = declination
        self.eq_time_min = eq_time_min

    @property
    def compass(self):
        return compass_point(self.azimuth)

    @property
    def usable(self):
        return SOLAR_ZENITH_MIN <= self.zenith <= SOLAR_ZENITH_MAX

    def advice(self):
        if self.elevation <= 0:
            return ("The sun is below the horizon (elevation %.1f deg). There is no "
                    "measurement to make." % self.elevation)
        if self.zenith < SOLAR_ZENITH_MIN:
            return ("Solar zenith %.1f deg: the sun is very high. Sun glint is hard to "
                    "avoid at any azimuth and rho is less reliable. Prefer a solar "
                    "zenith of %.0f-%.0f deg." % (self.zenith, SOLAR_ZENITH_MIN,
                                                  SOLAR_ZENITH_MAX))
        if self.zenith > SOLAR_ZENITH_MAX:
            return ("Solar zenith %.1f deg: the sun is low, so the signal is weak, the "
                    "atmospheric path is long, and rho grows. Usable but note it."
                    % self.zenith)
        return ("Solar zenith %.1f deg: inside the %.0f-%.0f deg window. Good."
                % (self.zenith, SOLAR_ZENITH_MIN, SOLAR_ZENITH_MAX))

    def __repr__(self):
        return ("SolarPosition(zenith=%.2f, azimuth=%.2f (%s), elevation=%.2f)"
                % (self.zenith, self.azimuth, self.compass, self.elevation))


class Pointing(object):
    """Where to aim, derived from the sun position and the chosen geometry."""

    def __init__(self, sun, view_zenith, relative_azimuth):
        self.sun = sun
        self.view_zenith = view_zenith
        self.relative_azimuth = relative_azimuth
        self.bearing_ccw = (sun.azimuth - relative_azimuth) % 360.0
        self.bearing_cw = (sun.azimuth + relative_azimuth) % 360.0

    @property
    def tilt_from_horizontal(self):
        """Instrument tilt below horizontal for the water scan.

        The water scan is specified as an angle from NADIR, but you aim relative to the
        horizon, so this is the number you actually set: 40 deg from nadir is 50 deg
        below horizontal.
        """
        return 90.0 - self.view_zenith

    @property
    def sky_elevation(self):
        """Elevation above the horizon for the sky scan. Mirror of the water scan."""
        return 90.0 - self.view_zenith

    def describe(self, declination=None):
        """Pointing instructions. Pass ``declination`` to also give phone-compass
        (magnetic) bearings, which is what the phone will actually read."""
        lines = [
            "SUN is at %.0f deg TRUE (%s), %.0f deg above the horizon."
            % (self.sun.azimuth, self.sun.compass, self.sun.elevation),
            "WATER: aim %.0f deg (%s) or %.0f deg (%s) TRUE, tilted %.0f deg BELOW "
            "horizontal."
            % (self.bearing_ccw, compass_point(self.bearing_ccw),
               self.bearing_cw, compass_point(self.bearing_cw),
               self.tilt_from_horizontal),
            "SKY:   same bearing you chose, tilted %.0f deg ABOVE horizontal."
            % self.sky_elevation,
            "PANEL: level, face up, viewed from straight down.",
        ]
        if declination is not None:
            lines.append(
                "PHONE COMPASS (magnetic, declination %+.1f deg): sun should read "
                "%.0f, water bearings %.0f or %.0f."
                % (declination,
                   magnetic_from_true(self.sun.azimuth, declination),
                   magnetic_from_true(self.bearing_ccw, declination),
                   magnetic_from_true(self.bearing_cw, declination)))
        return lines


def compass_point(azimuth):
    return _COMPASS[int((azimuth % 360.0) / 22.5 + 0.5) % 16]


def declination_from_sun_sighting(true_solar_azimuth, compass_reading_of_sun):
    """Calibrate a phone compass against the sun. Returns magnetic declination, degrees.

    A phone compass reads MAGNETIC bearing; the solar azimuth computed here is TRUE
    bearing. The two differ by the local magnetic declination, which is roughly -10 deg
    in Pennsylvania and varies by tens of degrees elsewhere. Pointing at a bearing that
    is 10 deg wrong is a real error in the 135 deg relative azimuth.

    You do not need a lookup table for this. Point the phone at the sun, read the
    magnetic bearing, and pass it here with the computed true solar azimuth:

        declination = true_solar_azimuth - compass_reading_of_sun

    That single sighting absorbs the local declination AND any constant offset in the
    phone's magnetometer, which a published declination value does not. Redo it if you
    move a long way, or after being near anything ferrous.

    Do not sight the sun through the instrument optics, and do not look at it directly.
    """
    return ((true_solar_azimuth - compass_reading_of_sun + 180.0) % 360.0) - 180.0


def true_from_magnetic(magnetic_bearing, declination):
    """Convert a phone-compass (magnetic) bearing to a true bearing."""
    return (magnetic_bearing + declination) % 360.0


def magnetic_from_true(true_bearing, declination):
    """Convert a true bearing to what the phone compass should read."""
    return (true_bearing - declination) % 360.0


def julian_day(year, month, day, hour_utc=0.0):
    """Julian Day for a Gregorian calendar date and fractional UTC hour."""
    if month <= 2:
        year -= 1
        month += 12
    a = math.floor(year / 100.0)
    b = 2 - a + math.floor(a / 4.0)
    jd = (math.floor(365.25 * (year + 4716))
          + math.floor(30.6001 * (month + 1))
          + day + b - 1524.5)
    return jd + hour_utc / 24.0


def solar_position(latitude, longitude, year, month, day, hour_utc):
    """Solar zenith and azimuth. NOAA Solar Calculator algorithm.

    latitude  : degrees North positive
    longitude : degrees EAST positive (so 77.86 W is -77.86)
    hour_utc  : fractional hours UTC, e.g. 14.5 for 14:30 UTC
    """
    jd = julian_day(year, month, day, hour_utc)
    t = (jd - 2451545.0) / 36525.0

    l0 = (280.46646 + t * (36000.76983 + t * 0.0003032)) % 360.0
    m = 357.52911 + t * (35999.05029 - 0.0001537 * t)
    e = 0.016708634 - t * (0.000042037 + 0.0000001267 * t)

    mrad = math.radians(m)
    c = (math.sin(mrad) * (1.914602 - t * (0.004817 + 0.000014 * t))
         + math.sin(2 * mrad) * (0.019993 - 0.000101 * t)
         + math.sin(3 * mrad) * 0.000289)

    true_long = l0 + c
    omega = 125.04 - 1934.136 * t
    app_long = true_long - 0.00569 - 0.00478 * math.sin(math.radians(omega))

    mean_obliq = 23.0 + (26.0 + ((21.448 - t * (46.815 + t * (0.00059 - t * 0.001813))))
                         / 60.0) / 60.0
    obliq_corr = mean_obliq + 0.00256 * math.cos(math.radians(omega))

    decl = math.degrees(math.asin(math.sin(math.radians(obliq_corr))
                                  * math.sin(math.radians(app_long))))

    y = math.tan(math.radians(obliq_corr / 2.0)) ** 2
    l0r = math.radians(l0)
    eq_time = 4.0 * math.degrees(
        y * math.sin(2 * l0r)
        - 2.0 * e * math.sin(mrad)
        + 4.0 * e * y * math.sin(mrad) * math.cos(2 * l0r)
        - 0.5 * y * y * math.sin(4 * l0r)
        - 1.25 * e * e * math.sin(2 * mrad))

    # True solar time in minutes; input is UTC so the timezone offset is zero.
    true_solar_time = (hour_utc * 60.0 + eq_time + 4.0 * longitude) % 1440.0
    hour_angle = (true_solar_time / 4.0) - 180.0
    if hour_angle < -180.0:
        hour_angle += 360.0

    latr = math.radians(latitude)
    declr = math.radians(decl)
    har = math.radians(hour_angle)

    cos_zen = (math.sin(latr) * math.sin(declr)
               + math.cos(latr) * math.cos(declr) * math.cos(har))
    cos_zen = max(-1.0, min(1.0, cos_zen))
    zenith = math.degrees(math.acos(cos_zen))
    elevation = 90.0 - zenith

    sin_zen = math.sin(math.radians(zenith))
    if abs(sin_zen) < 1e-9 or abs(math.cos(latr)) < 1e-9:
        azimuth = 180.0
    else:
        num = (math.sin(latr) * cos_zen) - math.sin(declr)
        den = math.cos(latr) * sin_zen
        arg = max(-1.0, min(1.0, num / den))
        az = math.degrees(math.acos(arg))
        azimuth = (az + 180.0) % 360.0 if hour_angle > 0 else (540.0 - az) % 360.0

    return SolarPosition(zenith, azimuth, elevation, decl, eq_time)


def pointing(latitude, longitude, year, month, day, hour_utc,
             view_zenith=DEFAULT_VIEW_ZENITH,
             relative_azimuth=DEFAULT_RELATIVE_AZIMUTH):
    """Sun position plus the two compass bearings that satisfy the relative azimuth."""
    sun = solar_position(latitude, longitude, year, month, day, hour_utc)
    return Pointing(sun, view_zenith, relative_azimuth)


def local_to_utc_hours(hour, minute, utc_offset_hours):
    """Local clock time to fractional UTC hours. Returns (hours_utc, day_shift).

    ``day_shift`` is -1, 0 or +1 and must be applied to the calendar date, which is the
    part people get wrong when the offset pushes across midnight.
    """
    h = hour + minute / 60.0 - utc_offset_hours
    shift = 0
    while h < 0:
        h += 24.0
        shift -= 1
    while h >= 24.0:
        h -= 24.0
        shift += 1
    return h, shift
