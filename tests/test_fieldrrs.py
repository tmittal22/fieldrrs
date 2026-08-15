"""fieldrrs tests. Standard library only, so they run wherever the package runs."""

import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fieldrrs import (  # noqa: E402
    RHO_MOBLEY1999, bin_spectrum, gaussian_resample, read_sed, residual_correction,
    rho_advice, rrs_from_sed, rrs_three_scan, write_batch_csv, write_rrs_csv,
)
from fieldrrs.sed import guess_role  # noqa: E402
from fieldrrs.solar import (  # noqa: E402
    compass_point, local_to_utc_hours, pointing, solar_position,
)

WL = [350.0 + i for i in range(651)]          # 350-1000 nm at 1 nm


def darwin_float(v, digits=12):
    s = "%.*e" % (digits, v)
    mant, _, exp = s.partition("e")
    return "%se%s%03d" % (mant, exp[0], int(exp[1:]))


def write_sed(path, wl, rad_ref, rad_target, comment="test", reflect_pct=True):
    lines = ["Version: 2.1", "Instrument: NaturaSpecPlus_SN0000",
             "Measurement: REFLECTANCE", "Date: 08/15/2026,08/15/2026",
             "Time: 12:00:00,12:00:30", "Latitude: 40.79", "Longitude: -77.86",
             "Comment: %s" % comment, "Channels: %d" % len(wl)]
    cols = ["Wvl", "Rad. (Ref.)", "Rad. (Target)"]
    data = [wl, rad_ref, rad_target]
    if reflect_pct:
        cols.append("Reflect. %")
        data.append([100.0 * t / r for t, r in zip(rad_target, rad_ref)])
    lines.append("Data:")
    lines.append("\t".join(cols))
    for row in zip(*data):
        lines.append("\t".join(darwin_float(v) for v in row))
    with open(path, "w", encoding="latin-1") as fh:
        fh.write("\n".join(lines))
    return path


def synthetic_station(tmp, rrs_true_fn, rho=RHO_MOBLEY1999, panel_r=0.99):
    """Build water/sky/panel .sed files from a KNOWN R_rs, so closure is checkable."""
    ed = [1.0 * math.exp(-((w - 550.0) ** 2) / (2 * 300.0 ** 2)) + 0.2 for w in WL]
    l_panel = [e * panel_r / math.pi for e in ed]
    l_sky = [0.02 * e / math.pi for e in ed]
    rrs_true = [rrs_true_fn(w) for w in WL]
    l_water = [rrs_true[i] * ed[i] + rho * l_sky[i] for i in range(len(WL))]
    w = write_sed(os.path.join(tmp, "station1_water.sed"), WL, l_panel, l_water, "water")
    s = write_sed(os.path.join(tmp, "station1_sky.sed"), WL, l_panel, l_sky, "sky")
    p = write_sed(os.path.join(tmp, "station1_panel.sed"), WL, l_panel, l_panel, "panel")
    return w, s, p, rrs_true


def clear_water(w):
    return 4e-3 * math.exp(-((w - 490.0) ** 2) / (2 * 60.0 ** 2))


class TestSedReader(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_reads_header_and_columns(self):
        w, _, _, _ = synthetic_station(self.tmp, clear_water)
        s = read_sed(w)
        self.assertEqual(s.header["Instrument"], "NaturaSpecPlus_SN0000")
        self.assertEqual(s.comment, "water")
        self.assertAlmostEqual(s.latitude, 40.79, places=4)
        self.assertEqual(len(s.wavelength), len(WL))
        self.assertIn("rad_target", s.columns)
        self.assertIn("rad_ref", s.columns)

    def test_percent_reflectance_becomes_fraction(self):
        p = os.path.join(self.tmp, "u.sed")
        with open(p, "w") as fh:
            fh.write("Data:\nWvl\tReflect. %\n400.0\t50.0\n401.0\t50.0\n")
        self.assertAlmostEqual(read_sed(p).reflectance[0], 0.5)

    def test_unit_scale_variant(self):
        p = os.path.join(self.tmp, "v.sed")
        with open(p, "w") as fh:
            fh.write("Data:\nWvl\tReflect. [1.0]\n400.0\t0.5\n401.0\t0.5\n")
        self.assertAlmostEqual(read_sed(p).reflectance[0], 0.5)

    def test_irr_alias(self):
        """DARWin writes 'Irr. (Ref.)'; parsers expecting 'Irrad.' drop the column."""
        p = os.path.join(self.tmp, "i.sed")
        with open(p, "w") as fh:
            fh.write("Data:\nWvl\tIrr. (Ref.)\n400.0\t1.5\n401.0\t1.6\n")
        self.assertIn("irr_ref", read_sed(p).columns)

    def test_missing_data_marker_is_explained(self):
        p = os.path.join(self.tmp, "bad.sed")
        with open(p, "w") as fh:
            fh.write("Version: 2.1\nWvl\tReflect. %\n400\t50\n")
        with self.assertRaises(ValueError) as cm:
            read_sed(p)
        self.assertIn("Data:", str(cm.exception))

    def test_role_guessing(self):
        w, s, p, _ = synthetic_station(self.tmp, clear_water)
        self.assertEqual(guess_role(read_sed(w)), "water")
        self.assertEqual(guess_role(read_sed(s)), "sky")
        self.assertEqual(guess_role(read_sed(p)), "panel")


class TestClosure(unittest.TestCase):
    """Build files from a known R_rs, read them back, check the physics inverts."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_three_file_closure(self):
        w, s, p, truth = synthetic_station(self.tmp, clear_water)
        res = rrs_from_sed(read_sed(w), read_sed(s), read_sed(p),
                           panel_reflectance=0.99, rho=RHO_MOBLEY1999)
        for got, want in zip(res.rrs, truth):
            self.assertAlmostEqual(got, want, places=9)

    def test_two_file_closure_uses_the_reference_column(self):
        """The DARWin-native workflow: panel is the REFERENCE, no third file."""
        w, s, _, truth = synthetic_station(self.tmp, clear_water)
        res = rrs_from_sed(read_sed(w), read_sed(s), None,
                           panel_reflectance=0.99, rho=RHO_MOBLEY1999)
        for got, want in zip(res.rrs, truth):
            self.assertAlmostEqual(got, want, places=9)

    def test_reflectance_path_matches_radiance_path(self):
        w, s, _, _ = synthetic_station(self.tmp, clear_water)
        a = rrs_from_sed(read_sed(w), read_sed(s), None, 0.99, source="radiance")
        b = rrs_from_sed(read_sed(w), read_sed(s), None, 0.99, source="reflectance")
        for x, y in zip(a.rrs, b.rrs):
            self.assertAlmostEqual(x, y, places=10)

    def test_control_wrong_rho_biases_low(self):
        """If rho did not matter, the closure test would not be testing the glint step."""
        w, s, _, truth = synthetic_station(self.tmp, clear_water)
        res = rrs_from_sed(read_sed(w), read_sed(s), None, 0.99, rho=0.05)
        self.assertFalse(all(abs(g - t) < 1e-9 for g, t in zip(res.rrs, truth)))
        self.assertTrue(all(g <= t + 1e-12 for g, t in zip(res.rrs, truth)))

    def test_panel_reflectance_scales_linearly(self):
        w, s, _, _ = synthetic_station(self.tmp, clear_water)
        a = rrs_from_sed(read_sed(w), read_sed(s), None, 0.99)
        b = rrs_from_sed(read_sed(w), read_sed(s), None, 0.50)
        for x, y in zip(a.rrs, b.rrs):
            self.assertAlmostEqual(y, x * 0.50 / 0.99, places=12)

    def test_material_over_subtraction_is_reported(self):
        """An absurd sky scan must be flagged."""
        res = rrs_three_scan(WL, [1.0] * len(WL), [100.0] * len(WL), [10.0] * len(WL))
        self.assertTrue(any("materially more than" in n for n in res.notes))

    def test_a_clean_measurement_produces_no_over_subtraction_warning(self):
        """Control, and the reason the threshold exists. In the far NIR the true signal
        is ~0, so L_w sits at +/- rounding; a bare `L_w < 0` test fires on a perfectly
        good spectrum and trains the operator to ignore warnings."""
        tmp = tempfile.mkdtemp()
        w, s, _, _ = synthetic_station(tmp, clear_water)
        res = rrs_from_sed(read_sed(w), read_sed(s), None, 0.99, rho=RHO_MOBLEY1999)
        self.assertFalse(any("materially more than" in n for n in res.notes),
                         "clean synthetic station should not warn: %s" % res.notes)

    def test_mismatched_grids_refused(self):
        with self.assertRaises(ValueError):
            rrs_three_scan([400.0, 401.0], [1.0, 1.0], [1.0], [1.0, 1.0])


class TestGlint(unittest.TestCase):
    def test_nir_zero_recovers_a_pure_offset(self):
        wl = [float(w) for w in range(400, 900)]
        true = [3e-3 * math.exp(-((w - 500) ** 2) / (2 * 50.0 ** 2)) if w <= 700 else 0.0
                for w in wl]
        off, method, notes = residual_correction(wl, [v + 5e-4 for v in true], "nir_zero")
        self.assertAlmostEqual(off, 5e-4, places=12)
        self.assertIn("INVALID in turbid", notes[0])

    def test_nir_zero_deletes_real_signal_in_turbid_water(self):
        """The documented failure mode, made explicit rather than left as a warning."""
        wl = [float(w) for w in range(400, 900)]
        off, _, _ = residual_correction(wl, [2e-3] * len(wl), "nir_zero")
        self.assertAlmostEqual(off, 2e-3, places=12)

    def test_nir_similarity_needs_the_bands(self):
        wl = [float(w) for w in range(400, 700)]
        with self.assertRaises(ValueError):
            residual_correction(wl, [1e-3] * len(wl), "nir_similarity")

    def test_none_is_a_no_op(self):
        self.assertEqual(residual_correction([400.0], [1.0], "none")[0], 0.0)

    def test_rho_advice_refuses_to_invent_a_high_wind_value(self):
        val, msg = rho_advice(12.0)
        self.assertIsNone(val)
        self.assertIn("NOT bundled", msg)
        val, msg = rho_advice(3.0)
        self.assertEqual(val, RHO_MOBLEY1999)


class TestResample(unittest.TestCase):
    def test_bin_reports_empty_bins_as_nan(self):
        wl = [float(w) for w in range(400, 500)]
        _, out, n = bin_spectrum(wl, [1.0] * len(wl), [450.0, 900.0], width=10.0)
        self.assertAlmostEqual(out[0], 1.0)
        self.assertNotEqual(out[1], out[1])          # NaN
        self.assertEqual(n[1], 0)

    def test_gaussian_preserves_a_constant(self):
        wl = [float(w) for w in range(400, 700)]
        _, out = gaussian_resample(wl, [2.5] * len(wl), [450.0, 550.0], fwhm=10.0)
        for v in out:
            self.assertAlmostEqual(v, 2.5, places=9)

    def test_gaussian_refuses_edge_centres(self):
        wl = [float(w) for w in range(400, 700)]
        _, out = gaussian_resample(wl, [1.0] * len(wl), [395.0], fwhm=10.0)
        self.assertNotEqual(out[0], out[0])


class TestExport(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_csv_carries_the_warnings(self):
        w, s, _, _ = synthetic_station(self.tmp, clear_water)
        res = rrs_from_sed(read_sed(w), read_sed(s), None, 0.99, rho=0.5)
        p = write_rrs_csv(os.path.join(self.tmp, "o.csv"), res)
        text = open(p).read()
        self.assertIn("# WARNING", text)
        self.assertIn("rho", text)

    def test_batch_csv_is_one_row_per_station(self):
        w, s, _, _ = synthetic_station(self.tmp, clear_water)
        res = rrs_from_sed(read_sed(w), read_sed(s), None, 0.99)
        p = write_batch_csv(os.path.join(self.tmp, "b.csv"),
                            [("st1", res), ("st2", res)],
                            centres=[443.0, 490.0, 555.0])
        rows = [ln for ln in open(p).read().splitlines() if ln.strip()]
        self.assertEqual(len(rows), 3)
        self.assertTrue(rows[0].startswith("station"))


class TestSolarPosition(unittest.TestCase):
    """Graded against the ANALYTIC extremes, not against another implementation.

    At local solar noon the solar zenith angle equals |latitude - declination|, and at
    the solstices the declination is +/-23.44 deg exactly. Those are geometry, so they
    are a real oracle rather than a second opinion.
    """

    LAT, LON = 40.7934, -77.86          # Penn State

    @staticmethod
    def noon_utc(lat, lon, y, m, d):
        h = 12.0 - lon / 15.0
        for _ in range(5):
            sp = solar_position(lat, lon, y, m, d, h)
            h = (720.0 - sp.eq_time_min - 4.0 * lon) / 60.0
        return h

    def test_june_solstice_noon_zenith(self):
        h = self.noon_utc(self.LAT, self.LON, 2026, 6, 21)
        sp = solar_position(self.LAT, self.LON, 2026, 6, 21, h)
        self.assertAlmostEqual(sp.zenith, self.LAT - 23.44, delta=0.05)
        self.assertAlmostEqual(sp.declination, 23.44, delta=0.02)

    def test_december_solstice_noon_zenith(self):
        h = self.noon_utc(self.LAT, self.LON, 2026, 12, 21)
        sp = solar_position(self.LAT, self.LON, 2026, 12, 21, h)
        self.assertAlmostEqual(sp.zenith, self.LAT + 23.44, delta=0.05)

    def test_overhead_at_tropic_of_cancer_on_june_solstice(self):
        h = self.noon_utc(23.44, 0.0, 2026, 6, 21)
        sp = solar_position(23.44, 0.0, 2026, 6, 21, h)
        self.assertLess(sp.zenith, 0.1)

    def test_equinox_declination_is_near_zero(self):
        h = self.noon_utc(0.0, 0.0, 2026, 3, 20)
        sp = solar_position(0.0, 0.0, 2026, 3, 20, h)
        self.assertLess(abs(sp.declination), 0.5)
        self.assertLess(sp.zenith, 0.5)

    def test_azimuth_is_south_at_northern_solar_noon(self):
        h = self.noon_utc(self.LAT, self.LON, 2026, 6, 21)
        sp = solar_position(self.LAT, self.LON, 2026, 6, 21, h)
        self.assertAlmostEqual(sp.azimuth, 180.0, delta=1.0)
        self.assertEqual(sp.compass, "S")

    def test_sun_moves_east_to_west_through_the_day(self):
        az = [solar_position(self.LAT, self.LON, 2026, 8, 15, float(h)).azimuth
              for h in (13, 15, 17, 19, 21)]
        self.assertEqual(az, sorted(az), "azimuth must increase E -> S -> W")

    def test_sun_below_horizon_is_reported(self):
        sp = solar_position(self.LAT, self.LON, 2026, 8, 15, 5.0)   # ~01:00 local
        self.assertLess(sp.elevation, 0)
        self.assertIn("below the horizon", sp.advice())

    def test_pointing_gives_two_bearings_135_either_side(self):
        p = pointing(self.LAT, self.LON, 2026, 8, 15, 16.0)
        for b in (p.bearing_ccw, p.bearing_cw):
            sep = abs((b - p.sun.azimuth + 180) % 360 - 180)
            self.assertAlmostEqual(sep, 135.0, delta=1e-6)

    def test_tilt_from_horizontal_is_the_complement(self):
        """40 deg from nadir is 50 deg below horizontal. This is the number you aim."""
        p = pointing(self.LAT, self.LON, 2026, 8, 15, 16.0, view_zenith=40.0)
        self.assertAlmostEqual(p.tilt_from_horizontal, 50.0)
        self.assertAlmostEqual(p.sky_elevation, 50.0)

    def test_usable_window_flags_high_and_low_sun(self):
        self.assertFalse(solar_position(self.LAT, self.LON, 2026, 8, 15, 12.0).usable)
        self.assertTrue(solar_position(self.LAT, self.LON, 2026, 8, 15, 16.0).usable)

    def test_local_to_utc_handles_the_midnight_wrap(self):
        h, shift = local_to_utc_hours(22, 30, -4.0)      # 22:30 EDT -> 02:30 UTC next day
        self.assertAlmostEqual(h, 2.5)
        self.assertEqual(shift, 1)
        h, shift = local_to_utc_hours(1, 0, 5.5)         # 01:00 IST -> 19:30 UTC prev day
        self.assertAlmostEqual(h, 19.5)
        self.assertEqual(shift, -1)

    def test_compass_points(self):
        self.assertEqual(compass_point(0), "N")
        self.assertEqual(compass_point(90), "E")
        self.assertEqual(compass_point(180), "S")
        self.assertEqual(compass_point(315), "NW")


class TestNoThirdPartyImports(unittest.TestCase):
    def test_package_is_pure_stdlib(self):
        """The whole point: this must run on a bare python.org install with no pip."""
        here = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "fieldrrs")
        banned = ("numpy", "scipy", "matplotlib", "pandas")
        for fn in os.listdir(here):
            if not fn.endswith(".py"):
                continue
            src = open(os.path.join(here, fn)).read()
            for b in banned:
                self.assertNotIn("import %s" % b, src, "%s imports %s" % (fn, b))
                self.assertNotIn("from %s" % b, src, "%s imports %s" % (fn, b))


if __name__ == "__main__":
    unittest.main(verbosity=2)
