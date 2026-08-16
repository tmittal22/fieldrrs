"""fieldrrs tests. Standard library only, so they run wherever the package runs."""

import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fieldrrs import (  # noqa: E402
    RHO_MOBLEY1999,
    bin_spectrum,
    cross_calibration_factor,
    ed_stability,
    gaussian_resample,
    integrated_irradiance,
    overcast_notes,
    par_from_ed,
    read_sed,
    residual_correction,
    rho_advice,
    rrs_from_sed,
    rrs_from_separate_ed,
    rrs_three_scan,
    write_batch_csv,
    write_rrs_csv,
)
from fieldrrs.sed import guess_role  # noqa: E402
from fieldrrs.solar import (  # noqa: E402
    ALT_RELATIVE_AZIMUTH, compass_point, declination_from_sun_sighting,
    local_to_utc_hours, magnetic_from_true, pointing, solar_position,
    true_from_magnetic,
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

    def test_utf8_bom_does_not_corrupt_the_first_header_key(self):
        """A BOM decoded as latin-1 glues three characters onto the first key, so
        'Version' becomes unfindable. Caught by running check_sed.py on a BOM file."""
        p = os.path.join(self.tmp, "bom.sed")
        body = "Version: 2.1\nComment: x\nData:\nWvl\tReflect. %\n400.0\t50.0\n401.0\t50.0\n"
        with open(p, "w", encoding="utf-8-sig") as fh:
            fh.write(body)
        s = read_sed(p)
        self.assertIn("Version", s.header)
        self.assertEqual(s.header["Version"], "2.1")

    def test_crlf_and_trailing_space_on_the_data_marker(self):
        """Real Windows exports have CRLF endings, and 'Data:' may carry trailing
        whitespace. Neither may break the parse."""
        p = os.path.join(self.tmp, "crlf.sed")
        with open(p, "w", newline="") as fh:
            fh.write("Comment: x\r\nData:  \r\nWvl\tReflect. %\r\n400.0\t50.0\r\n"
                     "401.0\t50.0\r\n")
        s = read_sed(p)
        self.assertEqual(len(s.wavelength), 2)
        self.assertAlmostEqual(s.reflectance[0], 0.5)

    def test_unknown_extra_columns_are_ignored_not_fatal(self):
        p = os.path.join(self.tmp, "extra.sed")
        with open(p, "w") as fh:
            fh.write("Data:\nWvl\tTgt./Ref. %\tReflect. %\tSomethingNew\n"
                     "400.0\t50.0\t50.0\t1.0\n401.0\t50.0\t50.0\t1.0\n")
        s = read_sed(p)
        self.assertAlmostEqual(s.reflectance[0], 0.5)
        self.assertIn("SomethingNew", s.raw_columns)

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


REAL_HEADER = (
    "Comment:\r\n"
    "Version: 2.4 [2.1.9284]\r\n"
    "File Name: C:\\Users\\W0168\\Documents\\SpectralEvolution\\x.sed\r\n"
    "<Metadata>\r\n"
    "USER_FIELD1: Range: 3.837m\r\n"
    "USER_FIELD2: Tilt (X): + 4.8\u00b0\r\n"
    "USER_FIELD3: Tilt (Y): +28.9\u00b0\r\n"
    "USER_FIELD4: Solar Angle: 32.27\u00b0\r\n"
    "</Metadata>\r\n"
    "Instrument: NaturaSpecPlus_SN25494G1 [3]\r\n"
    "Measurement: REFLECTANCE\r\n"
    "Date: 08/07/2025,08/07/2025\r\n"
    "Time: 05:56:02.07,06:07:42.07\r\n"
    "Foreoptic: FLENS8 {RADIANCE},FLENS8 {RADIANCE}\r\n"
    "Units: W/m^2/sr/nm\r\n"
    "Latitude: 40.89674\r\n"
    "Longitude: -78.20499\r\n"
    "GPS Time: 13:08:18\r\n"
    "Channels: 3\r\n"
    "Columns [4]:\r\n"
    "Data:\r\n"
    "Wvl\tRad. (Ref.)\tRad. (Target)\tReflect. %\r\n"
    "350.0\t1.119097e-002\t4.509495e-005\t0.3794\r\n"
    "351.0\t1.152954e-002\t4.620102e-005\t0.3751\r\n"
    "352.0\t1.160930e-002\t4.491070e-005\t0.3689\r\n"
)


class TestRealFileStructure(unittest.TestCase):
    """Pins everything a REAL NaturaSpecPlus_SN25494G1 export (2025-08-07) revealed.

    The reader was originally built from the DARWin2.exe string table, which did not
    expose the <Metadata> USER_FIELD block, the 'Columns [N]:' line, or the fact that
    the file is UTF-8 with degree signs. All four column labels were correctly inferred;
    these tests keep the rest from regressing.
    """

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.p = os.path.join(self.tmp, "NaturaSpecPlus_SN25494G1_00050.sed")
        with open(self.p, "w", encoding="utf-8", newline="") as fh:
            fh.write(REAL_HEADER)
        self.s = read_sed(self.p)

    def test_all_four_real_columns_are_recognised(self):
        for c in ("wavelength", "rad_ref", "rad_target", "reflectance"):
            self.assertIn(c, self.s.columns)
        self.assertEqual(self.s.reflectance_scale, 100.0)

    def test_degree_signs_survive_decoding(self):
        """UTF-8 decoded as latin-1 turns every degree sign into 'Â°'."""
        self.assertNotIn("\u00c2", self.s.header["USER_FIELD2"])
        self.assertIn("\u00b0", self.s.header["USER_FIELD2"])

    def test_metadata_block_is_parsed(self):
        uf = self.s.user_fields
        self.assertEqual(set(uf), {"Range", "Tilt (X)", "Tilt (Y)", "Solar Angle"})
        self.assertAlmostEqual(self.s.solar_elevation_deg, 32.27)
        self.assertAlmostEqual(self.s.tilt_x_deg, 4.8)
        self.assertAlmostEqual(self.s.tilt_y_deg, 28.9)
        self.assertAlmostEqual(self.s.range_m, 3.837)

    def test_foreoptic_fov(self):
        """FLENS8 is an 8 deg lens, inside the IOCCG <=20 deg guidance."""
        self.assertAlmostEqual(self.s.fov_deg, 8.0)

    def test_gps_time_parses_to_utc_hours(self):
        self.assertAlmostEqual(self.s.gps_time, 13 + 8 / 60.0 + 18 / 3600.0, places=9)

    def test_solar_angle_is_elevation_not_zenith(self):
        """The instrument's 'Solar Angle' matches computed ELEVATION, not zenith.

        At the logged GPS fix the computed elevation is 31.17 deg against the reported
        32.27 deg; the 1.10 deg gap is the 5.9 min between the fix and the scan. The
        zenith at that moment is 58.8 deg, nowhere near, which settles the ambiguity.
        """
        sp = solar_position(self.s.latitude, self.s.longitude, 2025, 8, 7,
                            self.s.gps_time)
        self.assertAlmostEqual(sp.elevation, self.s.solar_elevation_deg, delta=1.5)
        self.assertGreater(abs(sp.zenith - self.s.solar_elevation_deg), 20.0)

    def test_empty_comment_and_numeric_filename_give_no_role(self):
        """Real files have an empty Comment and a serial-numbered filename, so the
        auto-role guess CANNOT work. The operator must assign, and the GUI asks."""
        self.assertEqual(self.s.comment, "")
        self.assertEqual(guess_role(self.s), "unassigned")


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

    def test_sun_sighting_recovers_declination(self):
        """Point the phone at the sun, read magnetic, get declination. No lookup table."""
        p = pointing(self.LAT, self.LON, 2026, 8, 15, 16.0)
        true_dec = -11.0
        reading = magnetic_from_true(p.sun.azimuth, true_dec)
        self.assertAlmostEqual(
            declination_from_sun_sighting(p.sun.azimuth, reading), true_dec, places=9)

    def test_magnetic_true_round_trip(self):
        for dec in (-11.0, 0.0, +14.5):
            for b in (5.0, 175.0, 359.0):
                self.assertAlmostEqual(
                    true_from_magnetic(magnetic_from_true(b, dec), dec), b, places=9)

    def test_declination_wraps_the_short_way(self):
        """Sun at 5 deg true, compass reads 355: that is +10, not -350."""
        self.assertAlmostEqual(declination_from_sun_sighting(5.0, 355.0), 10.0, places=9)

    def test_describe_includes_magnetic_when_declination_given(self):
        p = pointing(self.LAT, self.LON, 2026, 8, 15, 16.0)
        self.assertFalse(any("PHONE COMPASS" in x for x in p.describe()))
        self.assertTrue(any("PHONE COMPASS" in x for x in p.describe(declination=-11.0)))

    def test_alternate_90_degree_azimuth_is_supported(self):
        """IOCCG v3.0 Fig 5.1 shows 40/90 as the geometry commonly applied from a
        structure, because 135 points back at the hull or its shadow."""
        p = pointing(self.LAT, self.LON, 2026, 8, 15, 16.0,
                     relative_azimuth=ALT_RELATIVE_AZIMUTH)
        for b in (p.bearing_ccw, p.bearing_cw):
            sep = abs((b - p.sun.azimuth + 180) % 360 - 180)
            self.assertAlmostEqual(sep, 90.0, delta=1e-6)

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

class TestOvercastAndMeasuredEd(unittest.TestCase):
    """Fully overcast work, and E_d measured rather than inferred from a panel."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _station_with_irradiance(self, rrs_fn, rho=RHO_MOBLEY1999):
        """Write scans that also carry a cosine-collector E_d column."""
        ed = [1.0 * math.exp(-((w - 550) ** 2) / (2 * 300.0 ** 2)) + 0.2 for w in WL]
        l_panel = [e * 0.99 / math.pi for e in ed]
        l_sky = [0.02 * e / math.pi for e in ed]
        truth = [rrs_fn(w) for w in WL]
        l_water = [truth[i] * ed[i] + rho * l_sky[i] for i in range(len(WL))]

        def write(path, tgt, irr):
            lines = ["Comment: x", "Data:",
                     "\t".join(["Wvl", "Rad. (Ref.)", "Rad. (Target)", "Irr. (Target)"])]
            for i, w in enumerate(WL):
                lines.append("\t".join(darwin_float(v) for v in
                                       (w, l_panel[i], tgt[i], irr[i])))
            open(path, "w").write("\n".join(lines))
            return path

        w = write(os.path.join(self.tmp, "w.sed"), l_water, ed)
        s = write(os.path.join(self.tmp, "s.sed"), l_sky, ed)
        return w, s, truth

    def test_measured_ed_recovers_the_known_rrs(self):
        """Closure through the irradiance path, independent of any panel."""
        w, s, truth = self._station_with_irradiance(clear_water)
        res = rrs_from_sed(read_sed(w), read_sed(s), None, source="irradiance")
        for got, want in zip(res.rrs, truth):
            self.assertAlmostEqual(got, want, places=9)

    def test_measured_ed_does_not_depend_on_panel_reflectance(self):
        """The whole point: the panel does not enter, so its calibration cannot bias
        the answer. The radiance path scales linearly with it; this must not."""
        w, s, _ = self._station_with_irradiance(clear_water)
        a = rrs_from_sed(read_sed(w), read_sed(s), None, panel_reflectance=0.99,
                         source="irradiance")
        b = rrs_from_sed(read_sed(w), read_sed(s), None, panel_reflectance=0.50,
                         source="irradiance")
        for x, y in zip(a.rrs, b.rrs):
            self.assertAlmostEqual(x, y, places=12)

    def test_missing_irradiance_column_is_explained(self):
        w, s, _, _ = synthetic_station(self.tmp, clear_water)   # radiance only
        with self.assertRaises(KeyError) as cm:
            rrs_from_sed(read_sed(w), read_sed(s), None, source="irradiance")
        self.assertIn("cosine", str(cm.exception))

    def test_overcast_relaxes_the_wind_limit(self):
        """Under a uniform sky rho stops being strongly wind-dependent, so a wind that
        would be refused under clear sky is accepted with a caveat instead."""
        clear_val, _ = rho_advice(12.0, sky="clear")
        over_val, over_msg = rho_advice(12.0, sky="overcast")
        self.assertIsNone(clear_val)
        self.assertEqual(over_val, RHO_MOBLEY1999)
        self.assertIn("BROKEN CLOUD", over_msg)

    def test_overcast_notes_carry_the_load_bearing_warnings(self):
        text = " ".join(overcast_notes())
        self.assertIn("NO SUN GLINT", text)
        self.assertIn("DO NOT apply the BRDF", text)
        self.assertIn("BROKEN CLOUD IS THE WORST CASE", text)
        self.assertIn("No satellite match-up", text)

    def test_unknown_source_names_the_valid_set(self):
        w, s, _, _ = synthetic_station(self.tmp, clear_water)
        with self.assertRaises(ValueError) as cm:
            rrs_from_sed(read_sed(w), read_sed(s), None, source="nonsense")
        self.assertIn("irradiance", str(cm.exception))

class TestAbsoluteRadiometry(unittest.TestCase):
    """Products that need a calibrated spectroradiometer, not a reflectance instrument."""

    def test_par_matches_the_closed_form_on_a_rectangular_spectrum(self):
        """For constant E_d over a narrow band the integral is analytic:
        PAR = E_d * dlambda * lambda_mid / 119.6. Graded against that, not against
        another implementation."""
        wl = [549.0, 550.0, 551.0]
        ed = [1.0, 1.0, 1.0]
        par, n = par_from_ed(wl, ed, lo=549.0, hi=551.0)
        self.assertEqual(n, 3)
        self.assertAlmostEqual(par, 1.0 * 2.0 * 550.0 / 119.6, places=6)

    def test_par_is_a_photon_flux_so_red_light_counts_more(self):
        """Equal ENERGY at 650 nm carries more photons than at 450 nm, in the ratio of
        the wavelengths. This is what makes PAR different from broadband irradiance."""
        blue = par_from_ed([449.0, 450.0, 451.0], [1.0] * 3, 449.0, 451.0)[0]
        red = par_from_ed([649.0, 650.0, 651.0], [1.0] * 3, 649.0, 651.0)[0]
        self.assertAlmostEqual(red / blue, 650.0 / 450.0, places=6)

    def test_par_magnitude_is_physically_sensible(self):
        """Full midday sun is about 2000 umol/m2/s. A model spectrum carrying ~420 W/m2
        across 400-700 nm must land near that or the constant is wrong."""
        wl = [400.0 + i for i in range(301)]
        ed = [1.5 * math.exp(-((w - 550) ** 2) / (2 * 250.0 ** 2)) for w in wl]
        par, _ = par_from_ed(wl, ed)
        self.assertTrue(1500 < par < 2500, par)

    def test_broadband_irradiance_integrates_to_the_rectangle(self):
        self.assertAlmostEqual(
            integrated_irradiance([500.0, 600.0], [2.0, 2.0], 500.0, 600.0), 200.0)

    def test_par_refuses_when_the_irradiance_channel_is_absent(self):
        with self.assertRaises(ValueError) as cm:
            par_from_ed([400.0], [1.0])
        self.assertIn("irradiance channel", str(cm.exception))

    def test_par_ignores_bands_outside_the_window(self):
        wl = [300.0, 500.0, 600.0, 900.0]
        ed = [99.0, 2.0, 2.0, 99.0]
        par, n = par_from_ed(wl, ed)
        self.assertEqual(n, 2)
        self.assertAlmostEqual(par, integrated_irradiance(wl, ed) * 550.0 / 119.6,
                               places=3)


class TestEdStability(unittest.TestCase):
    """Bracketing a station with two E_d scans: the QC a panel alone cannot give."""

    WL = [400.0 + i for i in range(301)]

    def test_stable_light_passes(self):
        r = ed_stability([1.0] * 301, [1.005] * 301, self.WL)
        self.assertTrue(r["stable"])
        self.assertAlmostEqual(r["worst_change"], 0.005, places=6)

    def test_a_cloud_crossing_is_caught_and_signed(self):
        r = ed_stability([1.0] * 301, [0.78] * 301, self.WL)
        self.assertFalse(r["stable"])
        self.assertAlmostEqual(r["after_over_before"], 0.78, places=6)
        self.assertIn("every R_rs here is suspect", r["verdict"])

    def test_only_the_requested_window_is_judged(self):
        """Junk outside 400-700 must not decide the verdict."""
        wl = [300.0, 500.0, 600.0, 900.0]
        r = ed_stability([1.0, 1.0, 1.0, 1.0], [99.0, 1.0, 1.0, 99.0], wl)
        self.assertEqual(r["n_bands"], 2)
        self.assertTrue(r["stable"])

    def test_mismatched_grids_refused(self):
        with self.assertRaises(ValueError):
            ed_stability([1.0, 1.0], [1.0])


class TestTwoInstrumentSetup(unittest.TestCase):
    """A SEPARATE hemispherical irradiance sensor beside the narrow-FOV radiometer.

    Better in the way that matters most, E_d logged simultaneously so the unchanged-light
    assumption disappears, at the cost of an error a single instrument does not have:
    R_rs divides a radiance from one instrument by an irradiance from another, so their
    relative calibration is a direct multiplicative bias.
    """

    RWL = [350.0 + i for i in range(651)]          # radiance instrument, 1 nm
    EWL = [400.0 + 3.3 * i for i in range(152)]    # irradiance instrument, its own grid

    @staticmethod
    def _ed(w):
        return 1.1 * math.exp(-((w - 550) ** 2) / (2 * 300.0 ** 2)) + 0.2

    def _scene(self, gain=1.0):
        panel = [self._ed(w) * 0.99 / math.pi for w in self.RWL]
        sky = [0.02 * self._ed(w) / math.pi for w in self.RWL]
        water = [clear_water(w) * self._ed(w) + RHO_MOBLEY1999 * s
                 for w, s in zip(self.RWL, sky)]
        ed_meas = [self._ed(w) * gain for w in self.EWL]
        return panel, sky, water, ed_meas

    def test_recovers_the_known_rrs_across_mismatched_grids(self):
        """The two instruments do not share a wavelength grid, and must not need to."""
        panel, sky, water, ed = self._scene()
        res = rrs_from_separate_ed(self.RWL, water, sky, ed, self.EWL)
        for i, w in enumerate(self.RWL):
            if 420 <= w <= 880:
                self.assertAlmostEqual(res.rrs[i], clear_water(w), places=5)

    def test_no_panel_enters_the_calculation(self):
        """The whole point: panel reflectance cannot bias this path because it is
        never used. Contrast the panel path, which scales linearly with it."""
        panel, sky, water, ed = self._scene()
        a = rrs_from_separate_ed(self.RWL, water, sky, ed, self.EWL)
        self.assertTrue(any("panel" in n.lower() and "absent" in n.lower()
                            for n in a.notes))

    def test_cross_calibration_detects_a_gain_offset(self):
        """A 6 % offset between instruments must show up as C = 1/1.06 and flat."""
        panel, sky, water, ed = self._scene(gain=1.06)
        c = cross_calibration_factor(self.RWL, panel, ed, self.EWL)
        self.assertAlmostEqual(c["mean"], 1.0 / 1.06, places=3)
        self.assertLess(c["spread"], 0.01)
        self.assertFalse(c["agree"])
        self.assertIn("spectrally flat", c["verdict"])

    def test_a_gain_offset_biases_every_rrs_by_the_same_fraction(self):
        """This is why cross-calibration matters: it is not noise, it is a bias."""
        _, sky, water, ed_ok = self._scene(gain=1.0)
        _, _, _, ed_off = self._scene(gain=1.06)
        a = rrs_from_separate_ed(self.RWL, water, sky, ed_ok, self.EWL)
        b = rrs_from_separate_ed(self.RWL, water, sky, ed_off, self.EWL)
        i = min(range(len(self.RWL)), key=lambda j: abs(self.RWL[j] - 443))
        self.assertAlmostEqual(b.rrs[i] / a.rrs[i], 1.0 / 1.06, places=6)

    def test_matched_instruments_are_declared_to_agree(self):
        panel, sky, water, ed = self._scene()
        c = cross_calibration_factor(self.RWL, panel, ed, self.EWL)
        self.assertTrue(c["agree"])
        self.assertIn("Use the measured E_d directly", c["verdict"])

    def test_applying_the_factor_removes_the_bias(self):
        panel, sky, water, ed = self._scene(gain=1.06)
        c = cross_calibration_factor(self.RWL, panel, ed, self.EWL)
        res = rrs_from_separate_ed(self.RWL, water, sky, ed, self.EWL,
                                   calibration_factor=c["factor"])
        i = min(range(len(self.RWL)), key=lambda j: abs(self.RWL[j] - 443))
        self.assertAlmostEqual(res.rrs[i], clear_water(443.0), places=5)

    def test_uncorrected_path_warns_about_the_combined_calibration(self):
        panel, sky, water, ed = self._scene()
        res = rrs_from_separate_ed(self.RWL, water, sky, ed, self.EWL)
        self.assertTrue(any("COMBINED absolute calibration" in n for n in res.notes))

    def test_bands_outside_the_irradiance_range_are_reported(self):
        """The radiance instrument reaches 350 and 1000 nm; the irradiance one does
        not. Those bands must be flagged, not silently extrapolated."""
        panel, sky, water, ed = self._scene()
        res = rrs_from_separate_ed(self.RWL, water, sky, ed, self.EWL)
        self.assertTrue(any("OUTSIDE the irradiance sensor" in n for n in res.notes))

    def test_mismatched_grid_without_ed_wavelength_is_refused(self):
        panel, sky, water, ed = self._scene()
        with self.assertRaises(ValueError) as cm:
            rrs_from_separate_ed(self.RWL, water, sky, ed)
        self.assertIn("ed_wavelength", str(cm.exception))


if __name__ == "__main__":
    unittest.main(verbosity=2)
