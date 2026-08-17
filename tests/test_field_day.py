"""Tests against the REAL 2026-08-16 Kotzebue dataset shipped in Data_NatureSpec/.

Every other test in this suite uses synthesised spectra, which verify arithmetic but
cannot catch anything about how the instrument actually writes files or how real water
actually looks. These are the counterpart: 60 genuine scans, and the properties asserted
are physical ones that a synthetic fixture could not have told us to expect.

They skip cleanly if the data folder is absent, so the suite still runs in a checkout
without it.
"""

import math
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

DATA = os.path.join(ROOT, "Data_NatureSpec", "2026_Aug_16")

if not os.path.isdir(DATA):
    raise unittest.SkipTest("no field data at %s" % DATA)

from fieldrrs.rrs import RHO_MOBLEY1999                      # noqa: E402
from fieldrrs.sed import read_sed                            # noqa: E402
from process_field_day import (                              # noqa: E402
    SKY, VEG, WATER, band, classify, load, process, station_key,
)
from survey_field_data import read_header_log                # noqa: E402


def _scans():
    if not hasattr(_scans, "cache"):
        _scans.cache = load(DATA)
    return _scans.cache


class TestRealFilesParse(unittest.TestCase):
    def test_all_sixty_scans_read(self):
        self.assertEqual(len(_scans()), 60)

    def test_they_are_direct_energy_radiance_files(self):
        s = _scans()[0]["spec"]
        self.assertEqual(s.header["Measurement"], "DIRECT_ENERGY")
        self.assertIn("RADIANCE", s.header["Radiometric Calibration"])

    def test_the_full_swir_grid_is_read(self):
        """2151 channels, 350-2500 nm. The visible-only assumption would silently
        truncate two thirds of the file."""
        s = _scans()[0]["spec"]
        self.assertEqual(len(s.wavelength), 2151)
        self.assertAlmostEqual(s.wavelength[0], 350.0)
        self.assertAlmostEqual(s.wavelength[-1], 2500.0)

    def test_both_radiance_columns_are_present(self):
        """Rad. (Ref.) is what makes panel=None legal for this dataset."""
        for s in _scans():
            self.assertTrue(s["spec"].has("rad_ref"))
            self.assertTrue(s["spec"].has("rad_target"))

    def test_metadata_block_gives_geometry_and_sun(self):
        s = _scans()[0]["spec"]
        self.assertTrue(30.0 <= s.solar_elevation_deg <= 40.0)
        self.assertGreater(s.range_m, 0.0)
        self.assertIsNotNone(s.tilt_y_deg)

    def test_position_is_arctic_alaska(self):
        s = _scans()[0]["spec"]
        self.assertAlmostEqual(s.latitude, 66.897, places=2)
        self.assertAlmostEqual(s.longitude, -162.603, places=2)


class TestHeaderLogIsMalformed(unittest.TestCase):
    """HeaderLog.csv has an UNQUOTED COMMA inside ScanType.

    'TARGET, ALLOW_SWIR2_UPDATE' is one logical field, so rows carry 23 values against
    22 headers and csv.DictReader shifts everything after it. The damage is silent and
    plausible-looking, which is why it gets its own test.
    """

    def test_rows_have_one_more_field_than_the_header(self):
        import csv
        with open(os.path.join(DATA, "HeaderLog.csv"), newline="") as fh:
            rows = list(csv.reader(fh))
        self.assertEqual(len(rows[1]), len(rows[0]) + 1)

    def test_dictreader_reports_a_temperature_as_the_latitude(self):
        """The control. If this ever stops failing, the file format was fixed and the
        right-anchored parser can be simplified."""
        import csv
        with open(os.path.join(DATA, "HeaderLog.csv"), newline="") as fh:
            row = next(iter(csv.DictReader(fh)))
        self.assertLess(float(row["Latitude"]), 0.0)      # -15.67, a detector temp

    def test_our_reader_recovers_the_real_values(self):
        log = read_header_log(DATA)
        first = log["NaturaSpecPlus_SN25494G1_00000"]
        self.assertEqual(first["ScanType"], "TARGET")
        self.assertRegex(first["GPSTime"], r"^\d{2}:\d{2}:\d{2}$")
        self.assertAlmostEqual(float(first["YTilt"]), 44.29, places=2)

    def test_the_log_agrees_with_the_sed_metadata(self):
        """Two independent records of the same tilt must match, which is what proves
        the right-anchored parse landed on the right columns."""
        log = read_header_log(DATA)
        for s in _scans()[:10]:
            self.assertAlmostEqual(float(log[s["name"]]["YTilt"]),
                                   s["spec"].tilt_y_deg, delta=0.15)


class TestRoleClassification(unittest.TestCase):
    def test_counts(self):
        roles = [s["role"] for s in _scans()]
        self.assertEqual(roles.count(SKY), 23)
        self.assertEqual(roles.count(WATER), 33)
        self.assertEqual(roles.count(VEG), 4)

    def test_tilt_cannot_separate_sky_from_water(self):
        """The justification for classifying on shape instead. If tilt DID separate
        them, the shape classifier would be unnecessary complexity."""
        sky = [s["spec"].tilt_y_deg for s in _scans() if s["role"] == SKY]
        wat = [s["spec"].tilt_y_deg for s in _scans() if s["role"] == WATER]
        self.assertLess(min(sky), max(wat))
        self.assertLess(min(wat), max(sky))          # the ranges overlap, both ways

    def test_sky_is_blue_rising_and_water_is_not(self):
        """Rayleigh goes as lambda^-4, so a clear sky must be blue-dominated. Water
        cannot be, because it is lit by that sky through an absorbing medium."""
        for s in _scans():
            if s["role"] == SKY:
                self.assertGreater(s["diag"]["blue_green"], 1.3)
            elif s["role"] == WATER:
                self.assertLess(s["diag"]["blue_green"], 1.0)

    def test_water_absorbs_the_nir_and_vegetation_does_the_opposite(self):
        wat = [s["diag"]["nir_vis"] for s in _scans() if s["role"] == WATER]
        veg = [s["diag"]["nir_vis"] for s in _scans() if s["role"] == VEG]
        self.assertLess(max(wat), min(veg))
        self.assertGreater(min(veg) / max(wat), 2.0, "margin to vegetation too thin")

    def test_the_four_vegetation_scans_are_the_expected_ones(self):
        self.assertEqual(sorted(s["n"] for s in _scans() if s["role"] == VEG),
                         ["00014", "00015", "00043", "00044"])


class TestStationGrouping(unittest.TestCase):
    def test_five_stations(self):
        self.assertEqual(len({s["station"] for s in _scans()}), 5)

    def test_every_station_has_sky_and_water(self):
        st = {}
        for s in _scans():
            st.setdefault(s["station"], []).append(s["role"])
        for k, roles in st.items():
            self.assertIn(SKY, roles, "station %s has no sky scan" % k)
            self.assertIn(WATER, roles, "station %s has no water scan" % k)

    def test_the_reference_is_constant_within_a_station(self):
        """The grouping IS the reference, so this must hold by construction; it is here
        to catch a rounding change in station_key silently merging two stations."""
        by = {}
        for s in _scans():
            by.setdefault(s["station"], []).append(s["spec"])
        for k, specs in by.items():
            vals = [band(sp, "rad_ref", 450, 650) for sp in specs]
            self.assertLess(max(vals) - min(vals), 1e-6)

    def test_panel_radiance_is_consistent_with_full_sun(self):
        """L_panel = E_d * R / pi. At 0.23-0.31 W m^-2 sr^-1 nm^-1 that implies E_d of
        order 0.7-1.0, which is full daylight. A shaded or mis-scanned panel would be
        an order of magnitude lower and would silently inflate every R_rs."""
        for s in _scans():
            lp = band(s["spec"], "rad_ref", 450, 650)
            ed = lp * math.pi / 0.99
            self.assertTrue(0.5 < ed < 1.5, "implied E_d %.2f from %s" % (ed, s["n"]))


class TestRetrievedRrs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = process(_scans(), RHO_MOBLEY1999, 0.99, "none")

    def test_all_five_stations_produce_rrs(self):
        self.assertEqual(len(self.res), 5)
        self.assertEqual(sum(len(st["rrs"]) for st in self.res), 33)

    def test_rrs_is_positive_across_the_visible(self):
        """Over-subtraction of skylight is the classic above-water failure and shows up
        as negative R_rs in the blue."""
        for st in self.res:
            for r in st["rrs"]:
                vis = [v for w, v in zip(r["res"].wavelength, r["res"].rrs)
                       if 400 <= w <= 700]
                self.assertGreater(min(vis), -1e-4, "negative R_rs in scan %s" % r["n"])

    def test_magnitude_is_sane_for_turbid_coastal_water(self):
        for st in self.res:
            for r in st["rrs"]:
                v = [x for w, x in zip(r["res"].wavelength, r["res"].rrs)
                     if 400 <= w <= 750]
                self.assertTrue(1e-3 < max(v) < 5e-2,
                                "peak R_rs %.4f in scan %s" % (max(v), r["n"]))

    def test_the_spectra_peak_in_the_green_not_the_blue(self):
        """Sediment-laden water. A blue peak here would mean clear oceanic water and
        would contradict everything else about this site."""
        for st in self.res:
            for r in st["rrs"]:
                wl, v = r["res"].wavelength, r["res"].rrs
                i = max((k for k in range(len(wl)) if 400 <= wl[k] <= 750),
                        key=lambda k: v[k])
                self.assertGreater(wl[i], 520.0, "scan %s peaks at %.0f nm"
                                   % (r["n"], wl[i]))

    def test_the_700nm_turbidity_peak_is_present(self):
        """a_w has a local minimum near 700 nm, so strongly backscattering water shows a
        distinct secondary peak there. Its presence is the signature of high SPM."""
        for st in self.res:
            mean = [sum(r["res"].rrs[i] for r in st["rrs"]) / len(st["rrs"])
                    for i in range(len(st["wl"]))]
            wl = st["wl"]
            near700 = max(v for w, v in zip(wl, mean) if 690 <= w <= 715)
            trough = min(v for w, v in zip(wl, mean) if 640 <= w <= 670)
            self.assertGreater(near700, trough,
                               "no 700 nm peak at station %.4f" % st["ref"])

    def test_rrs_collapses_beyond_750nm(self):
        """Liquid water absorption. If the NIR did not collapse, the sky subtraction or
        the panel would be wrong."""
        for st in self.res:
            mean = [sum(r["res"].rrs[i] for r in st["rrs"]) / len(st["rrs"])
                    for i in range(len(st["wl"]))]
            wl = st["wl"]
            green = max(v for w, v in zip(wl, mean) if 550 <= w <= 600)
            nir = sum(v for w, v in zip(wl, mean) if 890 <= w <= 910) / \
                sum(1 for w in wl if 890 <= w <= 910)
            self.assertLess(nir, 0.35 * green)

    def test_replicates_within_a_station_agree(self):
        """Random surface variability only. Large scatter would mean the water changed
        or a scan caught glint."""
        for st in self.res:
            if len(st["rrs"]) < 3:
                continue
            peaks = []
            for r in st["rrs"]:
                v = [x for w, x in zip(r["res"].wavelength, r["res"].rrs)
                     if 550 <= w <= 600]
                peaks.append(max(v))
            spread = (max(peaks) - min(peaks)) / (sum(peaks) / len(peaks))
            self.assertLess(spread, 0.6, "station %.4f replicate spread %.0f%%"
                            % (st["ref"], 100 * spread))


if __name__ == "__main__":
    unittest.main(verbosity=2)
