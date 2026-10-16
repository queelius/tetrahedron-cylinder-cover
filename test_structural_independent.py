"""Failure controls and exact identities for the independent structural checker."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

import structural_independent as check


class EnclosureTests(unittest.TestCase):
    def test_square_root_upper_enclosures(self):
        for value in [F(1), F(4), F(5, 4), F(2), F(1000001, 1000000)]:
            upper = check.sqrt_upper(value)
            self.assertGreaterEqual(upper * upper, value)
            self.assertLess(upper * upper - value, F(1, 10**11))

    def test_invalid_square_root_domain_rejected(self):
        for value in [F(0), F(-1), F(5)]:
            with self.assertRaises(ValueError):
                check.sqrt_upper(value)

    def test_zero_tilt_integral_and_impossible_target(self):
        sectors = check.make_sectors(8, (F(0),))
        bound = check.interval_bound(sectors, F(0), F(0), 4)
        self.assertEqual(bound, F(1, 2))
        check.require_bound(bound, F(1, 2))
        with self.assertRaises(check.VerificationError):
            check.require_bound(bound, F(501, 1000))

    def test_inverted_or_out_of_domain_interval_rejected(self):
        sectors = check.make_sectors(8, (F(0),))
        for lower, upper in [(F(1), F(0)), (F(-1), F(0)), (F(0), F(2))]:
            with self.assertRaises(ValueError):
                check.interval_bound(sectors, lower, upper, 4)

    def test_bad_geometric_family_rejected(self):
        with self.assertRaises(check.VerificationError):
            check.make_sectors(8, (F(100),))
        with self.assertRaises(ValueError):
            check.make_sectors(3, (F(0),))

    def test_test_point_and_other_sector_exclusion(self):
        # Finite identity controls, not a proof over all q,tau.
        positive = check.make_sectors(8, (F(-1, 200), F(27, 1250)))
        full = [check.Sector(-z.r, -z.l, z.a, -z.b)
                for z in reversed(positive)] + list(positive)
        for sector in positive:
            a, b = sector.a, sector.b
            for q in [sector.l + (sector.r-sector.l)/3,
                      sector.l + 2*(sector.r-sector.l)/3]:
                for tau in [F(0), F(1, 2), F(1)]:
                    d = 1 + tau*tau*a*(a-b*q)
                    radial = (1+tau*tau*a*a)/(2*d)
                    x = -q*(1+tau*a)/(2*d)
                    y = q*(1-tau*a)/(2*d)
                    t = radial + tau*b*x
                    self.assertLessEqual(abs(x), 1-t)
                    self.assertLessEqual(abs(y), t)
                    self.assertEqual(y-tau*a*x, q*radial)
                    self.assertEqual(t-tau*b*x, radial)
                    self.assertEqual(x+tau*a*y, -q*radial)
                    self.assertEqual(1-t-tau*b*y, radial)
                    lower, upper = [], []
                    for k, other in enumerate(full):
                        lo = other.l*t + tau*(other.a-other.b*other.l)*x
                        hi = other.r*t + tau*(other.a-other.b*other.r)*x
                        if lo <= y <= hi:
                            lower.append(k)
                        lo = other.l*(1-t)-tau*(other.a-other.b*other.l)*y
                        hi = other.r*(1-t)-tau*(other.a-other.b*other.r)*y
                        if lo <= x <= hi:
                            upper.append(k)
                    self.assertEqual(len(lower), 1)
                    self.assertEqual(len(upper), 1)
                    self.assertEqual(full[lower[0]], sector)
                    reflected = full[upper[0]]
                    self.assertEqual((reflected.l, reflected.r),
                                     (-sector.r, -sector.l))

    def test_reciprocal_radial_cutoff_rejects_linear_radius_shortcut(self):
        # Unequal endpoint radii: linear radius interpolation would accept
        # this point, but its exact barycentric weights put it outside.
        left_radius, right_radius, q, radial = F(1), F(3), F(1, 2), F(7, 4)
        reciprocal = (1-q)/left_radius + q/right_radius
        cutoff = 1/reciprocal
        linear_shortcut = (1-q)*left_radius + q*right_radius
        self.assertEqual(cutoff, F(3, 2))
        self.assertLessEqual(radial, linear_shortcut)
        self.assertGreater(radial*reciprocal, 1)


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = json.loads(Path(__file__).with_name(
            "structural_final_family.json").read_text())

    def test_complete_original_partition(self):
        intervals = check.read_partition(self.source)
        self.assertEqual(len(intervals), 210)
        self.assertEqual(intervals[0][0], 0)
        self.assertEqual(intervals[-1][1], 1)

    def test_missing_interval_rejected(self):
        data = copy.deepcopy(self.source)
        data["accepted_intervals"].pop(3)
        with self.assertRaises(check.VerificationError):
            check.read_partition(data)

    def test_overlap_rejected(self):
        data = copy.deepcopy(self.source)
        data["accepted_intervals"][1][0] = "0"
        with self.assertRaises(check.VerificationError):
            check.read_partition(data)

    def test_zero_width_interval_rejected(self):
        data = copy.deepcopy(self.source)
        data["accepted_intervals"][0][1] = "0"
        with self.assertRaises(check.VerificationError):
            check.read_partition(data)

    def test_wrong_scope_rejected(self):
        for field, value in [("n", 1024), ("max_tilt", "2"),
                             ("shape", ["0"]), ("q_subdivisions_per_sector", 0),
                             ("normalized_area_lower", "1/2")]:
            data = copy.deepcopy(self.source)
            data[field] = value
            with self.subTest(field=field), self.assertRaises(check.VerificationError):
                check.read_partition(data)

    def test_source_claimed_lower_values_are_not_trusted(self):
        data = copy.deepcopy(self.source)
        for row in data["accepted_intervals"]:
            row[2] = "999999"
        self.assertEqual(check.read_partition(data), check.read_partition(self.source))



    def test_malformed_exact_numbers_rejected(self):
        for value in [0.5, True, None, "1/0", "not-a-rational"]:
            with self.subTest(value=value), self.assertRaises(check.VerificationError):
                check.rational(value)

    def test_malformed_partition_and_endpoint_rejected(self):
        cases = []
        data = copy.deepcopy(self.source)
        data["accepted_intervals"] = []
        cases.append(data)
        data = copy.deepcopy(self.source)
        data["accepted_intervals"][0] = ["0", "1"]
        cases.append(data)
        data = copy.deepcopy(self.source)
        data["accepted_intervals"][-1][1] = "63/64"
        cases.append(data)
        data = copy.deepcopy(self.source)
        del data["shape"]
        cases.append(data)
        data = copy.deepcopy(self.source)
        data["n"] = "2048"
        cases.append(data)
        for data in cases:
            with self.assertRaises(check.VerificationError):
                check.read_partition(data)

if __name__ == "__main__":
    unittest.main()
