#!/usr/bin/env python3
"""Independent finite step for REPORT.md's fixed-family lower bound.

Standard-library exact arithmetic. Does not import structural_bound.py or trust
the saved interval lower values. Rational Newton upper bounds replace the
source's decimal-grid integer-square-root enclosures.
Alexander Richard Towell, with substantial OpenAI Codex assistance.
Predecessor construction and certificate: Ricardo Pinho; see LICENSE.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
from typing import Any

BITS = 48
NEWTON_STEPS = 8


class VerificationError(ValueError):
    """A proof-relevant condition was not established."""


@dataclass(frozen=True)
class Family:
    n: int
    shape: tuple[F, ...]
    subdivisions: int
    max_tilt: F
    target: F


FINAL_FAMILY = Family(2048, (F(-1, 200), F(27, 1250)), 4, F(1), F(4993, 10000))


@dataclass(frozen=True)
class Sector:
    l: F
    r: F
    a: F
    b: F


def rational(value: Any) -> F:
    if isinstance(value, bool) or not isinstance(value, (str, int)):
        raise VerificationError("Expected an exact rational string or integer")
    try:
        return F(value)
    except (ValueError, ZeroDivisionError) as error:
        raise VerificationError("Invalid rational") from error


def make_sectors(n: int, shape: tuple[F, ...]) -> tuple[Sector, ...]:
    if type(n) is not int or n <= 0 or n % 2:
        raise ValueError("n must be a positive even integer")

    def phi(q: F) -> F:
        return (1-q*q) * (1+sum(c*q**(2*i+2)
                                     for i, c in enumerate(shape))) / 4

    result = []
    for j in range(n//2, n):
        l, r = F(2*j, n)-1, F(2*(j+1), n)-1
        left, right = phi(l), phi(r)
        a, b = (r*left-l*right)/(r-l), (left-right)/(r-l)
        if not (left >= 0 and right >= 0 and 0 < a < 1
                and 0 <= b < 1 and 0 <= a-b <= 1-r):
            raise VerificationError(f"Geometric hypotheses fail in sector {j}")
        result.append(Sector(l, r, a, b))
    return tuple(result)


def sqrt_upper(value: F, bits: int = BITS) -> F:
    """Return a checked dyadic/rational upper bound for sqrt(value), 1<=value<=4."""
    if not 1 <= value <= 4 or type(bits) is not int or bits < 8:
        raise ValueError("Square-root domain is [1,4], with at least 8 bits")
    # Preserve exact rational square roots, in particular the zero-tilt endpoint.
    nr, dr = isqrt(value.numerator), isqrt(value.denominator)
    if nr*nr == value.numerator and dr*dr == value.denominator:
        return F(nr, dr)
    scale = 1 << bits
    upper = F(2)
    for _ in range(NEWTON_STEPS):
        newton = (upper + value/upper) / 2
        upper = F(-((-newton.numerator*scale)//newton.denominator), scale)
    if upper*upper < value:
        raise VerificationError("Square-root upper enclosure failed")
    return upper


def interval_bound(sectors: tuple[Sector, ...], lower: F, upper: F,
                   subdivisions: int, bits: int = BITS) -> F:
    """Lower enclosure for the full normalized area, for s in [lower,upper]."""
    if not 0 <= lower <= upper <= 1:
        raise ValueError("Squared-tilt interval must lie in [0,1]")
    if type(subdivisions) is not int or subdivisions <= 0 or not sectors:
        raise ValueError("Positive subdivision count and nonempty sectors required")
    if type(bits) is not int or bits < 8:
        raise ValueError("At least 8 rounding bits required")
    scale = 1 << bits
    total = F(0)
    for sector in sectors:
        l, r, a, b = sector.l, sector.r, sector.a, sector.b
        a2, ab = a*a, a*b
        root = sqrt_upper(1+upper*(a2+2*b*b), bits)
        numerator = 1+lower*a2
        subtotal = 0
        for k in range(subdivisions):
            q = l+(r-l)*F(k, subdivisions)
            radial = numerator/(2*(1+lower*(a2-ab*q)))
            value = radial*radial/root
            # Each rounded value is <= the integrand on its entire q/s cell.
            subtotal += (value.numerator*scale)//value.denominator
        total += (r-l)*F(subtotal, subdivisions*scale)
    return 2*total


def require_bound(bound: F, target: F) -> None:
    if bound < target:
        raise VerificationError(f"Lower enclosure {bound} does not reach {target}")


def read_partition(data: dict[str, Any]) -> tuple[tuple[F, F], ...]:
    """Validate fixed-family metadata and exact coverage; ignore claimed bounds."""
    try:
        if type(data["n"]) is not int or type(data["q_subdivisions_per_sector"]) is not int:
            raise VerificationError("Discrete parameters must be integers")
        family = Family(data["n"], tuple(rational(c) for c in data["shape"]),
                        data["q_subdivisions_per_sector"], rational(data["max_tilt"]),
                        rational(data["normalized_area_lower"]))
        if family != FINAL_FAMILY:
            raise VerificationError("Certificate is not for the stated final family")
        rows = data["accepted_intervals"]
        if not isinstance(rows, list) or not rows:
            raise VerificationError("Nonempty interval list required")
        intervals = []
        for row in rows:
            if not isinstance(row, list) or len(row) != 3:
                raise VerificationError("Each source interval must be a triple")
            lower, upper = rational(row[0]), rational(row[1])
            if not 0 <= lower < upper <= family.max_tilt**2:
                raise VerificationError("Empty, inverted or out-of-domain interval")
            # row[2] is deliberately not used: recompute the bound independently.
            intervals.append((lower, upper))
        intervals.sort()
        endpoint = F(0)
        for lower, upper in intervals:
            if lower != endpoint:
                raise VerificationError("Gap, overlap or duplicate in tilt partition")
            endpoint = upper
        if endpoint != family.max_tilt**2:
            raise VerificationError("Tilt partition does not reach its endpoint")
        return tuple(intervals)
    except (KeyError, TypeError) as error:
        raise VerificationError("Malformed certificate structure") from error


def verify(data: dict[str, Any]) -> dict[str, Any]:
    intervals = read_partition(data)
    family = FINAL_FAMILY
    sectors = make_sectors(family.n, family.shape)
    bounds = []
    for lower, upper in intervals:
        bound = interval_bound(sectors, lower, upper, family.subdivisions)
        require_bound(bound, family.target)
        bounds.append((lower, upper, bound))
    minimum = min(bound for _, _, bound in bounds)
    return {
        "verified": True,
        "scope": "fixed final secant axes; arbitrary shared even positive node radii; 0<=tau<=1",
        "n": family.n,
        "shape": [str(c) for c in family.shape],
        "max_tilt": str(family.max_tilt),
        "target": str(family.target),
        "sectors_checked": len(sectors),
        "q_subdivisions_per_sector": family.subdivisions,
        "tilt_intervals_checked": len(intervals),
        "q_tilt_cells_checked": len(intervals)*len(sectors)*family.subdivisions,
        "minimum_lower_enclosure": str(minimum),
        "minimum_margin_over_target": str(minimum-family.target),
        "square_root_method": "8 rational Newton steps, each rounded upward on a binary grid",
        "rounding_bits": BITS,
        "interval_bounds": [[str(l), str(u), str(b)] for l, u, b in bounds],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path, nargs="?",
                        default=Path(__file__).with_name("structural_final_family.json"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    raw = args.certificate.read_bytes()
    result = verify(json.loads(raw))
    result["source_certificate_sha256"] = hashlib.sha256(raw).hexdigest()
    result["checker_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k: v for k, v in result.items() if k != "interval_bounds"}, indent=2))


if __name__ == "__main__":
    main()
