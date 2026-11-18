"""Independently verify the original radian 11D root and its tree subdivision."""
import argparse
import json
import platform
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import isqrt
from pathlib import Path

import numpy as np


from rational_trig import PI, scale, trig

TARGET = F(836600003,10**9)
SOURCE_TAU = F(4183,5000)
PARTS = 2**11
BINARY_FIELDS = ("lo", "hi", "T")
IDENTITY_FIELDS = ("ns", "tau", "delta", "nparts", "split0", "N", "eps")
PRECISION = 10**21


def sqrt_bounds(x):
    radicand = x * PRECISION**2
    k = isqrt(radicand.numerator // radicand.denominator)
    assert k**2 <= radicand < (k+1)**2
    return F(k, PRECISION), F(k+1, PRECISION)


def lower_far_area(n, width):
    radius_lower = 1 / (2*trig(scale(PI,F(1,2*n)), "cos")[1])
    cos_lower = trig(scale(PI,F(1,n)), "cos")[0]
    s_floor = sqrt_bounds(radius_lower**2 + width**2 + 2*radius_lower*width*cos_lower)[0]
    z_lo, z_hi = (2*x for x in sqrt_bounds(s_floor**2 - F(1,4)))
    t = z_hi/(1 + 2*s_floor)
    assert 0 < t < F(1,2)
    # atan(z)=2*atan(z/(1+sqrt(1+z*z))) and sqrt(1+z*z)=2*s_floor.
    # The odd partial sum gives a strict upper bound on atan(t).
    terms = 25
    atan_upper = sum(((-1)**k * t**(2*k+1) / (2*k+1) for k in range(terms)), F(0))
    return (PI[0] + z_lo - 2*atan_upper)/4


def binary_values(meta):
    return {name: tuple(map(F, meta[name])) for name in BINARY_FIELDS}


def decimal_lower(value, digits=20):
    unit = 10**digits
    q = (value*unit).__floor__()
    whole, fraction = divmod(q,unit)
    return f"{whole}.{fraction:0{digits}d}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certdir", type=Path)
    args = parser.parse_args()
    headers = sorted(args.certdir.glob("part[0-9][0-9][0-9][0-9].json"))
    assert len(headers) == PARTS
    root = json.loads(headers[0].read_bytes())
    exact_root = binary_values(root)
    lo, hi, widths = (exact_root[k] for k in BINARY_FIELDS)
    assert len(lo) == len(hi) == 11 and root["split0"] == [2]*11
    assert root["ns"] == [3,5,7,9] and root["nparts"] == PARTS
    assert F(str(root["tau"])) == SOURCE_TAU
    mids = tuple((a+b)/2 for a, b in zip(lo, hi, strict=True))
    assert tuple(map(F, map(float, mids))) == mids
    assert all(F(float(hi[k])-float(lo[k])) == hi[k]-lo[k] for k in range(11))
    for k in range(11):
        # Reproduce the actual floating order of operations in the source decoder.
        a, b = float(lo[k]), float(hi[k])
        assert F(a + (b-a)/2) == mids[k]
        assert F(a + (b-a)) == hi[k]

    angle_indices = (2,5,8)
    angle_factors = (F(1,5),F(2,7),F(2,9))
    angle_gaps = {}
    for k, factor, body in zip(angle_indices, angle_factors, root["ns"][1:], strict=True):
        gap = hi[k]-factor*PI[1]
        assert lo[k] == 0 and gap > 0
        angle_gaps[body] = decimal_lower(gap)

    exterior = {}
    for k, (n, width) in enumerate(zip(root["ns"], widths, strict=True)):
        j = 3*k
        assert lo[j] == lo[j+1] == -width < 0 < width == hi[j] == hi[j+1]
        area = lower_far_area(n,width)
        assert area > TARGET
        exterior[n] = {
            "width": str(width),
            "strict_lower_area": decimal_lower(area),
            "strict_lower_margin_over_source_tau": decimal_lower(area-SOURCE_TAU),
            "strict_lower_margin_over_target": decimal_lower(area-TARGET),
        }

    root_sha = sha256()
    observed = set()
    trees = Counter()
    for part, path in enumerate(headers):
        assert path.name == f"part{part:04}.json"
        raw = path.read_bytes()
        root_sha.update(raw)
        meta = json.loads(raw)
        assert meta["part"] == part and binary_values(meta) == exact_root
        assert all(meta[field] == root[field] for field in IDENTITY_FIELDS)
        # Independent lexicographic binary coordinate extraction, last coordinate fastest.
        cell = tuple((part >> (10-k)) & 1 for k in range(11))
        assert cell == tuple(map(int,np.unravel_index(part,(2,)*11)))
        assert cell not in observed
        observed.add(cell)

        codes = np.load(path.with_suffix(".npy"), allow_pickle=False)
        assert codes.dtype == np.uint8 and codes.shape == (meta["boxes"],)
        branch = codes < 11
        assert np.all(branch | (codes >= 253))
        pending = 1 + np.cumsum(2*branch.astype(np.int64)-1)
        assert np.all(pending[:-1] > 0) and pending[-1] == 0
        trees["nodes"] += len(codes)
        for code, name in ((253,"multi"),(254,"bnd"),(255,"far")):
            observed_count = int(np.count_nonzero(codes == code))
            assert observed_count == meta[name]
            trees[name] += observed_count
    assert observed == set(product(range(2),repeat=11))
    assert trees["nodes"] == 2*sum(trees[name] for name in ("multi","bnd","far")) - PARTS
    print(json.dumps({
        "python": platform.python_version(),
        "numpy": np.__version__,
        "source_tau": str(SOURCE_TAU),
        "target_tau": str(TARGET),
        "header_concat_sha256": root_sha.hexdigest(),
        "angle_upper_padding_strict_lower_radians": angle_gaps,
        "exterior": exterior,
        "root_parts_checked": PARTS,
        "exact_root_grid_cover": True,
        "raw_preorder_structures": dict(trees),
        "scope": "Every consistent raw root header and complete preorder shape; independent rational root geometry. Leaf area arithmetic not rerun.",
    }, indent=2))


if __name__ == "__main__":
    main()
