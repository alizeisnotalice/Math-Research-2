#!/usr/bin/env python3
"""Independent fixed-total sub-envelope; stdlib Fraction/isqrt, no producer imports."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import isqrt
from pathlib import Path
import platform
import sys
import time


def fs(x):
    return str(x)


def sqrt_bounds(q, bits):
    den = 1 << bits
    k = isqrt(q.numerator * den * den // q.denominator)
    assert F(k, den) ** 2 <= q < F(k + 1, den) ** 2
    return F(k, den), F(k + 1, den)


def clipped_length(lo, hi, center_lo, center_hi):
    return max(F(0), min(hi, center_hi) - max(lo, center_lo))


def case(n, bits):
    assert n >= 3 and n % 2 == 1
    w = 1 - F(1, n * n)
    points = [F(i, n - 1) * w for i in range(n)]
    ends = sorted({p + shift for p in points for shift in (F(-1, 2), F(1, 2))})
    cells = []
    for lo, hi in zip(ends, ends[1:]):
        mid = (lo + hi) / 2
        count = sum(abs(mid - p) < F(1, 2) for p in points)
        assert count > 0
        cells.append((lo, hi, count))

    @lru_cache(None)
    def cutoff(count):
        lower, upper = sqrt_bounds(2 - F(count, n * n), bits)
        return lower - F(1, 2), upper - F(1, 2)

    lower_area = upper_area = F(0)
    proper_cells = full_cells = 0
    for xl, xh, cx in cells:
        for yl, yh, cy in cells:
            count = cx * cy
            if count == n * n:
                full_cells += 1
                continue
            proper_cells += 1
            hlo, hhi = cutoff(count)
            area = (xh - xl) * (yh - yl)
            # Failure c_S<=h is K=[w-h,h]^2. Larger h removes more area.
            small_failure = clipped_length(xl, xh, w - hlo, hlo) * clipped_length(yl, yh, w - hlo, hlo)
            large_failure = clipped_length(xl, xh, w - hhi, hhi) * clipped_length(yl, yh, w - hhi, hhi)
            lower_area += area - large_failure
            upper_area += area - small_failure
    assert 0 <= lower_area <= upper_area <= (1 + w) ** 2

    # Actual positive finite source-domain embedding: q=(33/32,w/2), mass epsilon.
    eps, qx = F(1, 1000), F(33, 32)
    # Any q-containing subset of mass>=1 must have its leftmost grid column
    # <=eps*n*w/(n-1)<=2eps. Its span then rules out eligibility.
    max_left = eps * n * w / (n - 1)
    assert max_left <= 2 * eps
    max_mass_exclusion_margin = (qx - 2 * eps) ** 2 - (1 + eps)
    bridge_margin = eps + F(1, n * n) - (qx - w) ** 2
    nonterminal_margin = qx ** 2 - F(35, 33) * (1 + eps)
    mst_cost = (n + 1) * w
    assert max_mass_exclusion_margin > 0 and bridge_margin > 0 and nonterminal_margin > 0
    assert w > F(3, 4) and mst_cost > F(5, 2)
    return {
        'n': n, 'old_sources': n * n, 'width': fs(w),
        'one_dimensional_cells': len(cells), 'proper_mass_cells': proper_cells,
        'full_mass_cells_excluded': full_cells, 'bits': bits,
        'area_lower': fs(lower_area), 'area_upper': fs(upper_area),
        'area_lower_float_display_only': float(lower_area),
        'area_upper_float_display_only': float(upper_area),
        'interval_width': fs(upper_area - lower_area),
        'strictly_above_7_over_2': lower_area > F(7, 2),
        'source_domain': {
            'unique_maximum_S_mass': '1', 'outside_source': [fs(qx), fs(w / 2), fs(eps)],
            'maximum_exclusion_margin': fs(max_mass_exclusion_margin),
            'bridge_margin': fs(bridge_margin), 'nonterminal_margin': fs(nonterminal_margin),
            'support_tree_cost': fs(mst_cost),
            'verification': 'column-count general proof plus exact scalar margins; not all 2^(n^2+1) subsets enumerated'
        }
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    started = time.perf_counter()
    spec = json.loads(args.input.read_text())
    result = {
        'schema': 'r070-independent-grid-fixed-total-sub-envelope-v1',
        'algorithm': 'rational rectangular cells; A=count/n^2; cutoff sqrt(2-A)-1/2; subtract central-square intersection',
        'scope': 'sub-envelope contained in full internal-mass/r_U clamp; not actual mixed halo or deletion area',
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
        'python': sys.version, 'platform': platform.platform(),
        'cases': [case(n, spec['bits']) for n in spec['odd_grid_sizes']],
        'elapsed_seconds': time.perf_counter() - started,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'cases': [{'n': c['n'], 'area_lower': c['area_lower_float_display_only'], 'above3.5': c['strictly_above_7_over_2']} for c in result['cases']], 'seconds': result['elapsed_seconds']}))


if __name__ == '__main__':
    main()
