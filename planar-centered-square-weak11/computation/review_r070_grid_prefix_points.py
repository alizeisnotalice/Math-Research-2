#!/usr/bin/env python3
"""Small independent replay: list every source; group actual 2D distances."""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import platform
import time


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--saved', type=Path, required=True)
    ap.add_argument('--producer', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--n', type=int, default=17)
    args = ap.parse_args()
    start = time.perf_counter()
    spec = json.loads(args.input.read_text())
    saved = json.loads(args.saved.read_text())
    assert saved['source_sha256'] == sha(args.producer) and saved['input_sha256'] == sha(args.input)
    n = args.n
    assert n in spec['grid_sizes'] and n >= 3
    k = isqrt(n)
    radius, mass, oldmass = F(k - 1, 2 * k), F((k - 1) ** 2, k * k), F(1, n * n)
    axis = [F(j * (n + 1), n * n) for j in range(n)]
    points = [(x, y) for x in axis for y in axis]
    old = {(F(row['x']), F(row['y'])): row for chunk in saved['chunks'] if chunk['n'] == n for row in chunk['rows']}
    rows, total_events = [], 0
    for nx in spec['coordinate_numerators']:
        for ny in spec['coordinate_numerators']:
            x, y = F(nx, 16), F(ny, 16)
            distances = [max(abs(x - px), abs(y - py)) for px, py in points]
            U = sum(d < radius for d in distances)
            remaining = Counter(d for d in distances if d >= radius)
            prefix, active, failed, minimum, first = 0, 0, 0, None, None
            for distance, count in sorted(remaining.items()):
                total_events += 1
                prefix += count
                M = oldmass * prefix
                if M + mass <= 1:
                    continue
                active += 1
                slack = (distance + radius) ** 2 - (M + mass)
                minimum = slack if minimum is None else min(minimum, slack)
                if slack <= 0:
                    failed += 1
                    if first is None:
                        first = dict(distance=str(distance), prefix_mass=str(M), squared_slack=str(slack))
            feasible = 0 < U < n * n and U * oldmass < mass
            row = dict(x=str(x), y=str(y), U_count=U, candidate_mass_exceeds_U=feasible,
                       active_prefix_count=active, violated_prefix_count=failed,
                       all_prefix_pass=feasible and failed == 0,
                       min_active_squared_slack=None if minimum is None else str(minimum),
                       first_violation=first)
            assert row == old[x, y], (x, y, row, old[x, y])
            rows.append(row)
    result = dict(status='PASS_SMALL_LITERAL_2D_SOURCE_PREFIX_REPLAY', n=n, sources=len(points),
                  candidate_mass=str(mass), receiver_points=len(rows),
                  all_prefix_pass_count=sum(row['all_prefix_pass'] for row in rows),
                  grouped_actual_distance_events=total_events, rows=rows,
                  source_sha256=sha(Path(__file__)), input_sha256=sha(args.input),
                  producer_sha256=sha(args.producer), saved_result_sha256=sha(args.saved),
                  python=platform.python_version(), elapsed_seconds=time.perf_counter() - start,
                  scope='Only requested n fully replayed by literal n^2 sources and actual 2D L-infinity distance groups; root larger n results source-bound only. Finite receiver points, not area or actual outside realizability.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ('status', 'n', 'sources', 'receiver_points', 'all_prefix_pass_count', 'grouped_actual_distance_events', 'elapsed_seconds')}))


if __name__ == '__main__':
    main()
