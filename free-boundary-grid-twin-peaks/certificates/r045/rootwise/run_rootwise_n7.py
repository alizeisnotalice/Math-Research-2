"""Exact root-pair stratification on the 1-Lipschitz widths 2..5, n=7."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import combinations
import json
import math
from pathlib import Path
import re
import subprocess
import sys

from run_parity_n8 import HERE, COUNTER, complement, digest, geometry, independent_sets, nonintersection, side_region

N = 7
TOTAL = N * N
REPO = HERE.parents[2]
NORMALIZER = REPO / 'computation/results/n7-exact-r005.json'
OLD_RAW = HERE.parents[1] / 'r043/n7-curves.json'


def curves(partial=''):
    if len(partial) == N:
        yield partial
        return
    for width in range(2, 6):
        if not partial or abs(width - int(partial[-1])) <= 1:
            yield from curves(partial + str(width))


def reflected_y(v):
    return (N - 1 - v // N) * N + v % N


def reflected_x(v):
    return v // N * N + N - 1 - v % N


def turned(v):
    return v % N * N + v // N


def read(path):
    return json.loads(path.read_text())


def timed_region(h, out):
    process = subprocess.run(['/usr/bin/time', '-l', str(COUNTER), str(N), h],
                             capture_output=True, text=True, check=True, timeout=90)
    data = json.loads(process.stdout)
    assert data['n'] == N and data['profile'] == h
    assert data['independent_sets'] == independent_sets(h)
    assert data['vertices'] == sum(map(int, h)) <= 36
    values = list(map(int, data['root_counts']))
    assert len(values) == TOTAL and sum(values) == int(data['rooted_count'])
    H = side_region(h)
    eligible = {v for v in H if v % N < (N + 3) // 4 and ADJ[v] <= H}
    assert eligible == {v for v, value in enumerate(values) if value > 0}
    rss = re.search(r'^\s*(\d+)\s+maximum resident set size$', process.stderr, re.M)
    assert rss, process.stderr
    (out / f'{h}.json').write_text(process.stdout)
    (out / f'{h}.resource').write_text(process.stderr)
    return h, values, int(rss[1])


ADJ = geometry(N)


def restore(physical, h):
    if h <= h[::-1]:
        return physical[h]
    return [physical[h[::-1]][reflected_y(v)] for v in range(TOTAL)]


def aggregate(profiles, physical):
    assigned = defaultdict(list)
    counts_by_class = defaultdict(lambda: defaultdict(int))
    event_counts = {}
    for h in profiles:
        parity = ''.join(str(int(w) % 2) for w in h)
        assigned[parity].append(h)
        left = restore(physical, h)
        right = restore(physical, complement(h))
        shuffle = math.comb(TOTAL, sum(map(int, h)))
        for p, a in enumerate(left):
            if not a:
                continue
            for q_mirror, b in enumerate(right):
                if not b:
                    continue
                q = reflected_x(q_mirror)
                pair = tuple(sorted((p, q)))
                rotated_pair = tuple(sorted((turned(p), turned(q))))
                counts_by_class[parity, 0][pair] += shuffle * a * b
                counts_by_class[parity, 1][rotated_pair] += shuffle * a * b
        event_counts[h] = shuffle * sum(left) * sum(right)
    assert set(profiles) == set(map(complement, profiles))
    assert all(event_counts[h] == event_counts[complement(h)] == event_counts[h[::-1]] for h in profiles)
    side = {h: side_region(h) for h in profiles}
    bridge_checks = 0
    for class_curves in assigned.values():
        for h, k in combinations(class_curves, 2):
            assert nonintersection(side[h], side[k], ADJ)
            bridge_checks += 1
    old_files = read(OLD_RAW)
    assert len(old_files) == 128
    for old in old_files:
        assert int(old['count']) == event_counts[old['curve']]
        assert int(old['left_rooted_count']) == sum(restore(physical, old['curve']))
        assert int(old['right_rooted_count']) == sum(restore(physical, complement(old['curve'])))
    all_pairs = set().union(*counts_by_class.values())
    rows = []
    for pair in sorted(all_pairs):
        strata = {key: counts[pair] for key, counts in counts_by_class.items() if pair in counts}
        best, peak_count = max(strata.items(), key=lambda item: (item[1], item[0]))
        assert peak_count > 0 and abs(pair[0] % N - pair[1] % N) + abs(pair[0] // N - pair[1] // N) >= 3
        rows.append({'peaks': list(pair), 'parity': best[0], 'transposed': bool(best[1]), 'count': str(peak_count)})
    single = {p: sum(weights.values()) for (p, t), weights in counts_by_class.items() if not t}
    A2 = int(read(NORMALIZER)['A2'])
    mass = sum(int(row['count']) for row in rows)
    with localcontext() as ctx:
        ctx.prec = 60
        lower = str(Decimal(mass) / Decimal(A2))
        best_class = str(Decimal(max(single.values())) / Decimal(A2))
    return {'profiles_per_orientation': len(event_counts), 'class_count_per_orientation': len(single),
            'verified_no_bridge_pairs_per_orientation': bridge_checks, 'A2': str(A2),
            'A2_source': str(NORMALIZER), 'A2_source_sha256': digest(NORMALIZER),
            'old_central_curve_crosschecks': len(old_files),
            'best_single_class_count': str(max(single.values())), 'best_single_class_probability': best_class,
            'rootwise_disjoint_count': str(mass), 'rootwise_disjoint_probability_exact': str(Fraction(mass, A2)),
            'rootwise_disjoint_probability_decimal': lower, 'root_pairs': len(rows), 'root_strata': rows,
            'single_class_counts': {p: str(count) for p, count in single.items()},
            'event_counts': {h: str(count) for h, count in event_counts.items()}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    parser.add_argument('--workers', type=int, default=12, choices=range(1, 13))
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    profiles = list(curves())
    canonical = sorted({min(h, h[::-1]) for h in profiles}, key=lambda h: (-independent_sets(h), h))
    assert len(profiles) == 1220 and len(canonical) == 644
    assert max(independent_sets(h) for h in profiles) == 3_729_091
    assert all(max(sum(map(int, h)), TOTAL - sum(map(int, h))) <= 36 for h in profiles)
    raw = args.output / 'raw'
    raw.mkdir()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        measured = list(pool.map(lambda h: timed_region(h, raw), canonical))
    physical = {h: values for h, values, _ in measured}
    data = aggregate(profiles, physical)
    data.update(python=sys.version, workers=args.workers, n=N,
                source_sha256=digest(HERE / 'cut_region_rootwise.cpp'), driver_sha256=digest(Path(__file__)),
                regions_waited=len(measured), max_rss_bytes=max(rss for _, _, rss in measured),
                regions_rss_bytes={h: rss for h, _, rss in measured})
    (args.output / 'summary.json').write_text(json.dumps(data, indent=2) + '\n')
    print('regions', data['regions_waited'], 'max rss', data['max_rss_bytes'], 'profiles', data['profiles_per_orientation'])
    print('best class', data['best_single_class_probability'])
    print('strict conditional lower', data['rootwise_disjoint_count'], data['rootwise_disjoint_probability_decimal'])


if __name__ == '__main__':
    main()
