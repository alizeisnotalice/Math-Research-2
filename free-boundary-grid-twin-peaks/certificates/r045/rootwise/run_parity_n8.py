"""Plan and count the full area-capped alternating parity batches on n=8."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import combinations
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys


HERE = Path(__file__).resolve().parent
COUNTER = Path(os.environ.get("A022_R045_REGION_BINARY", HERE / "cut_region_rootwise")).resolve()
N = 8
AREA_CAP = 36
SET_CAP = 12_000_000


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def alternating_curves(n, initial_parity, partial=''):
    if len(partial) == n:
        yield partial
        return
    parity = (initial_parity + len(partial)) % 2
    for k in range(1, n):
        if k % 2 == parity and (not partial or abs(k - int(partial[-1])) == 1):
            yield from alternating_curves(n, initial_parity, partial + str(k))


def complement(h):
    return ''.join(str(len(h) - int(x)) for x in h)


@lru_cache(maxsize=None)
def independent_sets(h):
    last_row = {0: 1}
    for width in map(int, h):
        last_row = {s: sum(count for t, count in last_row.items() if not s & t)
                    for s in range(1 << width) if not s & (s << 1)}
    return sum(last_row.values())


def get_plan():
    all_profiles = tuple(alternating_curves(N, 0))
    assert len(all_profiles) == 232
    eligible = [h for h in all_profiles if max(sum(map(int, h)), N * N - sum(map(int, h))) <= AREA_CAP]
    assert len(eligible) == 78 and set(map(complement, eligible)) == set(eligible)
    assert max(independent_sets(h) for h in eligible) == 10_610_055
    assert all(independent_sets(h) <= SET_CAP for h in eligible)
    assert set(map(lambda h: h[::-1], eligible)) == {
        h for h in alternating_curves(N, 1)
        if max(sum(map(int, h)), N * N - sum(map(int, h))) <= AREA_CAP}
    return sorted(eligible, key=lambda h: (-independent_sets(h), h))


def geometry(n):
    adj = []
    for v in range(n * n):
        x, y = v % n, v // n
        adj.append({u for u in range(n * n) if abs(u % n - x) + abs(u // n - y) == 1})
    return adj


def side_region(h):
    n = len(h)
    return {v for v in range(n * n) if v % n < int(h[v // n])}


def nonintersection(left, other_left, adjacent):
    common_left = left & other_left
    common_right = set(range(len(adjacent))) - (left | other_left)
    for diff in (left - other_left, other_left - left):
        while diff:
            pending = {diff.pop()}
            has_bridge = False
            while pending:
                v = pending.pop()
                has_bridge |= bool(adjacent[v] & common_left and adjacent[v] & common_right)
                reached = adjacent[v] & diff
                diff.difference_update(reached)
                pending.update(reached)
            if not has_bridge:
                return True
    return False


def audit_family(family):
    adjacency = geometry(N)
    sides = {h: side_region(h) for h in family}
    assert all(nonintersection(sides[h], sides[k], adjacency) for h, k in combinations(family, 2))
    return len(family) * (len(family) - 1) // 2


def timed_region(h, out):
    command = ['/usr/bin/time', '-l', str(COUNTER), str(N), h]
    process = subprocess.run(command, capture_output=True, text=True, check=True, timeout=120)
    record = json.loads(process.stdout)
    assert record['n'] == N and record['profile'] == h
    assert record['independent_sets'] == independent_sets(h)
    assert record['vertices'] == sum(map(int, h)) <= AREA_CAP
    roots = list(map(int, record['root_counts']))
    h_side = side_region(h)
    neighbors = geometry(N)
    eligible = {v for v in h_side if v % N < (N + 3) // 4 and neighbors[v] <= h_side}
    assert len(roots) == N * N and sum(roots) == int(record['rooted_count'])
    assert {v for v, a in enumerate(roots) if a > 0} == eligible
    rss = re.search(r'^\s*(\d+)\s+maximum resident set size$', process.stderr, re.M)
    assert rss, process.stderr
    (out / f'{h}.json').write_text(process.stdout)
    (out / f'{h}.resource').write_text(process.stderr)
    return h, record, int(rss[1])


def reflect_x(v):
    return v // N * N + N - 1 - v % N


def reflect_y(v):
    return (N - 1 - v // N) * N + v % N


def transpose(v):
    return v % N * N + v // N


def root_pair_counts(plan, regions):
    phase_zero = {}
    profiles = []
    for h in plan:
        left = list(map(int, regions[h]['root_counts']))
        right = list(map(int, regions[complement(h)]['root_counts']))
        shuffle = math.comb(N * N, regions[h]['vertices'])
        count = 0
        for p, a in enumerate(left):
            if not a:
                continue
            for q_in_mirror, b in enumerate(right):
                if not b:
                    continue
                q = reflect_x(q_in_mirror)
                pair = tuple(sorted((p, q)))
                phase_zero[pair] = phase_zero.get(pair, 0) + shuffle * a * b
                count += shuffle * a * b
        assert count == shuffle * int(regions[h]['rooted_count']) * int(regions[complement(h)]['rooted_count']) > 0
        profiles.append({'profile': h, 'reflected_profile': complement(h), 'count': str(count)})
    phase_one = {tuple(sorted(map(reflect_y, pair))): count for pair, count in phase_zero.items()}
    classes = [phase_zero, phase_one]
    classes += [{tuple(sorted(map(transpose, pair))): count for pair, count in phase.items()}
                for phase in classes[:]]
    stratum = set().union(*classes)
    rows = []
    for pair in sorted(stratum):
        weights = [phase.get(pair, 0) for phase in classes]
        assert max(weights) > 0
        best = max(range(len(weights)), key=lambda i: weights[i])
        rows.append({'peaks': list(pair), 'class_counts': list(map(str, weights)),
                     'chosen_class': best, 'chosen_count': str(weights[best])})
    total = sum(int(row['count']) for row in profiles)
    assert all(sum(phase.values()) == total for phase in classes)
    mass = sum(int(row['chosen_count']) for row in rows)
    assert total <= mass <= 4 * total
    return total, mass, profiles, rows


def old_family_gate(regions, old_path=None):
    m = N // 2
    # The old production used widths 3/5 on y even, whereas our independent
    # plan uses even width on y even. Row reflection identifies the two.
    family = [h for h in regions if all(int(w) == m if y % 2 == 0 else int(w) in (m - 1, m + 1)
                                        for y, w in enumerate(h))]
    assert len(family) == 16
    subtotal, _, terms, _ = root_pair_counts(family, regions)
    assert subtotal == 15497852066351452649557374277633220714870471956684105383340393957380983808
    if old_path is not None:
        old = json.loads(old_path.read_text())
        assert int(old['union_count']) == subtotal
        old_counts = {row['profile']: row['count'] for row in old['terms']}
        assert {row['profile'][::-1]: row['count'] for row in terms} == old_counts
    return subtotal


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    parser.add_argument('--workers', type=int, default=3, choices=range(1, 5))
    parser.add_argument('--plan-only', action='store_true')
    parser.add_argument('--old-summary', type=Path)
    options = parser.parse_args()
    plan = get_plan()
    options.output.mkdir(exist_ok=False)
    audit_pairs = audit_family(plan)
    max_by_side = max(independent_sets(h) for h in plan)
    plan_rows = [{'profile': h, 'area': sum(map(int, h)), 'independent_sets': independent_sets(h)} for h in plan]
    metadata = {'n': N, 'area_cap': AREA_CAP, 'independent_set_cap': SET_CAP,
                'retained_profiles_each_phase': len(plan), 'max_region_independent_sets': max_by_side,
                'checked_no_bridge_pairs_base_phase': audit_pairs, 'rootwise_limb_bits': 192,
                'physical_region_jobs': plan_rows, 'python': sys.version, 'workers': options.workers,
                'driver_sha256': digest(Path(__file__)),
                'counter_source_sha256': digest(HERE / 'cut_region_rootwise.cpp')}
    (options.output / 'plan.json').write_text(json.dumps(metadata, indent=2) + '\n')
    if options.plan_only:
        print('profiles per phase', len(plan), 'max independent sets', max_by_side, 'bridge exclusions', audit_pairs)
        return
    with ThreadPoolExecutor(max_workers=options.workers) as pool:
        measured = list(pool.map(lambda h: timed_region(h, options.output), plan))
    regions = {h: data for h, data, _ in measured}
    total, rootwise, event_rows, root_rows = root_pair_counts(plan, regions)
    gate_total = old_family_gate(regions, options.old_summary)
    normalizer = math.factorial(N * N)
    with localcontext() as context:
        context.prec = 60
        single_probability = str(Decimal(total) / Decimal(normalizer))
        rootwise_probability = str(Decimal(rootwise) / Decimal(normalizer))
    summary = dict(metadata, scope='Exact unconditioned finite counts. A2(8) has not been supplied.',
                   children_waited=len(measured), max_rss_bytes=max(x for _, _, x in measured),
                   regions_rss_bytes={h: rss for h, _, rss in measured},
                   counter_binary_sha256=digest(COUNTER),
                   old_family_gate_count=str(gate_total),
                   old_summary_sha256=digest(options.old_summary) if options.old_summary else None,
                   single_orientation_phase_union_count=str(total),
                   rootwise_disjoint_subunion_count=str(rootwise),
                   unconditional_single_probability_exact=str(Fraction(total, normalizer)),
                   unconditional_rootwise_probability_exact=str(Fraction(rootwise, normalizer)),
                   unconditional_single_probability_decimal=single_probability,
                   unconditional_rootwise_probability_decimal=rootwise_probability,
                   event_counts=event_rows, root_strata=root_rows)
    (options.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('union count', total, single_probability)
    print('rootwise count', rootwise, rootwise_probability)
    print('max RSS', summary['max_rss_bytes'], 'children waited', len(measured))


if __name__ == '__main__':
    main()
