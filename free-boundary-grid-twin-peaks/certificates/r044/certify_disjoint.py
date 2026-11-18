"""Verify and count the two alternating-row disjoint cut families."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
COUNTER = Path(os.environ.get('A022_R044_REGION_BINARY', HERE / 'cut_region')).resolve()
PROJECT = HERE.parents[1]
OLD_CURVES = HERE.parent / 'r043/n7-curves.json'
OLD_TOTALS = HERE / 'inputs'


def read(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def neighbors(n):
    vertices = range(n * n)
    return [{u for u in vertices if abs(v % n - u % n) + abs(v // n - u // n) == 1}
            for v in vertices]


def left_region(n, h):
    return {v for v in range(n * n) if v % n < int(h[v // n])}


def complement(h, n):
    return ''.join(str(n - int(x)) for x in h)


def profiles(n, active_parity):
    m = n // 2
    active = [y for y in range(n) if y % 2 == active_parity]
    for choices in itertools.product((m - 1, m + 1), repeat=len(active)):
        h = [m] * n
        for y, width in zip(active, choices):
            h[y] = width
        yield ''.join(map(str, h))


def no_bridge(n, h, k, adjacent):
    """Apply the round-42 shared-side bridge condition directly to grid sets."""
    left_h, left_k = left_region(n, h), left_region(n, k)
    common_left = left_h & left_k
    common_right = set(range(n * n)) - (left_h | left_k)
    for missing in (left_h - left_k, left_k - left_h):
        while missing:
            frontier = {missing.pop()}
            component = set(frontier)
            while frontier:
                reached = adjacent[frontier.pop()] & missing
                frontier.update(reached)
                component.update(reached)
                missing.difference_update(reached)
            if not any(adjacent[v] & common_left and adjacent[v] & common_right for v in component):
                return True
    return False


def region_counter(key):
    n, h = key
    program = [str(COUNTER), str(n), h]
    result = subprocess.run(program, capture_output=True, text=True, check=True, timeout=30)
    record = json.loads(result.stdout)
    assert result.stderr == '' and (record['n'], record['profile']) == key
    assert record['vertices'] == sum(map(int, h))
    return key, record


def compute_regions(keys, workers):
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return dict(pool.map(region_counter, sorted(keys)))


def weight(n, h, region_ledger):
    left = region_ledger[n, h]
    right = region_ledger[n, complement(h, n)]
    assert left['vertices'] + right['vertices'] == n * n
    value = math.comb(n * n, left['vertices']) * int(left['rooted_count']) * int(right['rooted_count'])
    return {'profile': h, 'left': left, 'right': right, 'count': str(value)}


def audit_old_n7(regions, old_files):
    for record in old_files:
        h = record['curve']
        fresh = weight(7, h, regions)
        assert fresh['count'] == record['count']
        for side in ('left', 'right'):
            for field in ('vertices', 'allowed_vertices', 'rooted_count', 'independent_sets', 'alpha'):
                assert fresh[side][field] == record[f'{side}_{field}'], (side, field, h)


def audit_n5(n5_results, ledger):
    records = {job['curves'][0]: int(job['count']) for job in ledger['jobs']
               if job['curves'][0] == job['curves'][1]}
    for record in n5_results:
        h = record['profile']
        assert h in ledger['curves']
        hc = complement(h, 5)
        key = min(h, hc, h[::-1], hc[::-1])
        assert int(record['count']) == records[key], (h, int(record['count']), records[key])


def summarize(n, h_family, region_ledger):
    adjacent = neighbors(n)
    pairs = list(itertools.combinations(h_family, 2))
    assert all(no_bridge(n, h, k, adjacent) for h, k in pairs)
    terms = [weight(n, h, region_ledger) for h in h_family]
    assert all(int(term['count']) > 0 for term in terms)
    evidence = OLD_TOTALS / f'parallelhash-generic-n{n}.json' if n < 7 else PROJECT / 'computation/results/n7-exact-r005.json'
    original = int(read(evidence)['A2'])
    total = sum(int(record['count']) for record in terms)
    with localcontext() as ctx:
        ctx.prec = 48
        ratio = str(Decimal(total) / Decimal(original))
    return {'n': n, 'events': len(terms), 'verified_no_bridge_pairs': len(pairs),
            'A2_source': str(evidence), 'A2_sha256': sha256(evidence), 'A2': str(original),
            'union_count': str(total), 'union_probability_exact': str(Fraction(total, original)),
            'union_probability_decimal': ratio, 'terms': terms}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=3)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    all_families = {(n, parity): tuple(profiles(n, parity)) for n in (4, 5, 6, 7) for parity in (0, 1)}
    new_keys = {(n, x) for (n, _), family in all_families.items() for h in family
                for x in (h, complement(h, n))}
    old_files = read(OLD_CURVES)
    old_keys = {(7, row['curve']) for row in old_files}
    assert len(old_keys) == 128
    regions = compute_regions(new_keys | old_keys, args.workers)
    audit_old_n7(regions, old_files)
    bundled = [dict(active_row_parity_zero_based=parity, **summarize(n, family, regions))
               for (n, parity), family in all_families.items()]
    n5_ledger_path = PROJECT / 'certificates/r042/n5-lipschitz-ledger.json'
    n5_results = [term for record in bundled if record['n'] == 5 for term in record['terms']]
    audit_n5(n5_results, read(n5_ledger_path))
    evidence = {'scope': 'Exact finite counts for disjoint families. No asymptotic extrapolation.',
                'python': sys.version, 'workers': args.workers,
                'region_counter_sha256': sha256(HERE / 'cut_region.cpp'),
                'region_calls_completed_and_waited': len(regions),
                'n5_crosscheck_terms': len(n5_results), 'n5_ledger_sha256': sha256(n5_ledger_path),
                'old_n7_crosscheck_terms': len(old_files), 'families': bundled}
    with args.out.open('x') as fp:
        fp.write(json.dumps(evidence, indent=2) + '\n')
    for record in bundled:
        print(record['n'], record['active_row_parity_zero_based'], record['events'],
              record['union_count'], record['union_probability_decimal'])


if __name__ == '__main__':
    main()
