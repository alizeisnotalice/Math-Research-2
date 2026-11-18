"""Replay the full n=8 disjoint cut family, with per-child resource receipts."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, localcontext
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import re
import subprocess
import sys

from certify_disjoint import HERE, COUNTER, complement, neighbors, no_bridge, profiles, sha256, weight


def independent_sets_by_row(h):
    """Separate row-transfer count, keeping only the last row's independent mask."""
    last = {0: 1}
    for width in map(int, h):
        last = {s: sum(count for t, count in last.items() if not s & t)
                for s in range(1 << width) if not s & (s << 1)}
    return sum(last.values())


def counter(h, output):
    run = subprocess.run(['/usr/bin/time', '-l', str(COUNTER), '8', h],
                         capture_output=True, text=True, timeout=90, check=True)
    profile = json.loads(run.stdout)
    assert (profile['n'], profile['profile'], profile['vertices']) == (8, h, sum(map(int, h)))
    assert profile['independent_sets'] == independent_sets_by_row(h)
    rss_match = re.search(r'^\s*(\d+)\s+maximum resident set size$', run.stderr, re.M)
    assert rss_match is not None, run.stderr
    (output / f'{h}.json').write_text(run.stdout)
    (output / f'{h}.resource').write_text(run.stderr)
    return h, profile, int(rss_match[1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    parser.add_argument('--workers', type=int, default=3, choices=range(1, 5))
    args = parser.parse_args()
    args.output.mkdir(parents=False, exist_ok=False)
    n = 8
    family = tuple(profiles(n, 0))
    mirrored = {h: complement(h, n) for h in family}
    assert set(mirrored.values()) == set(family)
    adjacent = neighbors(n)
    pairs = list(itertools.combinations(family, 2))
    assert all(no_bridge(n, h, k, adjacent) for h, k in pairs)
    planned = sorted(family, key=independent_sets_by_row, reverse=True)
    assert independent_sets_by_row(planned[0]) == 10_213_044
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        runs = list(pool.map(lambda h: counter(h, args.output), planned))
    regions = {(n, h): record for h, record, _ in runs}
    terms = [weight(n, h, regions) for h in family]
    event_counts = {term['profile']: int(term['count']) for term in terms}
    assert all(event_counts[h] == event_counts[mirrored[h]] > 0 for h in family)
    union_count = sum(event_counts.values())
    factorial = math.factorial(n * n)
    with localcontext() as ctx:
        ctx.prec = 60
        unconditional = str(Decimal(union_count) / Decimal(factorial))
    summary = {'scope': 'Exact n=8 unconditional union probability. No exact A2(8) supplied.',
               'python': sys.version, 'workers': args.workers, 'n': n,
               'counter_sha256': sha256(HERE / 'cut_region.cpp'),
               'binary_sha256': sha256(COUNTER),
               'driver_sha256': sha256(Path(__file__)),
               'utility_sha256': sha256(HERE / 'certify_disjoint.py'),
               'active_row_parity_zero_based': 0, 'events': len(family),
               'verified_no_bridge_pairs': len(pairs), 'children_completed_and_waited': len(runs),
               'max_rss_bytes': max(rss for _, _, rss in runs),
               'children_rss_bytes': {h: rss for h, _, rss in runs},
               'union_count': str(union_count), 'word_count': str(factorial),
               'unconditional_probability_exact': str(Fraction(union_count, factorial)),
               'unconditional_probability_decimal': unconditional, 'terms': terms}
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(union_count, unconditional, summary['max_rss_bytes'])


if __name__ == '__main__':
    main()
