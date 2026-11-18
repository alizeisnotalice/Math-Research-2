"""Exact finite A2 upper bounds from protected peaks in complementary regions."""
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import permutations
from math import comb, factorial
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess
import sys

from reference_crt import crt, independent_sizes, verified_primes

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mirrored(h):
    return ''.join(str(len(h) - int(x)) for x in h)


def neighborhoods(n):
    return [set(vv for vv in range(n*n) if abs(vv%n-v%n)+abs(vv//n-v//n)==1) for v in range(n*n)]


def local_data(h):
    n = len(h)
    H = [v for v in range(n*n) if v % n < int(h[v // n])]
    adj = neighborhoods(n)
    inside = set(H)
    neighbors = {v: adj[v] & inside for v in H}
    safe = {v for v in H if adj[v] <= inside}
    assert safe == {v for v in H if v % n < int(h[v // n]) - 1}
    return H, neighbors, safe


def brute_side(h):
    H, neighbors, safe = local_data(h)
    f, unsafe, bad_second = ([0] * (len(H)+1) for _ in range(3))
    for word in permutations(H):
        rank = dict(zip(word, range(len(H))))
        p = sum(all(rank[v] < rank[u] for u in neighbors[v]) for v in safe)
        f[p] += 1
        unsafe[p] += word[0] not in safe
        bad_second[p] += word[0] in safe and word[1] not in safe and word[1] not in neighbors[word[0]]
    return f, unsafe, bad_second


def brute_grid_n3(h):
    n = 3
    assert len(h) == n
    V = range(n*n)
    adj = neighborhoods(n)
    _, _, left = local_data(h)
    _, _, reflected_right = local_data(mirrored(h))
    right = {(v//n)*n+n-1-v%n for v in reflected_right}
    safe = left | right
    total, initial, refined, refined_two = 0, 0, 0, 0
    for word in permutations(V):
        rank = [0] * (n*n)
        for i, v in enumerate(word):
            rank[v] = i
        peaks = [v for v in V if all(rank[v] < rank[u] for u in adj[v])]
        count = sum(v in safe for v in peaks)
        total += len(peaks) == 2
        initial += count <= 2
        refined += count < 2 or (count == 2 and word[0] in safe)
        refined_two += count < 2 or (count == 2 and word[0] in safe
                                    and (word[1] in safe or word[1] in adj[word[0]]))
    return total, initial, refined, refined_two


def timed(job):
    h, ix, prime, out = job
    argv = ['/usr/bin/time', '-l', str(HERE / 'protected_mod'), h, str(ix), '--max-table-mib', '64']
    result = subprocess.run(argv, capture_output=True, text=True, check=True, timeout=20)
    record = json.loads(result.stdout)
    expected_sizes = independent_sizes(''.join(str(int(x)-1) for x in h))
    assert record['profile'] == h and record['prime_index'] == ix and int(record['prime']) == prime
    assert record['area'] == sum(map(int, h)) and record['protected_vertices'] == sum(map(int, h))-len(h)
    assert record['states_by_layer'] == record['processed_by_layer'] == expected_sizes
    assert record['independent_sets'] == sum(expected_sizes)
    size = re.search(r'^\s*(\d+)\s+maximum resident set size$', result.stderr, re.M)
    assert size
    raw = out / f'{h}-p{ix}.json'
    stat = out / f'{h}-p{ix}.resource'
    raw.write_text(result.stdout)
    stat.write_text(result.stderr)
    return (h, ix), record, {'profile': h, 'prime_index': ix, 'max_rss_bytes': int(size[1]),
                             'raw_sha256': digest(raw), 'resource_sha256': digest(stat)}


def region(h, records, primes):
    fields = ('peak_frequency_residues', 'unprotected_max_frequency_residues', 'bad_second_frequency_residues')
    freq, unsafe, bad_second = [[crt(column, primes) for column in zip(*(records[h, ix][field] for ix in range(3)))]
                    for field in fields]
    area = sum(map(int, h))
    assert sum(freq) == factorial(area)
    assert sum(unsafe) == len(h)*factorial(area-1)
    assert all(0 <= v <= u for u, v in zip(freq, unsafe))
    assert all(0 <= s <= u-v for u, v, s in zip(freq, unsafe, bad_second))
    if len(freq) < 3:
        return tuple(row + [0]*(3-len(freq)) for row in (freq, unsafe, bad_second))
    return freq, unsafe, bad_second


def merge(h, regional):
    a, b = h, mirrored(h)
    left, left_bad, left_second = regional[a]
    right, right_bad, right_second = regional[b]
    n = len(h)
    order_size = n*n
    size_left = sum(map(int, h))
    shuffle = comb(order_size, size_left)
    shuffle_left = comb(order_size-1, size_left-1)
    shuffle_right = comb(order_size-1, size_left)
    shuffle_left_left = comb(order_size-2, size_left-2)
    shuffle_right_right = comb(order_size-2, size_left)
    shuffle_split = comb(order_size-2, size_left-1)
    unfiltered = 0
    removed = 0
    removed_two = 0
    for i in range(3):
        for j in range(3-i):
            unfiltered += shuffle*left[i]*right[j]
            if i+j == 2:
                removed += shuffle_left*left_bad[i]*right[j] + shuffle_right*left[i]*right_bad[j]
                removed_two += shuffle_left_left*left_second[i]*right[j] + shuffle_right_right*left[i]*right_second[j]
                removed_two += shuffle_split*((left[i]-left_bad[i])*right_bad[j] + left_bad[i]*(right[j]-right_bad[j]))
    refined = unfiltered-removed
    refined_two = refined-removed_two
    assert 0 <= refined_two <= refined <= unfiltered <= factorial(order_size)
    return unfiltered, refined, refined_two


def decimal(r):
    with localcontext() as ctx:
        ctx.prec = 58
        return str(Decimal(r.numerator)/Decimal(r.denominator))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('out', type=Path)
    parser.add_argument('--workers', type=int, default=3)
    cli = parser.parse_args()
    cli.out.mkdir(exist_ok=False)
    audit = ('111', '112', '122', '222', '1122', '1112')
    cuts = ('112', '22233', '3333444', '44444444', '444445555', '444555544', '454545454', '232345678')
    sides = sorted(set(audit) | {u for h in cuts for u in (h, mirrored(h))})
    primes, nodes = verified_primes(HERE/'pratt-61bit.json')
    with ThreadPoolExecutor(max_workers=cli.workers) as pool:
        observed = list(pool.map(timed, [(h, i, p, cli.out) for h in sides for i, p in enumerate(primes)]))
    counts = {key: result for key, result, _ in observed}
    regional = {h: region(h, counts, primes) for h in sides}
    for h in audit:
        exact = brute_side(h)
        compact = regional[h]
        assert all(x == y + [0]*(len(x)-len(y)) for x, y in zip(exact, compact))
    exact_n3, upper_n3, refined_n3, refined_two_n3 = brute_grid_n3('112')
    assert exact_n3 == 143808 and (upper_n3, refined_n3, refined_two_n3) == merge('112', regional)
    totals = {int(n): int(value) for n, value in json.loads((HERE/'calibration-counts.json').read_text()).items()}
    upper = []
    for h in cuts:
        n = len(h)
        baseline, refined, refined_two = merge(h, regional)
        fact = factorial(n*n)
        record = {'n': n, 'profile': h, 'reflected_other_side': mirrored(h),
                  'both_safe_population': n*(n-2), 'protected_convolution_upper_A2': str(baseline),
                  'protected_max_refined_upper_A2': str(refined),
                  'unfiltered_upper_probability': decimal(Fraction(baseline, fact)),
                  'refined_upper_probability': decimal(Fraction(refined, fact)),
                  'protected_top_two_refined_upper_A2': str(refined_two),
                  'top_two_refined_upper_probability': decimal(Fraction(refined_two, fact)),
                  'unprotected_global_max_removed_words': str(baseline-refined),
                  'unprotected_nonadjacent_global_second_removed_words': str(refined-refined_two)}
        if n in totals:
            assert totals[n] <= refined_two
            record['known_A2'] = str(totals[n])
            record['overestimate_factor'] = decimal(Fraction(refined_two, totals[n]))
        upper.append(record)
        print(n, h, record['top_two_refined_upper_probability'], record.get('overestimate_factor', ''))
    ledger = {'scope': 'Rigorous finite upper counts for original K=2', 'primes': list(map(str, primes)),
              'pratt_nodes': nodes, 'python': sys.version, 'workers': cli.workers,
              'counter_sha256': digest(HERE/'protected_mod.cpp'), 'driver_sha256': digest(Path(__file__)),
              'jobs_waited': len(observed), 'max_worker_rss_bytes': max(m['max_rss_bytes'] for _,_,m in observed),
              'independent_small_side_gates': list(audit), 'n3_full_permutations_checked': factorial(9),
              'upper_bounds': upper, 'regions': {h: {'counts': list(map(str, f)), 'unsafe_max': list(map(str, g)),
                                                        'unsafe_second_after_safe_nonadjacent_max': list(map(str, j))}
                                               for h, (f,g,j) in regional.items()},
              'receipts': [m for _,_,m in observed]}
    (cli.out/'summary.json').write_text(json.dumps(ledger, indent=2)+'\n')


if __name__ == '__main__':
    main()
