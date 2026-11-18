"""Independent exact small rooted counts and finite root-stratified gate."""
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, localcontext
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys

from run_parity_n8 import HERE, COUNTER, complement, digest, geometry, independent_sets, nonintersection, side_region

REPO = HERE.parents[2]
DIAGONAL = REPO / 'certificates/r042/n5-lipschitz-ledger.json'
NORMALIZER = HERE.parents[1] / 'r044/inputs/parallelhash-generic-n5.json'


def binary(n, h):
    work = subprocess.run([str(COUNTER), str(n), h], capture_output=True,
                          text=True, check=True, timeout=40)
    result = json.loads(work.stdout)
    assert work.stderr == '' and result['profile'] == h and result['independent_sets'] == independent_sets(h)
    return list(map(int, result['root_counts']))


def subset_completion(n, h):
    """Count exact rooted growth orders by traversing all subsets backward."""
    H = side_region(h)
    vertices = sorted(H)
    adjacent = geometry(n)
    to_bit = {v: 1 << i for i, v in enumerate(vertices)}
    neighbors = [sum(to_bit[u] for u in adjacent[v] & H) for v in vertices]
    full = (1 << len(vertices)) - 1
    frontier = [0] * (full + 1)
    for s in range(1, full + 1):
        low = s & -s
        frontier[s] = frontier[s ^ low] | neighbors[low.bit_length() - 1]
    dp = [0] * (full + 1)
    dp[full] = 1
    for s in range(full - 1, 0, -1):
        allowed_next = frontier[s] & (full ^ s)
        while allowed_next:
            low = allowed_next & -allowed_next
            dp[s] += dp[s | low]
            allowed_next ^= low
    roots = [0] * (n * n)
    for v in vertices:
        if v % n < (n + 3) // 4 and adjacent[v] <= H:
            roots[v] = dp[to_bit[v]]
    return roots


def application_gate():
    n = 5
    ledger = json.loads(DIAGONAL.read_text())
    profiles = ledger['curves']
    known = {job['curves'][0]: int(job['count']) for job in ledger['jobs'] if job['curves'][0] == job['curves'][1]}
    with ThreadPoolExecutor(max_workers=3) as pool:
        rooted = dict(zip(profiles, pool.map(lambda h: binary(n, h), profiles)))
    groups = defaultdict(list)
    classes = defaultdict(lambda: defaultdict(int))
    for h in profiles:
        parity = ''.join(str(int(x) % 2) for x in h)
        groups[parity].append(h)
        hc = complement(h)
        total = math.comb(n * n, sum(map(int, h)))
        for p, a in enumerate(rooted[h]):
            if not a:
                continue
            for reflected_q, b in enumerate(rooted[hc]):
                if not b:
                    continue
                q = reflected_q // n * n + (n - 1 - reflected_q % n)
                pair = tuple(sorted((p, q)))
                turned = tuple(sorted(v % n * n + v // n for v in pair))
                classes[parity, 0][pair] += total * a * b
                classes[parity, 1][turned] += total * a * b
        representative = min(h, h[::-1], hc, hc[::-1])
        assert known[representative] == total * sum(rooted[h]) * sum(rooted[hc])
    adjacency = geometry(n)
    positions = {h: side_region(h) for h in profiles}
    pairs = [(h, k) for members in groups.values() for h, k in itertools.combinations(members, 2)]
    assert all(nonintersection(positions[h], positions[k], adjacency) for h, k in pairs)
    all_pairs = set().union(*classes.values())
    rootwise = sum(max(cls.get(pair, 0) for cls in classes.values()) for pair in all_pairs)
    best_class = max(sum(cls.values()) for cls in classes.values())
    assert rootwise > best_class
    a2 = int(json.loads(NORMALIZER.read_text())['A2'])
    frac = Fraction(rootwise, a2)
    with localcontext() as ctx:
        ctx.prec = 55
        decimal = str(Decimal(rootwise) / Decimal(a2))
    return {'n': n, 'profiles_each_orientation': len(profiles), 'classes_both_orientations': len(classes),
            'compatible_root_pairs': len(all_pairs), 'checked_no_bridge_pairs_vertical': len(pairs),
            'small_children_waited': len(rooted), 'best_single_class_count': str(best_class),
            'rootwise_disjoint_count': str(rootwise), 'normalizer_A2': str(a2),
            'rootwise_conditional_lower_exact': str(frac), 'rootwise_conditional_lower_decimal': decimal,
            'r042_sha256': digest(DIAGONAL), 'A2_sha256': digest(NORMALIZER)}


def main():
    gates = {4: ('1111', '2222', '1232', '2323', '3333'),
             5: ('11111', '12121', '22334', '22233', '34343', '32434', '43434')}
    checked = []
    for n, profiles in gates.items():
        for h in profiles:
            expected = subset_completion(n, h)
            result = binary(n, h)
            assert result == expected
            checked.append({'n': n, 'profile': h, 'roots': result})
    full = application_gate()
    output = {'python': sys.version, 'source_sha256': digest(HERE / 'cut_region_rootwise.cpp'),
              'independent_subset_root_recurrence': checked, 'n5_r042_gate': full}
    with Path(sys.argv[1]).open('x') as stream:
        stream.write(json.dumps(output, indent=2) + '\n')
    print(len(checked), full)


if __name__ == '__main__':
    main()
