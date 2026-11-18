"""Exact fixed-parity disjoint sums from the n=5 ledger and fresh n=6,7 counts."""
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import sys

from certify_disjoint import HERE, OLD_TOTALS, PROJECT, complement, compute_regions, neighbors, no_bridge, read, sha256, weight


def all_curves(n, widths, p, h=''):
    if len(h) == n:
        yield h
        return
    for width in widths:
        if p is not None and width % 2 != p[len(h)]:
            continue
        if not h or abs(width - int(h[-1])) <= 1:
            yield from all_curves(n, widths, p, h + str(width))


def parity_key(h):
    return ''.join(str(int(x) % 2) for x in h)


def count_summary(n, rows, A2):
    adjacent = neighbors(n)
    h_family = [row['profile'] for row in rows]
    pairs = list(combinations(h_family, 2))
    assert len({parity_key(h) for h in h_family}) == 1
    assert all(no_bridge(n, h, k, adjacent) for h, k in pairs)
    total = sum(int(row['count']) for row in rows)
    with localcontext() as ctx:
        ctx.prec = 55
        ratio = str(Decimal(total) / Decimal(A2))
    return {'n': n, 'width_parity': parity_key(h_family[0]), 'profiles': len(rows),
            'nonempty_events': sum(int(row['count']) > 0 for row in rows),
            'verified_no_bridge_pairs': len(pairs), 'A2': str(A2), 'union_count': str(total),
            'union_probability_decimal': ratio, 'union_probability_exact': str(Fraction(total, A2)),
            'terms': rows}


def main():
    output = Path(sys.argv[1])
    ledger_path = PROJECT / 'certificates/r042/n5-lipschitz-ledger.json'
    old = read(ledger_path)
    old_diagonal = {job['curves'][0]: int(job['count']) for job in old['jobs']
                    if job['curves'][0] == job['curves'][1]}
    schemes = defaultdict(list)
    all_n5_curves = set(all_curves(5, range(1, 5), None))
    assert all_n5_curves == set(old['curves']) and len(all_n5_curves) == 178
    for h in sorted(all_n5_curves):
        hc = complement(h, 5)
        representative = min(h, h[::-1], hc, hc[::-1])
        schemes[parity_key(h)].append({'profile': h, 'count': str(old_diagonal[representative])})
    n5_sum = read(OLD_TOTALS / 'parallelhash-generic-n5.json')['A2']
    n5 = [count_summary(5, rows, int(n5_sum)) for rows in schemes.values()]
    n5.sort(key=lambda row: int(row['union_count']), reverse=True)

    records = []
    n7_source = PROJECT / 'computation/results/n7-exact-r005.json'
    n6_source = OLD_TOTALS / 'parallelhash-generic-n6.json'
    configurations = [(6, range(1, 6), n6_source), (7, range(2, 6), n7_source)]
    families = []
    for n, widths, source in configurations:
        for start in (0, 1):
            parity = tuple((start + y) % 2 for y in range(n))
            family = tuple(all_curves(n, widths, parity))
            families.append((n, widths, source, family))
    keys = {(n, x) for n, _, _, family in families for h in family for x in (h, complement(h, n))}
    regions = compute_regions(keys, 3)
    for n, widths, source, family in families:
        result = count_summary(n, [weight(n, h, regions) for h in family], int(read(source)['A2']))
        records.append(dict(width_range=list(widths), A2_source=str(source), A2_sha256=sha256(source), **result))
    summary = {'scope': 'Exact finite fixed-parity unions', 'python': sys.version,
               'counter_source_sha256': sha256(HERE / 'cut_region.cpp'),
               'n5_ledger_sha256': sha256(ledger_path), 'region_children_completed_and_waited': len(regions),
               'n5_classes': n5, 'n6_n7_classes': records}
    with output.open('x') as fp:
        fp.write(json.dumps(summary, indent=2) + '\n')
    for block in n5[:3] + records:
        print(block['n'], block['width_parity'], block['profiles'], block['union_count'], block['union_probability_decimal'])


if __name__ == '__main__':
    main()
