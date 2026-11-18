"""Check complete disjoint coverage and agreement of two exact enumerators."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

FIELDS = ('raw', 'canonical', 'excluded', 'fair', 'unknown')


def compare(directory, totals, cap):
    summary = json.loads((directory / 'summary.json').read_text())
    ledger = (directory / 'chunks.jsonl').read_bytes()
    assert summary['status'] == 'complete' and summary['totals'] == totals and summary['cap'] == cap
    assert hashlib.sha256(ledger).hexdigest() == summary['ledger_sha256']
    records = sorted((json.loads(line) for line in ledger.splitlines()), key=lambda row: (row['m'], row['n'], row['begin']))
    assert len(records) == summary['chunks']
    expected = {(m, t-m) for t in totals for m in range(1, t//2+1)}
    layers = {}
    seen = set()
    count = lambda d: sum(math.comb(d-1, k) for k in range(min(d-1, cap-2)+1))
    for pair, rows in itertools.groupby(records, lambda row: (row['m'], row['n'])):
        assert pair in expected and pair not in seen
        seen.add(pair)
        m, n = pair
        outer, inner = count(m), count(n)
        end = 0
        current = {key: 0 for key in FIELDS}
        masks = []
        for row in rows:
            assert row['cap'] == cap and row['outer'] == outer and row['inner'] == inner
            assert row['begin'] == end and end < row['end'] <= outer
            assert row['raw'] == (row['end']-end)*inner
            end = row['end']
            assert row['canonical'] == row['excluded'] + row['fair'] + row['unknown']
            assert len(row['unknown_masks']) == row['unknown']
            masks.extend(tuple(x) for x in row['unknown_masks'])
            for k in FIELDS:
                current[k] += row[k]
        assert end == outer and len(set(masks)) == len(masks)
        saved = summary['pairs'][f'{m},{n}']
        assert all(saved[k] == current[k] for k in FIELDS)
        assert saved['unknown_masks'] == [list(x) for x in masks]
        layer = layers.setdefault(m+n, {key: 0 for key in FIELDS})
        for k in FIELDS:
            layer[k] += current[k]
    assert seen == expected and {tuple(map(int, key.split(','))) for key in summary['pairs']} == expected
    assert {k: sum(row[k] for row in layers.values()) for k in FIELDS} == summary['counts']
    return summary, records, layers


def main():
    if not __debug__:
        raise RuntimeError('run without -O')
    p = argparse.ArgumentParser()
    p.add_argument('first', type=Path)
    p.add_argument('second', type=Path)
    p.add_argument('--totals', nargs='+', type=int, required=True)
    p.add_argument('--cap', type=int, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    x, rows_x, totals = compare(a.first, a.totals, a.cap)
    y, rows_y, totals_y = compare(a.second, a.totals, a.cap)
    assert rows_x == rows_y and totals == totals_y
    assert x['source_sha256'] != y['source_sha256']
    report = {'status':'verified_identical_full_chunk_records', 'cap':a.cap, 'layers':totals,
              'chunks':len(rows_x), 'first_source':x['source_sha256'], 'second_source':y['source_sha256'],
              'first_ledger':x['ledger_sha256'], 'second_ledger':y['ledger_sha256'],
              'first_seconds':x['seconds'], 'second_seconds':y['seconds']}
    with a.out.open('x') as f:
        f.write(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
