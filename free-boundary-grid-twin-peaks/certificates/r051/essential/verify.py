"""Verify essential contact tables against direct words and earlier marginals."""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
from math import comb
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
CERT = HERE.parents[1]
sys.path.insert(0, str(CERT / 'r048/merge-repair'))
from audit import graph, merge, peaks, shell


def integer_rows(file):
    return (tuple(map(int, line.split())) for line in file.read_text().splitlines())


def brute():
    neighbors = graph(3)
    result = Counter()
    for word in permutations(range(9)):
        position = {v: i for i, v in enumerate(word, 1)}
        roots = peaks(neighbors, position)
        if len(roots) != 2:
            continue
        p = word[0]
        q, = roots - {p}
        w, branch = merge(neighbors, word, q)
        rank = position[w]
        boundary = shell(neighbors, branch) - {w}
        dependent = [
            position[v] for v in boundary
            if all(position[v] < position[u] for u in neighbors[v] if u not in branch)
        ]
        reached = {q}
        length = 0
        while not any(w in neighbors[v] for v in reached):
            reached |= {u for v in reached for u in neighbors[v] if u in branch}
            length += 1
        distance = abs(p % 3 - q % 3) + abs(p // 3 - q // 3)
        result[rank, len(branch), distance, length, min(dependent, default=10) - rank] += 1
    return result


def audit(n):
    table = HERE / f'n{n}-wide.tsv'
    new = {}
    essential = defaultdict(Counter)
    contact = defaultdict(Counter)
    for i, b, d, length, h, mass in integer_rows(table):
        key = (i, b, d, length, h)
        assert key not in new
        assert b + 2 <= i < n*n and 1 <= b <= 8 and 2 <= d <= 2*(n-1)
        assert 0 <= length < b and 1 <= h <= n*n+1-i and mass > 0
        new[key] = mass
        essential[i,b,d][h] += mass
    for i, b, h, d, outside, mass in integer_rows(CERT/'r049/low-branch'/f'n{n}-contact.tsv'):
        contact[i,b,d][h] += mass
    previous = CERT/'r050/delete-branch'/('n3.tsv' if n == 3 else f'n{n}-b8.tsv')
    assert essential.keys() == contact.keys()
    assert essential.keys() == {row[:3] for row in integer_rows(previous)}
    comparisons = 0
    totals = defaultdict(Counter)
    for i, b, d, count, deleted in integer_rows(previous):
        ess = essential[i,b,d]
        raw = contact[i,b,d]
        sentinel = n*n+1-i
        assert sum(raw.values()) == sum(ess.values()) == count
        assert ess[sentinel] == deleted
        sr, se = count, count
        for cut in range(1, sentinel+1):
            sr -= raw[cut]
            se -= ess[cut]
            assert 0 <= sr <= se
            comparisons += 1
        for scope in ('all', str(b), f'd={d}'):
            acc = totals[scope]
            acc['words'] += count
            acc['deletable'] += deleted
            acc['raw_h1'] += raw[1]
            acc['essential_h1'] += ess[1]
            acc['raw_full_gap'] += sum(w for h,w in raw.items() if h >= i-b)
            acc['essential_full_gap'] += sum(w for h,w in ess.items() if h >= i-b)
    if n == 3:
        assert brute() == new
    if n == 5:
        for lower, upper, filename, resource, source in (
            (2,2,'mac-b2-finished.json','n5-mac-b2.resource','essential_b2_frozen.cpp'),
            (6,8,'mac-b6-8-finished.json','n5-mac-b6-8.resource','essential.cpp'),
        ):
            receipt = json.loads((HERE/filename).read_text())
            raw = ''.join(' '.join(map(str, row))+'\n' for row in integer_rows(table)
                          if lower <= row[1] <= upper)
            assert receipt['waited_reaped'] and receipt['returncode'] == 0
            assert receipt['table_sha256'] == sha256(raw.encode()).hexdigest()
            assert receipt['resource_sha256'] == sha256((HERE/resource).read_bytes()).hexdigest()
            assert receipt['source_sha256'] == sha256((HERE/source).read_bytes()).hexdigest()
    for (i,b,d,length,h), mass in new.items():
        for scope in ('all', str(b), f'd={d}'):
            totals[scope]['root_at_contact'] += mass * (length == 0)
            totals[scope]['root_at_or_next_to_contact'] += mass * (length <= 1)
    for acc in totals.values():
        acc['deletion_rate'] = str(Fraction(acc['deletable'], acc['words']))
    return {'n': n, 'nonzero_cells':len(new), 'groups':len(essential), 'tail_comparisons':comparisons,
            'table_sha256':sha256(table.read_bytes()).hexdigest(), 'counts':totals}


if __name__ == '__main__':
    print(json.dumps([audit(n) for n in (3,4,5)], indent=2))
