#!/usr/bin/env python3
"""Independent exhaustive checker of support ledger and exclusion witnesses.

Does not import the searcher. Enumeration uses integer masks; certificates use
direct convolution counts. The uncompressed ledger is reproducible, not Git data.
"""
import argparse
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import time


def support_list(d, cap):
    result = []
    for bits in range(1 << (d - 1)):
        if bits.bit_count() <= cap - 2:
            result.append((0, *(k for k in range(1, d) if bits & (1 << (k-1))), d))
    return sorted(result, key=lambda p: (len(p), p))


def is_representative(I, J):
    g = 0
    for i in I + J:
        g = math.gcd(g, i)
    if g != 1:
        return False
    def key(p):
        return len(p[0]), p[0], len(p[1]), p[1]
    options = [(I,J), (tuple(sorted(I[-1]-i for i in I)), tuple(sorted(J[-1]-j for j in J)))]
    if I[-1] == J[-1]:
        options += [(b,a) for a,b in options[:]]
    return key((I,J)) == min(map(key, options))


def verify_record(record):
    il, jl, status, witness = record
    I, J = set(il), set(jl)
    counts = [sum(k-i in J for i in I) for k in range(max(I)+max(J)+1)]
    ui, uj = {0,max(I)}, {0,max(J)}
    for i in I:
        for j in J:
            if counts[i+j] == 1:
                ui.add(i)
                uj.add(j)
    if status == 'excluded':
        k,i,j = witness
        assert i in ui and j in uj and k == i+j and counts[k] >= 2, record
    elif status == 'fair':
        assert witness is None and all(v<=1 for v in counts), record
    elif status == 'unknown':
        assert witness is None and any(v>1 for v in counts), record
        assert all(counts[i+j] == 1 for i in ui for j in uj), record
    else:
        raise AssertionError(status)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('directory', type=Path)
    args = ap.parse_args()
    report = json.loads((args.directory/'summary.json').read_text())
    assert report['complete'], 'Partial run cannot certify full coverage'
    c = report['config']
    by_degree = {d:support_list(d,c['max_terms']) for d in range(c['min_degree'],c['total']-c['min_degree']+1)}
    counts = dict(raw=0, canonical=0, fair=0, excluded=0, unknown=0)
    digest = hashlib.sha256()
    start = time.monotonic()
    with gzip.open(args.directory/'ledger.jsonl.gz','rt') as ledger:
        for total in range(2*c['min_degree'],c['total']+1):
            for m in range(c['min_degree'],total//2+1):
                for I,J in itertools.product(by_degree[m],by_degree[total-m]):
                    counts['raw'] += 1
                    if not is_representative(I,J):
                        continue
                    line = next(ledger)
                    digest.update(line.encode())
                    record = json.loads(line)
                    assert tuple(record[0]) == I and tuple(record[1]) == J, (I,J,record)
                    verify_record(record)
                    counts['canonical'] += 1
                    counts[record[2]] += 1
        assert next(ledger,None) is None, 'Unexpected extra ledger record'
    assert counts == report['counts'], (counts,report['counts'])
    assert digest.hexdigest() == report['ledger_sha256_uncompressed']
    result = {'verified':True,'counts':counts,'ledger_sha256_uncompressed':digest.hexdigest(),
              'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'elapsed_seconds':time.monotonic()-start}
    (args.directory/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
