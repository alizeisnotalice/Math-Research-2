#!/usr/bin/env python3
"""Deterministic support enumeration; exact singleton-unit certificates, no floats."""
import argparse
import gzip
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import sys
import time


def supports(degree, max_terms):
    for size in range(2, min(degree + 1, max_terms) + 1):
        for interior in itertools.combinations(range(1, degree), size - 2):
            yield (0, *interior, degree)


def canonical(I, J):
    if (len(I), I) > (len(J), J) and I[-1] == J[-1]:
        return False
    if math.gcd(*I, *J) != 1:
        return False
    pair = (I, J)
    rev = (tuple(I[-1] - i for i in reversed(I)),
           tuple(J[-1] - j for j in reversed(J)))
    if I[-1] == J[-1] and (len(rev[0]), rev[0]) > (len(rev[1]), rev[1]):
        rev = (rev[1], rev[0])
    # Compare using the same total ordering used for swapping.
    key = lambda p: ((len(p[0]), p[0]), (len(p[1]), p[1]))
    return key(pair) <= key(rev)


def classify(I, J):
    diagonals = {}
    for i in I:
        for j in J:
            diagonals.setdefault(i + j, []).append((i, j))
    unitI, unitJ = {0, I[-1]}, {0, J[-1]}
    for edges in diagonals.values():
        if len(edges) == 1:
            i, j = edges[0]
            unitI.add(i)
            unitJ.add(j)
    collisions = False
    for k in sorted(diagonals):
        edges = diagonals[k]
        if len(edges) > 1:
            collisions = True
            for i, j in edges:
                if i in unitI and j in unitJ:
                    return "excluded", [k, i, j]
    return ("unknown" if collisions else "fair"), None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--total', type=int, default=24)
    ap.add_argument('--min-degree', type=int, default=6)
    ap.add_argument('--max-first-degree', type=int)
    ap.add_argument('--max-terms', type=int, default=6)
    ap.add_argument('--seconds', type=float, default=14400)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.min_degree < 1 or args.total < 2*args.min_degree or args.max_terms < 2 or args.seconds <= 0:
        ap.error('Require positive degrees/time, total>=2*min-degree and max-terms>=2')
    if args.max_first_degree is not None and args.max_first_degree < args.min_degree:
        ap.error('max-first-degree must be at least min-degree')
    if not Path('/Volumes/data').is_mount():
        raise SystemExit('data is not mounted')
    args.out.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    config = {k: v for k, v in vars(args).items() if k != 'out'}
    config['version'] = 2
    degree_pairs = [(m,t-m) for t in range(2*args.min_degree,args.total+1)
                    for m in range(args.min_degree,min(t//2,args.max_first_degree or t//2)+1)]
    counts = dict(raw=0, canonical=0, excluded=0, fair=0, unknown=0)
    digest = hashlib.sha256()
    shards = []
    complete = True
    with gzip.open(args.out/'ledger.jsonl.gz', 'wt') as ledger, open(args.out/'unknown.jsonl', 'w') as unknown:
        for total in range(2 * args.min_degree, args.total + 1):
            for m in range(args.min_degree, min(total // 2,args.max_first_degree or total//2) + 1):
                n = total - m
                shard_counts = dict(raw=0, canonical=0, excluded=0, fair=0, unknown=0)
                for I in supports(m, args.max_terms):
                    for J in supports(n, args.max_terms):
                        counts['raw'] += 1
                        shard_counts['raw'] += 1
                        if not canonical(I, J):
                            continue
                        status, witness = classify(I, J)
                        record = [I, J, status, witness]
                        line = json.dumps(record, separators=(',', ':')) + '\n'
                        ledger.write(line)
                        digest.update(line.encode())
                        for key in ('canonical', status):
                            counts[key] += 1
                            shard_counts[key] += 1
                        if status == 'unknown':
                            unknown.write(json.dumps({'I': I, 'J': J}) + '\n')
                    if time.monotonic() - start > args.seconds:
                        complete = False
                        break
                shards.append({'m': m, 'n': n, 'complete': complete, **shard_counts})
                state = {'config': config, 'complete': complete, 'counts': counts,
                         'shards': shards, 'ledger_sha256_uncompressed': digest.hexdigest(),
                         'elapsed_seconds': time.monotonic() - start,
                         'python': sys.version, 'pid': os.getpid(),
                         'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
                (args.out/'checkpoint.json').write_text(json.dumps(state, indent=2)+'\n')
                print(json.dumps(shards[-1]), flush=True)
                if not complete:
                    break
            if not complete:
                break
    state['complete'] = complete and len(shards) == len(degree_pairs)
    (args.out/'summary.json').write_text(json.dumps(state, indent=2)+'\n')
    print(json.dumps({'complete': state['complete'], 'counts': counts, 'seconds': time.monotonic()-start}), flush=True)


if __name__ == '__main__':
    main()
