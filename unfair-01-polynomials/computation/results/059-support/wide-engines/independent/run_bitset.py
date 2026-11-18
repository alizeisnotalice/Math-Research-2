"""Run and verify coverage of a compact, freshly recomputed support scan."""
import argparse
import concurrent.futures
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time


def support_count(degree, cap):
    return sum(math.comb(degree - 1, k) for k in range(min(degree - 1, cap - 2) + 1))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if not __debug__:
        raise RuntimeError('run without -O')
    parser = argparse.ArgumentParser()
    parser.add_argument('--total', type=int, nargs='+', required=True)
    parser.add_argument('--cap', type=int, required=True)
    parser.add_argument('--m', type=int, nargs='+')
    parser.add_argument('--workers', type=int, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    binary = source / 'support_bitset'
    args.out.mkdir(parents=True)
    stamp = time.monotonic()
    jobs = []
    for total in sorted(set(args.total)):
        for m in args.m or range(1, total // 2 + 1):
            n = total - m
            assert 1 <= m <= n
            outer = support_count(m, args.cap)
            inner = support_count(n, args.cap)
            chunk = max(1, 12_000_000 // inner)
            for start in range(0, outer, chunk):
                jobs.append((m, n, start, min(outer, start + chunk), outer, inner))
    assert len(set(jobs)) == len(jobs)
    # Schedule expensive jobs first, then verify full, disjoint outer coverage.
    jobs.sort(key=lambda row: (row[3] - row[2]) * row[5], reverse=True)
    fields = ('raw', 'canonical', 'excluded', 'fair', 'unknown')

    def run(job):
        m, n, start, end, outer, inner = job
        cmd = [str(binary), *map(str, (m, n, args.cap, start, end))]
        data = json.loads(subprocess.check_output(cmd, text=True))
        assert tuple(data[key] for key in ('m', 'n', 'begin', 'end', 'outer', 'inner')) == job
        assert data['cap'] == args.cap and data['raw'] == (end - start) * inner
        assert data['canonical'] == data['excluded'] + data['fair'] + data['unknown']
        assert len(data['unknown_masks']) == data['unknown']
        return data

    results = []
    with (args.out / 'chunks.jsonl').open('x') as ledger, concurrent.futures.ThreadPoolExecutor(args.workers) as pool:
        for future in concurrent.futures.as_completed([pool.submit(run, job) for job in jobs]):
            result = future.result()
            results.append(result)
            print(json.dumps(result, separators=(',', ':')), file=ledger, flush=True)
    pairs = {}
    for row in sorted(results, key=lambda x: (x['m'], x['n'], x['begin'])):
        key = f"{row['m']},{row['n']}"
        pair = pairs.setdefault(key, {**dict.fromkeys(fields, 0), 'end': 0, 'outer': row['outer'], 'unknown_masks': []})
        assert pair['end'] == row['begin']
        pair['end'] = row['end']
        for k in fields:
            pair[k] += row[k]
        pair['unknown_masks'].extend(row['unknown_masks'])
    assert len(results) == len(jobs) and all(v['end'] == v['outer'] for v in pairs.values())
    out = {'status': 'complete', 'totals': args.total, 'cap': args.cap,
           'workers': args.workers, 'seconds': time.monotonic() - stamp,
           'chunks': len(results), 'counts': {k: sum(x[k] for x in pairs.values()) for k in fields},
           'pairs': pairs, 'source_sha256': sha(source / 'support_bitset.cpp'), 'binary_sha256': sha(binary),
           'ledger_sha256': sha(args.out / 'chunks.jsonl'), 'python': sys.version,
           'scope': 'Enumeration of the listed finite layers. Nonempty unknowns require further exact adjudication.'}
    (args.out / 'summary.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: out[k] for k in ('status', 'totals', 'counts', 'chunks', 'seconds')}), flush=True)


if __name__ == '__main__':
    main()
