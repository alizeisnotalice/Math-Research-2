#!/usr/bin/env python3
"""Audit layer metadata against independent binomial coverage counts.

This audits the scope and checked ledger hashes; it does not replace re-running
verify_support_layer.py on the full regenerated ledgers.
"""
from collections import Counter
import hashlib
import json
from math import comb
from pathlib import Path
import sys


def audit(path):
    report = json.loads(Path(path).read_text())
    base = Path(__file__).resolve().parent
    for name, digest in report['source_sha256'].items():
        assert hashlib.sha256((base/name).read_bytes()).hexdigest() == digest
    output = []
    for group in report['groups']:
        low, high, cap = group['min_total'],group['total'],group['max_terms']
        pairs = set()
        counts = Counter()
        for job in group['jobs']:
            s,v = job['summary'],job['verification']
            assert s['complete'] and v['verified']
            c = s['config']
            assert (c['min_total'],c['total'],c['max_terms']) == (low,high,cap)
            assert c['min_degree'] == c['max_first_degree']
            assert s['counts'] == v['counts']
            assert s['ledger_sha256_uncompressed'] == v['ledger_sha256_uncompressed']
            assert s['script_sha256'] == report['source_sha256']['support_search_layer.py']
            assert v['verifier_sha256'] == report['source_sha256']['verify_support_layer.py']
            subtotal = Counter()
            for shard in s['shards']:
                m,n = shard['m'],shard['n']
                assert shard['complete'] and m == c['min_degree']
                assert m <= n and low <= m+n <= high and (m,n) not in pairs
                pairs.add((m,n))
                def number(d):
                    return sum(comb(d-1,j) for j in range(min(d-1,cap-2)+1))
                assert shard['raw'] == number(m)*number(n)
                assert shard['canonical'] == sum(shard[k] for k in ('fair','excluded','unknown'))
                subtotal.update({k:shard[k] for k in s['counts']})
            assert dict(subtotal) == s['counts']
            assert len(job['unknown']) == s['counts']['unknown']
            counts.update(subtotal)
        expected = {(m,t-m) for t in range(low,high+1) for m in range(1,t//2+1)}
        assert pairs == expected
        assert dict(counts) == group['counts']
        output.append({'min_total':low,'total':high,'max_terms':cap,
                       'degree_pairs':len(pairs),'counts':dict(counts)})
    return {'coverage_metadata_verified':True,'groups':output,
            'note':'Exact record verification requires full ledger replay, separately recorded.'}


if __name__ == '__main__':
    print(json.dumps(audit(sys.argv[1]),indent=2))
