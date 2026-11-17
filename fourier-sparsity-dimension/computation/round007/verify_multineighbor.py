#!/usr/bin/env python3
"""Reconcile two complete quotient scans and certify all returned coset slices."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'round006'))
from verify_switch_closure import echelon, fwht
from analyze_minimum_switches import span


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--audit', type=Path, required=True)
    ap.add_argument('--results', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--candidates', type=Path, required=True)
    a = ap.parse_args()
    row = json.loads(a.audit.read_text())['rows'][1]
    assert Path(row['source_file']).name == 'corrected-results.json' and row['row_index'] == 1
    here = Path(__file__).resolve().parent
    header = (here / 'neighbor_input.hpp').read_text()
    rel_text = header.split('INPUT_REL={')[1].split('}')[0]
    assert [int(s, 16) for s in re.findall(r'0x[0-9a-f]+', rel_text)] == row['relative_generators']
    assert [int(x) for x in header.split('INPUT_PERM={')[1].split('}')[0].split(',')] == row['permutation']
    assert int(header.split('INPUT_H_BITS=')[1].split(';')[0]) == row['h_bits']
    C = [(1 << 64)-1] + [sum(((x >> j) & 1) << x for x in range(64)) for j in range(6)]
    rel = row['relative_generators'][4:]
    G = C + rel
    piv = echelon(G)
    assert len(piv) == 11

    def rem(x):
        for b in sorted(piv, reverse=True):
            if x >> b & 1:
                x ^= piv[b][0]
        return x

    receipt = a.results / 'remote-multineighbor-execution.json'
    remote = json.loads(receipt.read_text())
    r = remote['result']
    assert remote['ssh_exit_code'] == 0
    assert remote['runner_sha256'] == sha(here / 'run_remote_multineighbor.py')
    for n, h in r['source_hashes'].items():
        assert sha(here / n) == h
    shards = sorted(r['shards'], key=lambda x: x['summary']['start'])
    at = 0
    for shard in shards:
        s = shard['summary']
        assert shard['exit_code'] == 0 and s['status'] == 'passed'
        assert s['start'] == at and s['count'] == s['tested']
        assert s['W_rank'] == 11 and s['quotient_rank'] == 31
        at += s['count']
        hits = [h for h in r['found'] if s['start'] <= h['index'] < at]
        assert len(hits) == s['hits']
    assert at == r['total_tested'] == 1 << 31 and len(shards) == 14
    assert len({h['index'] for h in r['found']}) == len(r['found']) == 464
    rh = {rem(int(h['relative_phase_hex'], 16)) for h in r['found']}
    lh = set()
    summaries = []
    for i in range(2):
        summary = json.loads((a.results / f'local-independent-{i}.summary.json').read_text())
        summaries.append(summary)
        assert summary['start'] == i * (1 << 30) and summary['count'] == 1 << 30
        assert summary['W_rank'] == 11 and summary['quotient_rank'] == 31
        assert summary['criterion'] == 'all_63_derivatives_balanced'
        records = [json.loads(l) for l in (a.results / f'local-independent-{i}.found.jsonl').read_text().splitlines()]
        assert len(records) == summary['hits']
        assert len({h['index'] for h in records}) == len(records)
        for h in records:
            assert summary['start'] <= h['index'] < summary['start'] + summary['count']
            value = rem(int(h['relative_phase_hex'], 16))
            assert value == int(h['class_normal_form_hex'], 16)
            lh.add(value)
    assert rh == lh and len(rh) == 464
    q = [(((x & 7) & row['permutation'][x >> 3]).bit_count() & 1)
         ^ ((row['h_bits'] >> (x >> 3)) & 1) for x in range(64)]
    seed = sum(bit << x for x, bit in enumerate(q))
    for h in r['found']:
        assert int(h['phase_hex'], 16) == seed ^ int(h['relative_phase_hex'], 16)
    # Recheck canonical representatives; all affine C translates preserve bentness.
    for v in rh:
        for z in span(rel):
            assert all(abs(x) == 8 for x in fwht([1-2*(q[x] ^ (((z ^ v) >> x) & 1)) for x in range(64)]))
    periods = sorted(v for v in rh if {x ^ v for x in rh} == rh)
    candidates = {'W_generators': G, 'K': sorted(rh), 'period': periods,
                  'original_rels': row['relative_generators']}
    a.candidates.write_text(json.dumps(candidates, indent=2)+'\n')
    result = {'status': 'passed', 'arithmetic': 'integer and F2', 'tolerance': 0,
              'remote_tested': at, 'local_tested': sum(s['count'] for s in summaries),
              'quotient_classes': len(rh), 'returned_coset_slices_verified': len(rh)*16,
              'quotient_span_rank': len(echelon(rh)), 'period_rank': len(echelon(periods)),
              'scope': 'fixed dimension11 W, extensions inside RM(3,6); not general A28',
              'independent_criteria': ['integer Walsh', 'all 63 nonzero derivatives balanced'],
              'local_summaries': summaries, 'remote_receipt_sha256': sha(receipt),
              'audit_sha256': sha(a.audit), 'candidate_sha256': sha(a.candidates),
              'source_sha256': sha(Path(__file__)), 'python': sys.version,
              'platform': platform.platform()}
    a.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
