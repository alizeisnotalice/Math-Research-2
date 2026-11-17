#!/usr/bin/env python3
"""Check each recorded 8-block partition exactly; no search is performed."""
import argparse
import hashlib
import json
from pathlib import Path

FULL = (1 << 64) - 1

def c_basis():
    return [FULL] + [sum(((x >> j) & 1) << x for x in range(64)) for j in range(6)]

parser = argparse.ArgumentParser()
parser.add_argument('partitions', type=Path)
parser.add_argument('audit', type=Path)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
rawp, rawa = args.partitions.read_bytes(), args.audit.read_bytes()
parts, audit = json.loads(rawp), json.loads(rawa)
lookup = {(r['source_file'], r['row_index']): r for r in audit['rows']}
checked = []
for row in parts['rows']:
    blocks = [int(x, 16) for x in row['partition_masks_hex']]
    if len(blocks) != 8 or any(x.bit_count() != 8 for x in blocks): raise AssertionError('weight')
    union = 0
    for x in blocks:
        if union & x: raise AssertionError('not disjoint')
        union |= x
    if union != FULL: raise AssertionError('not cover')
    if row['partition_codewords_hex'] is not None:
        cert = lookup[(row['source_file'], row['row_index'])]
        generators = c_basis() + cert['relative_generators']
        for mask, code in zip(blocks, [int(x, 16) for x in row['partition_codewords_hex']]):
            rebuilt = 0
            for j, g in enumerate(generators):
                if code >> j & 1: rebuilt ^= g
            if rebuilt != mask: raise AssertionError('codeword')
    checked.append({'source_file': row['source_file'], 'row_index': row['row_index'],
                    'method': row['method'], 'partition_verified': True})
args.out.write_text(json.dumps({'status':'passed','partitions_sha256':hashlib.sha256(rawp).hexdigest(),
                                'audit_sha256':hashlib.sha256(rawa).hexdigest(),'rows':checked},indent=2)+'\n')
