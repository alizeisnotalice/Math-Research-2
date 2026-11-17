#!/usr/bin/env python3
"""Certify all 20 switch partitions, membership, and joint normal rank exactly."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time


def echelon(vectors):
    pivots = {}
    for j, x in enumerate(vectors):
        label = 1 << j
        while x:
            p = x.bit_length() - 1
            if p in pivots:
                y, code = pivots[p]
                x ^= y
                label ^= code
            else:
                pivots[p] = (x, label)
                break
    return pivots


def coordinates(x, pivots):
    label = 0
    while x:
        p = x.bit_length() - 1
        if p not in pivots:
            raise AssertionError('not in span')
        y, code = pivots[p]
        x ^= y
        label ^= code
    return label


def fwht(values):
    a = list(values)
    step = 1
    while step < len(a):
        for start in range(0, len(a), 2 * step):
            for j in range(start, start + step):
                x, y = a[j], a[j + step]
                a[j], a[j + step] = x + y, x - y
        step *= 2
    return a


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('audit', type=Path)
    ap.add_argument('partitions', type=Path)
    ap.add_argument('--source-dir', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    begun = time.monotonic()
    rawa, rawp = args.audit.read_bytes(), args.partitions.read_bytes()
    audit, parts = json.loads(rawa), json.loads(rawp)
    assert parts['audit_sha256'] == sha(rawa)
    lookup = {(r['source_file'], r['row_index']): r for r in audit['rows']}
    assert len(lookup) == len(parts['rows']) == 20
    full = (1 << 64) - 1
    C = [full] + [sum(((x >> j) & 1) << x for x in range(64)) for j in range(6)]
    results = []
    for item in parts['rows']:
        row = lookup[(item['source_file'], item['row_index'])]
        source = (args.source_dir / Path(row['source_file']).name).read_bytes()
        assert sha(source) == row['source_sha256']
        original = json.loads(source)['rows'][row['row_index']]
        assert original['affine_flat_basis_masks'] == row['relative_generators']
        assert original['permutation'] == row['permutation'] and original['h_bits'] == row['h_bits']
        G = C + row['relative_generators']
        pivots = echelon(G)
        assert len(G) == len(pivots) == 15
        blocks = [int(x, 16) for x in item['partition_masks_hex']]
        assert len(blocks) == 8
        labels, union = [], 0
        for block in blocks:
            assert block.bit_count() == 8 and not (union & block)
            union |= block
            labels.append(coordinates(block, pivots))
        assert union == full and len(echelon(blocks)) == 8
        q = [(((x & 7) & row['permutation'][x >> 3]).bit_count() & 1)
             ^ ((row['h_bits'] >> (x >> 3)) & 1) for x in range(64)]
        all_normals, block_outputs, directions = [], [], []
        for block in blocks:
            points = [x for x in range(64) if block >> x & 1]
            b = points[0]
            dirs = [x ^ b for x in points]
            ubasis = [x for x, _ in echelon(dirs).values()]
            assert len(ubasis) == 3 and len(set(dirs)) == 8
            directions.append(sorted(dirs))
            for w in ubasis:
                slopes = sum((((g >> (b ^ w)) ^ (g >> b)) & 1) << j for j, g in enumerate(G))
                all_normals.append(slopes | (w << 15))
            for values in [q] + [[(g >> x) & 1 for x in range(64)] for g in G]:
                for w in ubasis:
                    slope = values[b ^ w] ^ values[b]
                    assert all((values[x ^ w] ^ values[x]) == slope for x in points)
            out = fwht([(1 - 2 * q[x]) if block >> x & 1 else 0 for x in range(64)])
            assert set(out) <= {-8, 0, 8} and sum(v != 0 for v in out) == 8
            block_outputs.append(out)
        assert all(sum(out[x] != 0 for out in block_outputs) == 1 for x in range(64))
        joint_rank = len(echelon(all_normals))
        slope_rank = len(echelon([x & ((1 << 15) - 1) for x in all_normals]))
        assert joint_rank == slope_rank == 7
        results.append({'source_name': Path(row['source_file']).name, 'row_index': row['row_index'],
                        'dimV': 15, 's': 64, 'blocks': 8, 'all_blocks_in_V': True,
                        'block_coordinates_hex': [hex(x) for x in labels],
                        'distinct_block_directions': len({tuple(x) for x in directions}),
                        'joint_normal_rank': joint_rank, 'slope_rank': slope_rank,
                        'dual_blocks_partition': True, 'closure_bound': 15})
    result = {'status': 'passed', 'arithmetic': 'integer and F2', 'tolerance': 0,
              'python': sys.version, 'platform': platform.platform(), 'elapsed_seconds': time.monotonic()-begun,
              'audit_sha256': sha(rawa), 'partitions_sha256': sha(rawp),
              'source_sha256': sha(Path(__file__).read_bytes()), 'rows': results}
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'certificates': len(results),
                      'all_partitions_membership_checked': True,
                      'nonparallel_partitions': sum(r['distinct_block_directions'] > 1 for r in results)}))


if __name__ == '__main__':
    main()
