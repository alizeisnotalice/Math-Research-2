"""Select original preorder leaf positions in a rational outer cover of rho_5 <= pi/15."""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import arb, ctx
from verify_tomar_fan_table import PINNED_MANIFEST

CUT = Fraction(20943951024, 10**11)


def leaf_mask(codes, part, angle_high):
    """Values 1, 2, 3 denote inside, crossing, strictly outside. Internal nodes are zero."""
    assert codes.dtype == np.uint8 and codes.ndim == 1
    # Scale by the numerator of angle_high and CUT's denominator. The upper
    # comparison endpoint doubles on each split in the angle coordinate.
    stride = angle_high.numerator * CUT.denominator
    threshold = 2 * CUT.numerator * angle_high.denominator
    root_bit = (part >> 8) & 1  # bit for coordinate 2, C-order over eleven axes
    stack = [(root_bit * stride, (root_bit + 1) * stride, threshold)]
    mask = np.zeros(len(codes), dtype=np.uint8)
    for pos, code in enumerate(codes.tolist()):
        lo, hi, bound = stack.pop()
        if code >= 253:
            mask[pos] = 1 if hi <= bound else 3 if lo > bound else 2
        elif code == 2:
            mid = lo + hi
            stack.append((mid, 2 * hi, 2 * bound))
            stack.append((2 * lo, mid, 2 * bound))
        else:
            assert 0 <= code < 11
            stack.extend(((lo, hi, bound), (lo, hi, bound)))
    assert not stack
    return mask


def run(source, output, masks):
    ctx.prec = 128
    assert arb.pi() / 15 < arb(CUT.numerator) / CUT.denominator
    manifest = (source / 'SHA256SUMS').read_bytes()
    assert hashlib.sha256(manifest).hexdigest() == PINNED_MANIFEST
    pinned = {n.removeprefix('./'): h for h, n in (line.split() for line in manifest.decode().splitlines())}
    if masks:
        masks.mkdir(parents=True, exist_ok=True)
    rows = []
    for part in range(2048):
        stem = f'part{part:04}'
        header = source / (stem + '.json')
        tree = source / (stem + '.npy')
        for path in (header, tree):
            assert hashlib.sha256(path.read_bytes()).hexdigest() == pinned[path.name]
        meta = json.loads(header.read_bytes())
        assert meta['part'] == part and meta['nparts'] == 2048 and meta['split0'] == [2]*11
        assert meta['ns'] == [3,5,7,9] and meta['lo'][2] == 0
        codes = np.load(tree, allow_pickle=False)
        mask = leaf_mask(codes, part, Fraction(meta['hi'][2]))
        internal, inside, crossing, outside = map(int, np.bincount(mask, minlength=4))
        assert internal + 1 == inside + crossing + outside
        rows.append([part, inside, crossing, outside, hashlib.sha256(mask.tobytes()).hexdigest()])
        if masks:
            np.save(masks / (stem + '.npy'), mask, allow_pickle=False)
    counts = [sum(row[k] for row in rows) for k in (1,2,3)]
    assert sum(counts) == 16193402
    result = dict(cut=str(CUT), manifest=PINNED_MANIFEST, counts=counts,
                  active_roots=sum(r[1]+r[2]>0 for r in rows),
                  columns=['part','inside','crossing','outside','mask_sha256'], rows=rows)
    output.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    print(counts, result['active_roots'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--masks', type=Path)
    args = parser.parse_args()
    run(args.source, args.output, args.masks)
