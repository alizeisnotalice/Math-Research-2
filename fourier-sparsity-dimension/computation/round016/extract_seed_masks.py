#!/usr/bin/env python3
"""Make a minimal, auditable numeric seed file from R16's JSON input."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--receipt', type=Path, required=True)
    args = ap.parse_args()
    data = json.loads(args.input.read_text())
    masks = data['relative_bent_phase_masks']
    if data.get('seed_count') != 896 or len(masks) != 896:
        raise ValueError('expected exactly 896 masks')
    if any(type(x) is not int or not 0 <= x < (1 << 16) for x in masks):
        raise ValueError('mask outside four-variable truth-table range')
    if len(set(masks)) != 896:
        raise ValueError('duplicate seed mask')
    args.out.write_text(''.join(f'{x}\n' for x in masks))
    receipt = {'schema': 'a28-r016-seed-extraction-v1', 'count': len(masks),
               'json_sha256': sha(args.input), 'text_sha256': sha(args.out),
               'extractor_sha256': sha(__file__), 'minimum': min(masks),
               'maximum': max(masks), 'contains_zero': 0 in masks}
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
