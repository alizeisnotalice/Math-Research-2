#!/usr/bin/env python3
"""Audit the local four-word averaging inequality for monotone path h values."""
from __future__ import annotations

import hashlib
import itertools
import json
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
CERT_ROOT = Path(os.environ.get('A022_CERTIFICATES_ROOT', HERE.parent.parent))
INPUT = CERT_ROOT / 'r026' / 'r026-monotone-prefix-scan' / 'monotone-prefix-ie-scan.json'


def words(k: int):
    return map(''.join, itertools.product('RD', repeat=k))


def main() -> None:
    raw = INPUT.read_bytes()
    source = json.loads(raw)
    out = {
        'kind': 'exact-four-word-turn-convexity-audit',
        'input_sha256': hashlib.sha256(raw).hexdigest(),
        'definition': 'Delta=h(pRDq)+h(pDRq)-h(pRRq)-h(pDDq)',
        'scope': 'Hash-guarded integer post-processing of existing n=5,6 all-path completion data; no completion count is re-run.',
        'runs': [],
    }
    for run in source['runs']:
        values = {r['word']: int(r['hS']) for r in run['rows']}
        ell = run['n'] - 1
        rows = []
        for i in range(ell - 1):
            for p in words(i):
                for q in words(ell - i - 2):
                    delta = values[p + 'RD' + q] + values[p + 'DR' + q] - values[p + 'RR' + q] - values[p + 'DD' + q]
                    rows.append({'position': i, 'prefix': p, 'suffix': q, 'delta': str(delta)})
        plus = [z for z in rows if int(z['delta']) > 0]
        minus = [z for z in rows if int(z['delta']) < 0]
        assert plus and minus
        out['runs'].append({'n': run['n'], 'comparisons': len(rows),
                            'positive_count': len(plus), 'negative_count': len(minus),
                            'zero_count': len(rows) - len(plus) - len(minus),
                            'positive_example': plus[0], 'negative_example': minus[0]})
    text = json.dumps(out, indent=2) + '\n'
    HERE.joinpath('four-word-turn-convexity-audit.json').write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
