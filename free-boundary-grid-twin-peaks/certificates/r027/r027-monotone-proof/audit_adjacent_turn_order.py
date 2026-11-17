#!/usr/bin/env python3
"""Hash-guarded exact audit of the adjacent-turn monotonicity candidate.

Input is the pre-existing serial integer-IE output for all corner monotone
prefixes at n=5,6.  We do not re-run its completion counter here.  The audit
tests whether swapping a local RD to DR has one fixed direction for h(S).
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
CERT_ROOT = Path(os.environ.get('A022_CERTIFICATES_ROOT', HERE.parent.parent))
INPUT = CERT_ROOT / 'r026' / 'r026-monotone-prefix-scan' / 'monotone-prefix-ie-scan.json'


def main() -> None:
    raw = INPUT.read_bytes()
    data = json.loads(raw)
    report = {
        'kind': 'exact-adjacent-turn-order-audit',
        'input_certificate_path': str(INPUT.relative_to(CERT_ROOT)),
        'input_sha256': hashlib.sha256(raw).hexdigest(),
        'scope': 'Exact integer h values already obtained by the cited serial IE runs. This audit only compares them; it is not an all-n calculation.',
        'runs': [],
    }
    for run in data['runs']:
        val = {row['word']: int(row['hS']) for row in run['rows']}
        comparisons = []
        for word, old in sorted(val.items()):
            for i in range(len(word) - 1):
                if word[i:i + 2] != 'RD':
                    continue
                moved = word[:i] + 'DR' + word[i + 2:]
                new = val[moved]
                comparisons.append({'word_RD': word, 'word_DR': moved,
                                    'h_RD': str(old), 'h_DR': str(new),
                                    'sign_h_DR_minus_h_RD': (new > old) - (new < old)})
        plus = [x for x in comparisons if x['sign_h_DR_minus_h_RD'] > 0]
        minus = [x for x in comparisons if x['sign_h_DR_minus_h_RD'] < 0]
        assert plus and minus
        report['runs'].append({
            'n': run['n'], 'comparisons': len(comparisons),
            'DR_larger_count': len(plus), 'DR_smaller_count': len(minus),
            'explicit_DR_larger': plus[0], 'explicit_DR_smaller': minus[0],
            'no_ties': all(x['sign_h_DR_minus_h_RD'] for x in comparisons),
        })
    out = json.dumps(report, indent=2) + '\n'
    HERE.joinpath('adjacent-turn-order-audit.json').write_text(out)
    print(out, end='')


if __name__ == '__main__':
    main()
