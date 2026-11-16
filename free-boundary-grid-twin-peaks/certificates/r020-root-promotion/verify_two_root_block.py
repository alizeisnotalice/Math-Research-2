#!/usr/bin/env python3
"""Cross-check full n=3 enumeration, n=3 IE, and n=4 IE identities."""
import json
from pathlib import Path

base = Path(__file__).parent
full = json.loads((base / 'two-root-block-identity-n3.json').read_text())
ie3 = json.loads((base / 'two-root-block-ie-n3.json').read_text())
ie4 = json.loads((base / 'two-root-block-ie-n4.json').read_text())

assert int(full['N_falling_2_A1']) == int(ie3['N_falling_2_A1'])
assert int(full['RHS']) == int(ie3['RHS'])
for key in ('weighted_by_peak_count', 'K2_weighted_by_manhattan_distance'):
    for item in set(full[key]) | set(ie3[key]):
        assert int(full[key].get(item, 0)) == int(ie3[key].get(item, 0))
assert int(ie4['N_falling_2_A1']) == int(ie4['RHS'])
print(json.dumps({'status': 'passed',
                  'n3_full_vs_IE': True,
                  'n4_identity': True}, indent=2))
