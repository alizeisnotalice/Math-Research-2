#!/usr/bin/env python3
"""Match all unknown ledger entries to independently rechecked certificates."""
import json
from pathlib import Path
from diagonal_domination import verify_certificate
from check_layer14_exceptions import check

base=Path(__file__).resolve().parent/'results'
ledger=json.loads((base/'014-layer-ledgers.json').read_text())
cases=json.loads((base/'014-domination.json').read_text())['cases']
key=lambda c:(tuple(c['I']),tuple(c['J']))
unknown=[key(c) for g in ledger['groups'] for j in g['jobs'] for c in j['unknown']]
assert len(set(unknown))==len(unknown)
assert sorted(map(key,cases))==sorted(unknown)
special={key(c) for c in check()['cases']}
seen=set()
direct=0
for c in cases:
    if c['certificate'] is not None:
        assert c['status']=='excluded'
        verify_certificate(c['I'],c['J'],c['certificate'])
        direct+=1
    else:
        assert key(c) in special
        seen.add(key(c))
assert seen==special
print(json.dumps({'all_unknowns_exactly_covered':True,'unknowns':len(unknown),
                  'direct_certificates':direct,'special_certificates':len(seen)}))
