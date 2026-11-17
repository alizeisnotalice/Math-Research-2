#!/usr/bin/env python3
"""Parent cross-check: coordinate-derived D4 and direct cofactor CRT.

Prime and production certification belongs to n7-global-h45/source/.
This intentionally small check is a second arithmetic reconstruction only.
"""
import hashlib
import json
import math
from pathlib import Path

root=Path(__file__).resolve().parent/'n7-global-h45'
totals=[0]*5
multiplicities={}
for v in (0,7,8,14,15,16,21,22,23,24):
    rows=[json.loads((root/f'raw/v{v}/n7-v{v}-p{i}.json').read_text()) for i in range(4)]
    moduli=[int(d['prime']) for d in rows]
    modulus=math.prod(moduli)
    assert modulus>math.factorial(48)
    for i,d in enumerate(rows):
        assert (d['n'],d['v'],d['prime_index'])==(7,v,i)
    values=[sum(int(d[f'h{k}_residue'])*(modulus//m)*pow(modulus//m,-1,m) for d,m in zip(rows,moduli))%modulus for k in range(1,6)]
    assert all(0<=x<=math.factorial(48) for x in values)
    r,c=divmod(v,7)
    orbit={(r,c),(c,6-r),(6-r,6-c),(6-c,r),(r,6-c),(6-r,c),(c,r),(6-c,6-r)}
    multiplicities[str(v)]=len(orbit)
    totals=[a+len(orbit)*b for a,b in zip(totals,values)]
assert sum(multiplicities.values())==49
expected=json.loads((root/'reports/n7-global-h45-crt-audit.json').read_text())['global']['h']
assert totals==[expected[f'h{i}'] for i in range(1,6)]
assert totals[2]-3*totals[3]+2*totals[4]==0
print(json.dumps({'status':'passed','multiplicities':multiplicities,'h':totals,'R45':0,'scope':'Finite second CRT and D4 arithmetic check only','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
