"""Turn two-sided endpoint certificates and a verified inner disk into area bounds."""
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

root = Path(__file__).resolve().parent
src = {name: (root / name).read_bytes() for name in ('secants.json','area-result.json','point60/finite-result.json')}
records = json.loads(src['secants.json'])
branches = records['results']
geometric = json.loads(src['area-result.json'])
assert F(geometric['radius_host_upper']) == F(59,100)
finite = json.loads(src['point60/finite-result.json'])
assert set(finite['rows']) == {f'{h}-{i}' for h in ('r155','r175') for i in ('r151','r216')}
sampling_loss = F('.001371001')
rows = {}
for host in ('r155','r175'):
    for body in ('r151','r216','sparse','five-points'):
        d = min(F(branches[f'{host}-{body}-{sign}']['lower']) for sign in (1,-1))
        escapes = {body:d}
        if body in ('r151','r216'):
            escapes[body+'-35'] = d - sampling_loss
        for name, depth in escapes.items():
            assert depth > 0 and (F('.59')+depth)**2 - F(1,4) < F(8,25)**2
            rows[f'{host}-{name}'] = {'support_lower':str(depth),'area_gain_lower':str(F(25,16)*depth**2)}
result = {'diameter_squared':records['diameter_squared'],'bounds':rows,
          'source_sha256':{n:sha256(b).hexdigest() for n,b in src.items()}}
(root/'strict-gains.json').write_text(json.dumps(result,indent=2)+'\n')
