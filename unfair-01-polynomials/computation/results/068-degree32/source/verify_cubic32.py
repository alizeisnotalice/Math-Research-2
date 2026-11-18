"""Verify that both midpoint choices produce their signed cubic in row 24."""
from argparse import ArgumentParser
import json
from pathlib import Path
import sys

from sympy import ZZ, __version__ as sympy_version
from sympy.polys.rings import ring

sys.dont_write_bytecode=True
from census32 import descriptor, sha
from prefix_29_30 import model, slots


def verify(shape):
    specification=descriptor(shape)
    positions=slots(32,specification)
    R,t,*u=ring(['t',*[f'u{k}' for k in positions]],ZZ)
    rows=dict(model(t,32,specification,u)[0])
    assert all(sum(m[1:])<=2 for k,f in rows.items() if k<24 for m in f)
    cubic={m:c for m,c in rows[24].items() if sum(m[1:])>=3}
    d=1 if shape[-1]=='A' else -1
    assert cubic=={(0,3)+(0,)*(len(u)-1):d}
    return {'shape':shape,'variable_indices':positions,'row24_cubic_coefficient':d}


def main():
    if not __debug__:
        raise SystemExit('assertions required')
    p=ArgumentParser(description=__doc__)
    p.add_argument('--census',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    rows=json.loads(a.census.read_bytes())['ledger']['candidates']
    result=[]
    for mid in 'AB':
        name=next(x['shape'] for x in rows if x['shape'][0]=='X' and x['shape'][-1]==mid)
        result.append(verify(name))
    output={'status':'verified','rows':result,'census_sha256':sha(a.census),'source_sha256':sha(Path(__file__)),
            'python':sys.version,'sympy':sympy_version}
    with a.out.open('x') as f:
        f.write(json.dumps(output,sort_keys=True,indent=2)+'\n')
    print(result,flush=True)


if __name__=='__main__':
    main()
