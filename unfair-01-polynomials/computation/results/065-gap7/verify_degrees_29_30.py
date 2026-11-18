"""Specialized late-word unit-ideal check over two finite residue fields."""
from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json,sys,time
sys.dont_write_bytecode=True

import sympy as sp
from sympy.polys import distributedmodules as dm
from sympy.polys.orderings import grevlex
from modular_core_29_30 import equations, forward_eliminate
from prefix_29_30 import late_codes,census,canon,slots,exact_system
from verify_degrees_25_26 import MODS, pack, unpack, check_modulus


def modvalue(word,mod):
    p,_,w=mod
    return sum(a*b for a,b in zip(word,w,strict=True))%p


@lru_cache(maxsize=24)
def at_shape(j,shape,mod):
    return equations(j,shape,mod)


def individual(j,shape,mod,z,witness=None,track=False):
    R,vs,initial=at_shape(j,shape,mod)
    E=initial[:-1]+[initial[-1]-vs[-1]+R(z)]
    result=None
    if witness is not None:
        result=[unpack(f,R) for f in witness]
    else:
        simplified=forward_eliminate(R,vs,E,track)
        if not simplified:return False,None
        unit=next(((f,h) for f,h in simplified if f.is_ground),None)
        if unit is not None:
            if not track:return True,None
            f,h=unit
            result=[a/f.const() for a in h]
        else:
            symbols=tuple(v.as_expr() for v in vs[:-1])
            expressions=[f.as_expr() for f,h in simplified]
            if track:
                basis,change=dm.sdm_groebner([dm.sdm_from_vector([f],grevlex,R.domain,gens=symbols)
                                              for f in expressions],dm.sdm_nf_mora,grevlex,R.domain,extended=True)
                for pol,mul in zip(basis,change,strict=True):
                    if any(any(e[1:]) for e,c in pol):continue
                    c=pol[0][1]
                    reorders=[R.from_expr(x)/c for x in dm.sdm_to_vector(mul,symbols,R.domain,n=len(simplified))]
                    result=[sum(s*each[k] for s,(_,each) in zip(reorders,simplified,strict=True)) for k in range(len(E))]
                    break
                if result is None:return False,None
            else:
                G=sp.groebner(expressions,*symbols,order='grevlex',modulus=mod[0])
                return len(G)==1 and G[0]==1,None
    assert result and len(result)==len(E) and sum(f*g for f,g in zip(E,result,strict=True))==R.one
    if track:return True,[pack(f,mod[0]) for f in result]
    return True,None


def exact_remaining(j,shape,bits):
    data=[]
    for k,coeff in exact_system(j,shape):
        data.append([k,[[a,b.coefficients] for a,b in sorted(coeff.items())]])
    return {'early_shape':shape,'late_words':bits,'var_indices':slots(j,shape),'early_equations_Qt':data}


def process(j,shapes,track,past=None,start_index=1):
    previous={} if past is None else {canon(r['shape']):r for r in past}
    counts=Counter()
    all_rows=[]
    remaining=[]
    for i,entry in enumerate(shapes,start_index):
        shape=canon(entry['shape'])
        byword=late_codes(j,shape,entry['A'],entry['B'],counts)
        counts['shapes']+=1
        counts['words']+=len(byword)
        result={}
        live=[]
        for word in byword:
            rejected=False
            for mod in MODS:
                p=mod[0]
                z=modvalue(word,mod)
                key=f'{p}:{z}'
                if key in result:
                    ok,ack=result[key]
                else:
                    supplied=previous[shape]['witnesses'].get(key) if past is not None else None
                    if past is not None and supplied is None:
                        ok,ack=False,None
                    else:
                        ok,ack=individual(j,shape,mod,z,supplied,track=track)
                    result[key]=(ok,ack)
                    counts[f'mod{p}_systems_checked']+=1
                if ok:
                    counts[f'mod{p}_words_rejected']+=1
                    rejected=True
                    break
            if not rejected:
                live.append(word)
        if live:
            counts['shapes_live']+=1
            counts['words_live']+=len(live)
            remaining.append(exact_remaining(j,shape,live))
            print(j,i,'live',shape,live,flush=True)
        if track or live:
            all_rows.append({'shape':shape,'survivors':live,
                             'witnesses':{key:ack for key,(ok,ack) in result.items() if ok and ack is not None}})
        if i%100==0:
            print(j,i,'progress, live',counts.get('shapes_live',0),flush=True)
    assert counts['words']+counts['tail_interval_rejected']==counts['support_word_count']
    assert sum(counts[f'mod{p}_words_rejected'] for p,*_ in MODS)+counts['words_live']==counts['words']
    return dict(counts),all_rows,remaining


def main():
    if not __debug__:raise SystemExit('run with assertions enabled')
    ap=ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--witness',type=Path)
    ap.add_argument('--scout',type=Path)
    ap.add_argument('--only',type=int,choices=(29,30))
    ap.add_argument('--offset',type=int,default=0)
    args=ap.parse_args()
    assert args.witness is None or args.scout is None
    if args.out.exists():raise SystemExit('refuse existing output')
    stamp=time.monotonic()
    past=None if args.witness is None else {r['degree']:r['modular_certificate'] for r in json.loads(args.witness.read_text())['degrees']}
    cache=None if args.scout is None else {r['degree']:r for r in json.loads(args.scout.read_text())}
    degrees=[]
    for j in (29,30) if args.only is None else (args.only,):
        for p,t,w in MODS: check_modulus(p,t,w)
        base=census(j) if cache is None else cache[j]
        print(j,'early',base['counts'],flush=True)
        c,r,u=process(j,base['late'][args.offset:],cache is None,None if past is None else past[j],args.offset+1)
        data={'degree':j,'early_counts':base['counts'],'modular_counts':c,
              'modular_certificate':r,'unresolved':u}
        degrees.append(data)
        print(j,c,flush=True)
    here=Path(__file__).parent
    src=('verify_degrees_17_18.py','verify_degree_19.py','verify_degrees_21_22.py','verify_degrees_23_24.py','verify_degrees_25_26.py',
         'prefix_29_30.py','modular_core_29_30.py',Path(__file__).name)
    report={'status':'complete','python':sys.version,'sympy':sp.__version__,'offset':args.offset,
            'mode':'scout' if cache is not None else 'verify' if past is not None else 'generate',
            'degrees':degrees,'source_sha256':{n:sha256((here/n).read_bytes()).hexdigest() for n in src},
            'seconds':time.monotonic()-stamp}
    if args.witness is not None:report['witness_hash']=sha256(args.witness.read_bytes()).hexdigest()
    with args.out.open('x') as f:f.write(json.dumps(report,sort_keys=True,separators=(',',':'))+'\n')


if __name__=='__main__':main()
