#!/usr/bin/env python3
"""Combine six D4 puncture representatives into the exact full n=5 L tail."""
from __future__ import annotations
import argparse, json
from fractions import Fraction
from pathlib import Path

HERE=Path(__file__).resolve().parent
FILES={0:('tail-full-n5-v0.json',4),1:('tail-full-n5-v1.json',8),2:('tail-full-n5-v2.json',4),
       6:('tail-full-n5-v6.json',4),7:('tail-full-n5-v7.json',4),12:('tail-full-n5-center.json',1)}

def one(path:Path):
    d=json.loads(path.read_text()); assert d['n']==5 and d['u128_checked'] and len(d['holes'])==1
    return d['holes'][0]

def main(table: Path, input_dir: Path, output: Path):
    ref=json.loads(table.read_text())
    T=None; A=B=0; classes=[]
    for v,(fn,mult) in FILES.items():
        x=one(input_dir/fn); assert x['v']==v
        h=list(map(int,x['histogram_L'])); a=sum((i+1)*z for i,z in enumerate(h)); b=sum(h)
        r=ref['per_vertex'][v]
        assert a==r['A2_v'] and b==r['B_v'],(v,a,b,r)
        if T is None:T=[0]*len(h)
        assert len(T)==len(h)
        for i,z in enumerate(h):T[i]+=mult*z
        A+=mult*a;B+=mult*b
        tails=[sum(h[i:]) for i in range(len(h))]
        ratios=[Fraction(tails[i+1],tails[i]) for i in range(len(tails)-1) if tails[i]]
        classes.append({'representative_v':v,'xy':x['xy'],'orbit_size':mult,'A2_v':str(a),'B_v':str(b),
                        'B_over_A':str(Fraction(b,a)),'max_local_tail_ratio':str(max(ratios,default=Fraction(0))),
                        'local_tail_contraction_4_5':all(q<=Fraction(4,5) for q in ratios)})
    assert A==ref['A2_total'] and B==ref['B_total'],(A,B,ref['A2_total'],ref['B_total'])
    tails=[sum(T[i:]) for i in range(len(T))]
    ratios=[Fraction(tails[i+1],tails[i]) for i in range(len(tails)-1) if tails[i]]
    out={'n':5,'status':'passed','source':'six D4 representative exact C++ runs; totals cross-checked against r027 inclusion-exclusion table',
         'A2_total':str(A),'B_total':str(B),'B_over_A2':str(Fraction(B,A)),
         'L_histogram_global':list(map(str,T)),'tail_T_ge_t':list(map(str,tails)),
         'max_global_consecutive_tail_ratio':str(max(ratios,default=Fraction(0))),
         'global_tail_contraction_4_5':all(q<=Fraction(4,5) for q in ratios),'classes':classes}
    output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--table',type=Path,required=True); ap.add_argument('--input-dir',type=Path,default=HERE); ap.add_argument('--output',type=Path,default=HERE/'n5-tail-summary.json')
    a=ap.parse_args(); main(a.table,a.input_dir,a.output)
