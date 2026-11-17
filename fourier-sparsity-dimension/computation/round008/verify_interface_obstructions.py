#!/usr/bin/env python3
"""Two explicit exact counterexamples to invalid closure/partial-rank steps."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys


def walsh(values):
    a=list(values);step=1
    while step<len(a):
        for b in range(0,len(a),2*step):
            for j in range(b,b+step):a[j],a[j+step]=a[j]+a[j+step],a[j]-a[j+step]
        step*=2
    return a


def span(gens):
    out=[0]
    for g in gens:out += [x^g for x in out]
    return out


def rank(vectors):
    p={}
    for x in vectors:
        while x:
            b=x.bit_length()-1
            if b in p:x^=p[b]
            else:p[b]=x;break
    return len(p)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    # x = u0 + 2u1 + 4v0 + 8v1.
    f=sum((((x&3)&(x>>2)).bit_count()&1)<<x for x in range(16))
    def transform(mask):return walsh([1-2*((f^mask)>>x&1) for x in range(16)])
    A=sum(1<<x for x in range(16) if x>>2==0)
    B=sum(1<<x for x in range(16) if not(x&1) and not(x&8))
    assert (A,B)==(15,85)
    assert all(set(map(abs,transform(v)))=={4} for v in span([A,B]))
    union=A|B;bad=transform(union)
    assert union==95 and union.bit_count()==6 and set(map(abs,bad))=={0,4,8}
    assert bad==[-8,0,4,4,0,0,4,-4,-8,0,-4,-4,0,0,-4,4]
    C=[65535]+[sum(((x>>j)&1)<<x for x in range(16)) for j in range(4)]
    G=C+[A];assert rank(G)==6
    assert all(set(map(abs,transform(v)))=={4} for v in span(G))
    slopes=[sum((((v>>d)^(v>>0))&1)<<j for j,d in enumerate([1,2])) for v in G]
    assert rank(slopes)==2
    result={'status':'passed','arithmetic':'integer and F2','tolerance':0,
            'union_closure_counterexample':{'seed_phase':f,'generators':[A,B],
                'all_four_slices_bent':True,'union_mask':union,'union_walsh':bad,'union_legal':False},
            'partial_slope_counterexample':{'generators':G,'dimension':6,'all_64_slices_bent':True,
                'selected_block':A,'selected_blocks':1,'slope_rank':2,'incorrect_full_partition_bound':0},
            'python':sys.version,'platform':platform.platform(),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','counterexamples':2}))


if __name__=='__main__':main()
