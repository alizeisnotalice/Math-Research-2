#!/usr/bin/env python3
"""Independent output-domain FWHT, rather than 32768 separate bent slices."""
import argparse,hashlib,json,sys,time
from pathlib import Path

def parity(x):return bin(x).count('1')%2
def fwht(a):
    h=1
    while h<len(a):
        for p in range(0,len(a),2*h):
            for i in range(p,p+h):a[i],a[i+h]=a[i]+a[i+h],a[i]-a[i+h]
        h*=2
    return a
def rank(rows):
    piv={}
    for y in rows:
        while y:
            b=y.bit_length()-1
            if b in piv:y^=piv[b]
            else:piv[b]=y;break
    return len(piv)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();beg=time.monotonic()
    raw=args.input.read_bytes();input_data=json.loads(raw);rows=[]
    C=[2**64-1]+[sum(((x>>j)&1)<<x for x in range(64)) for j in range(6)]
    for row in input_data['rows']:
        basis=C+row['relative_generators'];assert rank(basis)==15
        S=[sum(((b>>x)&1)<<j for j,b in enumerate(basis)) for x in range(64)]
        assert len(set(S))==64 and rank(S)==15
        phases=[parity((x&7)&row['permutation'][x>>3])^((row['h_bits']>>(x>>3))&1) for x in range(64)]
        spectrum=[0]*32768
        for a,b in zip(S,phases):spectrum[a]=1-2*b
        H=fwht(spectrum[:]);assert {abs(y) for y in H}=={8}
        # A one-coefficient sign mutation must be detected.
        bad=spectrum[:];bad[S[0]]*=-1;assert {abs(y) for y in fwht(bad)}!={8}
        SS=set(S);periods=[a^S[0] for a in S if {s^a^S[0] for s in S}==SS]
        rows.append({'source_file':row['source_file'],'source_sha256':row['source_sha256'],'row_index':row['row_index'],
                     'support':S,'coefficient_numerators_over_8':[1-2*b for b in phases],
                     'linear_rank':15,'affine_rank':rank([a^S[0] for a in S]),
                     'periods':sorted(periods),'period_dimension':rank(periods),
                     'all_32768_output_values_boolean':True,'one_coefficient_mutation_detected':True})
    assert len(rows)==20 and sum(r['period_dimension']==0 for r in rows)==10
    result={'status':'passed','algorithm':'one length-32768 integer FWHT per output function, plus a negative control',
            'input_sha256':hashlib.sha256(raw).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'python':sys.version,'rows':rows,'elapsed_seconds':time.monotonic()-beg}
    args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','P0_count':10,'elapsed_seconds':result['elapsed_seconds']}))
if __name__=='__main__':main()
