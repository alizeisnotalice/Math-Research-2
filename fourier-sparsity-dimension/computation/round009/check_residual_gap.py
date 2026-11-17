#!/usr/bin/env python3
"""Exact finite checks for the degree/second-weight and Boolean grid steps."""
import argparse, hashlib, itertools, json, sys, time
from pathlib import Path

def walsh(a):
    a=list(a);h=1
    while h<len(a):
        for i in range(0,len(a),2*h):
            for j in range(i,i+h):
                x,y=a[j],a[j+h];a[j],a[j+h]=x+y,x-y
        h*=2
    return a

def weight_scan_rm24():
    monomials=[m for m in range(16) if bin(m).count('1')<=2]
    masks=[sum((int(x&m==m)<<x) for x in range(16)) for m in monomials]
    counts={};cur=0
    for z in range(1<<len(masks)):
        if z:cur^=masks[(z&-z).bit_length()-1]
        w=bin(cur).count('1');counts[w]=counts.get(w,0)+1
    assert sorted(w for w in counts if w)[:2]==[4,6]
    return {'monomial_masks':monomials,'codewords':sum(counts.values()),'weight_counts':counts,'first_two_nonzero_weights':[4,6]}

def bool3_grid():
    counts={};cap_pass=0
    for mask in range(256):
        w=walsh([1-2*((mask>>x)&1) for x in range(8)])
        assert sum(a*a for a in w)==64 and all(a%2==0 for a in w)
        maximum=max(abs(a) for a in w);counts[maximum]=counts.get(maximum,0)+1
        cap_pass+=int(maximum<=3)
    assert cap_pass==0
    return {'boolean_functions':256,'max_unnormalized_Walsh_histogram':counts,
            'functions_with_all_normalized_coefficients_le_3_over_8':cap_pass}

def degree_check():
    bent=0;deg={};faces=0
    for mask in range(65536):
        a=walsh([1-2*((mask>>x)&1) for x in range(16)])
        if any(abs(v)!=4 for v in a):continue
        bent+=1;coef=[(mask>>x)&1 for x in range(16)]
        for b in range(4):
            for x in range(16):
                if x>>b&1:coef[x]^=coef[x^(1<<b)]
        d=max((bin(x).count('1') for x in range(16) if coef[x]),default=0)
        assert d<=2;deg[d]=deg.get(d,0)+1
        for fixed in range(4):
            for rhs in range(2):
                wt=sum((mask>>x)&1 for x in range(16) if (x>>fixed&1)==rhs)
                assert wt%2==0;faces+=1
    assert bent==896
    return {'n':4,'truth_tables':65536,'bent_functions':bent,'degree_counts':deg,'coordinate_three_face_parity_checks':faces}

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();start=time.monotonic()
    result={'status':'passed','arithmetic':'integer/F2; zero tolerance','scope':'finite checks only; general residual theorem uses analytic proof and cited RM second-weight theorem',
            'RM_2_4':weight_scan_rm24(),'three_variable_grid':bool3_grid(),'bent_degree':degree_check(),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'elapsed_seconds':time.monotonic()-start}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
