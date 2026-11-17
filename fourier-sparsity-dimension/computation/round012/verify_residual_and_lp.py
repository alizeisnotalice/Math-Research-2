#!/usr/bin/env python3
"""Exact t4 controls and an independently assembled full MacWilliams LP check."""
import argparse,hashlib,json,math,sys,time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
def par(x):return bin(x).count('1')%2
def rank(vs):
    piv={}
    for v in vs:
        while v:
            i=v.bit_length()-1
            if i in piv:v^=piv[i]
            else:piv[i]=v;break
    return len(piv)
def kraw(n,j,w):return sum((-1)**i*math.comb(w,i)*math.comb(n-w,j-i) for i in range(max(0,j-n+w),min(j,w)+1))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--certificate',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();beg=time.monotonic()
    upper=F(79,1024);lower=F(61,768)
    assert F(1,12)>upper and lower-upper==F(7,3072)
    f0=lambda x:((x&1)*((x>>3)&1))^(((x>>1)&1)*((x>>2)&1))
    v=lambda x:((x&1)*((x>>2)&1))^(((x>>1)&1)*((x>>2)&1))^(((x>>1)&1)*((x>>3)&1))
    R=[p for p in range(64) if ((p>>4)&1)==v(p&15)]
    # All 31 codimension-one subcodes of the known six-dimensional residual,
    # containing 1_R; four basis directions are chosen in its five-dimensional quotient.
    summaries=[]
    for normal in range(1,32):
        ker=[]
        for y in range(1,32):
            if par(normal&y)==0 and rank(ker+[y])>len(ker):ker.append(y)
        assert len(ker)==4
        atom={}
        for p in R:
            label=sum(par(u&(p&31))<<j for j,u in enumerate(ker))
            atom.setdefault(label,[]).append(p)
        items=sorted(atom.items());r=len(items)
        P=[[sum((-1)**(f0(p&15)^par(xi&p)) for p in pts) for _,pts in items] for xi in range(64)]
        assert all(all(z==0 for z in row) for row in P[32:])
        for xi in range(32):
            assert all(sum(z*((-1)**par(h*1&label)) for (label,_),z in zip(items,P[xi])) in [-8,8] for h in range(16))
            if r<16:
                assert all(z%2==0 for z in P[xi]) and max(map(abs,P[xi]))>=4
        for i,(_,pts) in enumerate(items):assert sum(row[i]**2 for row in P)==64*len(pts)
        weights=[]
        for h in range(16):
            wt=sum(len(pts)*par(h&label) for label,pts in items);weights.extend([wt,32-wt])
        assert 8 not in weights and all(w in [0,12,16,20,32] for w in weights)
        summaries.append({'normal':normal,'atom_count':r,'atom_sizes':dict(Counter(len(pts) for _,pts in items)),'code_dimension':5})
    assert Counter(x['atom_count'] for x in summaries)=={12:30,16:1}
    raw=a.certificate.read_bytes();cert=json.loads(raw);assert cert['status']=='passed'
    A=[F(0)]*65;A[0]=A[64]=F(1)
    for w,z in zip(cert['weights'],cert['rational_A_weights']):A[w]=A[64-w]=F(z)
    assert sum(A)==65536 and all(x>=0 for x in A)
    assert all(A[i]==0 for i in list(range(1,8))+list(range(57,64)))
    rows=[]
    for j in range(65):
        B=sum(A[w]*kraw(64,j,w) for w in range(65))/65536
        T=F(kraw(64,j,28)+kraw(64,j,36),2);plus=(B+T)/2;minus=(B-T)/2
        assert B>=0 and plus>=0 and minus>=0,(j,B,T)
        if j%4==0:assert plus>=A[j]
        if j%4==2:assert minus>=A[j]
        if j==2:assert B==0
        rows.append({'weight':j,'A':str(A[j]),'dual_B':str(B),'signed_positive':str(plus),'signed_negative':str(minus)})
    assert any(x.denominator!=1 for x in A)
    fourth=[]
    for q in [8,16,32,64]:
        s=q*q;val=kraw(s,4,(s-q)//2);assert val==-s*(s-1)//12
        need=1+(-val//(s//4) if q==8 else -val);defect=(need-1).bit_length()
        fourth.append({'q':q,'K4':val,'quotient_deficiency_at_least':defect,'code_dimension_upper_bound':(s-defect)//2})
    out={'status':'passed','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,
         't4_scalar_bounds':{'upper':str(upper),'r13_big_atom_lower':str(lower),'strict_gap':str(lower-upper)},
         'real_residual_subcodes':summaries,'LP':{'q':8,'D':16,'certificate_sha256':hashlib.sha256(raw).hexdigest(),'all_65_weight_constraints_checked':True,'integer_counts':False,'rows':rows},
         'fourth_defect':fourth,'elapsed_seconds':time.monotonic()-beg,
         'scope':'Finite controls; rational LP is only a continuous relaxation, not a code or a counterexample.'}
    a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'passed','real_subcodes':31,'LP_integer_counts':False,'elapsed_seconds':out['elapsed_seconds']}))
if __name__=='__main__':main()
