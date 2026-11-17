#!/usr/bin/env python3
"""Integer controls for the general granularity-two and saturated-residual proofs."""
import argparse, hashlib, json, platform, sys, time
from collections import Counter
from math import comb
from pathlib import Path

def parity(x): return bin(x).count('1') % 2
def fwht(v):
    a=list(v);h=1
    while h<len(a):
        for p in range(0,len(a),2*h):
            for i in range(p,p+h): a[i],a[i+h]=a[i]+a[i+h],a[i]-a[i+h]
        h*=2
    return a
def rank(rows):
    piv={}
    for y in rows:
        while y:
            i=y.bit_length()-1
            if i in piv:y^=piv[i]
            else:piv[i]=y;break
    return len(piv)
def K(n,j,w):
    return sum((-1)**i*comb(w,i)*comb(n-w,j-i)
               for i in range(max(0,j-(n-w)),min(j,w)+1))
def bits_word(values):return sum(int(b)<<i for i,b in enumerate(values))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();start=time.monotonic()
    # Root and the two child selectors are independent bits 0,1,2.
    vals=[]
    for x in range(128):
        root=x&1; child=(x>>(1+root))&1;leaf=2*root+child
        vals.append(1-2*((x>>(3+leaf))&1))
    W=fwht(vals);S=[i for i,w in enumerate(W) if w]
    assert len(S)==16 and {abs(W[i]) for i in S}=={32}
    assert rank(S)==7 and rank([i^S[0] for i in S])==6
    assert K(16,4,6)==K(16,4,10)==-20
    assert [K(16,2,w) for w in [0,6,8,10,16]]==[120,0,-8,0,120]
    # Actual sharp six-dimensional atom code from round010, not a scalar relaxation.
    f=lambda x: ((x&1)*((x>>3)&1))^(((x>>1)&1)*((x>>2)&1))
    v=lambda x: ((x&1)*((x>>2)&1))^(((x>>1)&1)*((x>>2)&1))^(((x>>1)&1)*((x>>3)&1))
    basis=[65535]+[bits_word((x>>j)&1 for x in range(16)) for j in range(4)]+[bits_word(v(x) for x in range(16))]
    assert rank(basis)==6
    code=[]
    for coeff in range(64):
        word=0
        for j,b in enumerate(basis):
            if coeff>>j&1:word^=b
        code.append(word)
    wt=lambda w:bin(w).count('1')
    hist=dict(Counter(wt(w) for w in code));assert hist=={0:1,16:1,8:30,6:16,10:16}
    assert all(parity(x&y)==0 for x in code for y in code)
    assert all(sum((1-2*f(i))*((w>>i)&1) for i in range(16))%4==0 for w in code)
    six=[w for w in code if wt(w)==6]
    assert all(wt(x&y)==(6 if i==j else 2) for i,x in enumerate(six) for j,y in enumerate(six))
    assert rank(six)<=6  # binary rank differs from the positive real Gram rank 16.
    # Multiplication-module failure: all 8 affine (z,w) shifts preserve bentness;
    # multiplying z by 1+w destroys it.
    family=[]
    for coeff in range(8):
        ph=[((x&1)*((x>>1)&1))^(((x>>2)&1)*((x>>3)&1))^
            (coeff&1)^(((coeff>>1)&1)*((x>>3)&1))^(((coeff>>2)&1)*((x>>2)&1)) for x in range(16)]
        family.append({abs(y) for y in fwht([1-2*b for b in ph])}=={4})
    bad=[((x&1)*((x>>1)&1))^(((x>>2)&1)*((x>>3)&1))^(((x>>2)&1)*(1^((x>>3)&1))) for x in range(16)]
    bad_abs=sorted({abs(y) for y in fwht([1-2*b for b in bad])})
    assert all(family) and bad_abs!=[4]
    out={'status':'passed','arithmetic':'exact integers/F2','python':sys.version,'platform':platform.platform(),
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'tree_sharpness':{'support':S,'coefficient_numerators_over_4':[W[i]//32 for i in S],'linear_rank':7,'affine_rank':6,'truth_table_checks':128},
         'K4_weight6_and_10':[-20,-20],'saturated_dim6_atom_code':{'basis':basis,'weights':hist,'self_orthogonal_pairs_checked':4096,'six_weight_gram':'4 I_16 + 2 J_16','real_gram_determinant':36*4**15},
         'module_countercontrol':{'all_8_members_bent':True,'multiplied_direction_walsh_absolute_values':bad_abs},
         'elapsed_seconds':time.monotonic()-start,
         'scope':'Finite arithmetic controls; all-dimension theorems are proved separately.'}
    a.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k in ['status','elapsed_seconds']}))
if __name__=='__main__':main()
