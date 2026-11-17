#!/usr/bin/env python3
"""Independent fiber-count histogram and exact coordinate-span certificates."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path
import platform
import sys
import time


def rank(vectors):
    p={}
    for x in vectors:
        while x:
            b=x.bit_length()-1
            if b in p:x^=p[b]
            else:p[b]=x;break
    return len(p)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--receipt',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    start=time.monotonic();receipt=json.loads(a.receipt.read_text());r=receipt['result']
    assert receipt.get('execution_exit_code',receipt.get('ssh_exit_code'))==0 and r['status']=='passed'
    here=Path(__file__).resolve().parent
    assert r['source_sha256']==hashlib.sha256((here/'triangular_weight_scan.cpp').read_bytes()).hexdigest()
    shards=sorted(r['shards'],key=lambda x:x['result']['start']);at=0
    hist=[0]*257;activity=[0]*257
    for shard in shards:
        s=shard['result'];assert shard['exit_code']==0 and s['status']=='passed'
        assert s['start']==at and s['k']==4 and s['q']==16 and s['dimV']==31
        assert len(s['histogram'])==len(s['parameter_activity_by_weight'])==257
        assert sum(s['histogram'])==s['count'];at+=s['count']
        for w in range(257):
            hist[w]+=s['histogram'][w];activity[w]|=s['parameter_activity_by_weight'][w]
    assert at==sum(hist)==1<<31
    k=4;q=16;slots=[(j,p) for j in range(k) for p in range(1<<j)]
    cylinders=[sum(1<<v for v in range(q) if (v&((1<<j)-1))==p) for j,p in slots]
    active=[0]*(1<<len(slots));by_size=Counter()
    for z in range(1<<len(slots)):
        if z:
            low=z&-z;active[z]=active[z^low]|cylinders[low.bit_length()-1]
        by_size[active[z].bit_count()]+=1
        image=[]
        for v in range(q):
            b=v
            for i,(j,p) in enumerate(slots):
                if z>>i&1 and (v&((1<<j)-1))==p:b^=1<<j
            image.append(b)
        assert len(set(image))==q
    independent=[0]*257
    for occupied,number in by_size.items():
        for ones in range(q-occupied+1):
            w=q*occupied//2+q*ones
            independent[w]+=number*(1<<occupied)*comb(q-occupied,ones)
    assert hist==independent and hist[q]==3*q
    C=[(1<<q)-1]+[sum(((v>>j)&1)<<v for v in range(q)) for j in range(k)]
    for j in range(k):C.append(sum(1<<(q+i) for i,(jj,p) in enumerate(slots) if jj==j))
    assert rank(C)==2*k+1
    basis_weights=[q]*q+[q*(1<<(k-j-1)) for j,p in slots]
    rows=[]
    for L in [Fraction(1),Fraction(3,2),Fraction(2),Fraction(3),Fraction(4),Fraction(7),Fraction(8)]:
        bound=int(L*q);observed=0
        for w in range(bound+1):observed|=activity[w]
        allowed=sum(1<<j for j,w in enumerate(basis_weights) if w<=bound)
        # Observed OR gives an upper bound on the span. Individually eligible
        # coordinate basis vectors supply the matching lower bound exactly.
        assert observed==allowed
        low_basis=[1<<j for j,w in enumerate(basis_weights) if w<=bound]
        J=max(0,k-1-(L.numerator//L.denominator).bit_length()+1)
        assert len(low_basis)==2*q-(1<<J)
        residual=31-rank(low_basis+C)
        assert residual==(1<<J)-1-J
        rows.append({'L':str(L),'J':J,'low_span_rank':len(low_basis),
                     'V_mod_C_plus_low':residual,'activity_mask_hex':hex(observed),
                     'eligible_words_including_zero':sum(hist[:bound+1])})
    result={'status':'passed','arithmetic':'integers, F2 and Fraction','tolerance':0,
            'complete_codewords':at,'independent_triangular_permutations':len(active),
            'minimum_words':hist[q],'histogram_nonzero':{str(w):n for w,n in enumerate(hist) if n},
            'slope_support_histogram':dict(sorted(by_size.items())),'thresholds':rows,
            'elapsed_seconds':time.monotonic()-start,'python':sys.version,'platform':platform.platform(),
            'receipt_sha256':hashlib.sha256(a.receipt.read_bytes()).hexdigest(),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','complete_codewords':at,'thresholds':rows}))


if __name__=='__main__':main()
