#!/usr/bin/env python3
"""Continuous parameter rejection via monotone one-point Green bounds; 128-bit Arb."""
import os,sys,json,hashlib,time
from fractions import Fraction as Q
from pathlib import Path
import flint
from flint import arb,ctx
ctx.prec=128
ROOT=Path(__file__).resolve().parents[1]
TH=Q(18,5)

def A(q):return arb(q.numerator)/q.denominator

def poslog(x):return x.log().max(arb(0))

def logs(s,b):
    u=1/s;au=A(u);ab=A(b)
    plus=arb(0) if u>=1+b else poslog((1+(1+4*ab*au).sqrt())/(2*au))
    minus=arb(0) if u>=1-b else poslog((1+(1-4*ab*au).sqrt())/(2*au))
    d=(1+16*ab*ab*au*au).sqrt()
    imag=poslog((1+d+2*((d+1)/2).sqrt()).sqrt()/(2*au))
    return plus,minus,imag

def bound(box):
    sl,sh,bl,bh=box
    lp,_,li=logs(sh,bh)
    _,lm,_=logs(sh,bl)
    denom=lp+lm+2*li
    if denom==0:return None
    return arb.pi()*A(sl)**2*(1+A(bl)**2)/denom

if __name__=='__main__':
    assert os.path.ismount('/Volumes/data')
    st=time.monotonic()
    tail_small=9*arb.pi()/(64*(arb(9)/8).log())
    tail_large=49*arb.pi()/(16*(arb(21)/4).log())
    assert tail_small>A(TH) and tail_large>A(TH)
    queue=[(Q(3,4),Q(7,2),Q(0),Q(1,2))];leaves=[];visited=0
    while queue:
        box=queue.pop();visited+=1
        if visited>20000:raise RuntimeError('bounded cell budget exhausted')
        v=bound(box)
        if v is None or v>A(TH):
            leaves.append({'box':[str(x) for x in box],'status':'infeasible_zero_vertex' if v is None else 'cost_gt_18_over_5','lower_enclosure':None if v is None else str(v)})
            continue
        sl,sh,bl,bh=box
        if (sh-sl)/Q(11,4)>=(bh-bl)/Q(1,2):
            mid=(sl+sh)/2;queue.extend([(sl,mid,bl,bh),(mid,sh,bl,bh)])
        else:
            mid=(bl+bh)/2;queue.extend([(sl,sh,bl,mid),(sl,sh,mid,bh)])
    out={'status':'PASS_continuous_parameter_vertex_cost_floor','scope':'s>0,0<=b<=1/2; cost floor >18/5 for D4-averaged F(w)=s(w+bw^2) killedGreen majorants. Proof uses radial nesting and b monotonicity described separately; no D boundary grid is a proof.','python':sys.version,'flint':flint.__version__,'precision_bits':128,'workers':1,'threshold':'18/5','small_s':'s<=3/4; infeasible if s<=2/3, otherwise disk-containment bound','small_tail_bound':str(tail_small),'large_s':'s>=7/2; disk-containment bound','large_tail_bound':str(tail_large),'initial_box':['3/4','7/2','0','1/2'],'visited':visited,'leaf_count':len(leaves),'leaves':leaves,'elapsed_seconds':time.monotonic()-st,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dest=ROOT/'computation/results/asymmetric-green-average-vertex-certificate.json';dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='leaves'},indent=2))
