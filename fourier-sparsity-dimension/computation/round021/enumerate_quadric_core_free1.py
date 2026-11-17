#!/usr/bin/env python3
"""Exact local 6-point-quadric x F2 test for 0,+/-8 Walsh support 6."""
import json, platform, sys, hashlib, time
from pathlib import Path

def fwht(a):
    a=list(a); d=1
    while d<len(a):
        for i in range(0,len(a),2*d):
            for j in range(i,i+d): a[j],a[j+d]=a[j]+a[j+d],a[j]-a[j+d]
        d*=2
    return a

def main():
    started=time.monotonic()
    # Q=u0u1+u2u3; its 1-set has six points.
    B=[u for u in range(16) if (((u>>0)&1)&((u>>1)&1)) ^ (((u>>2)&1)&((u>>3)&1))]
    assert len(B)==6
    good=0; common_free_character=0; witness=None
    # coordinates are (u,z), with z as bit 4.
    for mask in range(1<<12):
        h=[0]*32
        for j,u in enumerate(B):
            for z in range(2): h[u|(z<<4)]= -1 if ((mask>>(2*j+z))&1) else 1
        W=fwht(h)
        if sum(w!=0 for w in W)!=6 or any(abs(w)!=8 for w in W if w): continue
        good+=1
        ratios=[]
        for j,u in enumerate(B): ratios.append(h[u]*h[u|16])
        # ratio is independent of u iff sign is r(u)*chi_b(z).
        if len(set(ratios))==1: common_free_character+=1
        elif witness is None:
            witness={'mask_hex':hex(mask),'ratios_z1_over_z0':ratios,'walsh_support':[i for i,w in enumerate(W) if w],'walsh_values':[w for w in W if w]}
    out={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-started,'python':sys.version,'platform':platform.platform(),'B':B,'sign_tables':1<<12,'local_plateaued_tables':good,'common_free_character_tables':common_free_character,'non_affine_free_fiber_witness':witness}
    json.dump(out,Path(__file__).with_name('quadric_core_free1.json').open('w'),indent=2,sort_keys=True);print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
