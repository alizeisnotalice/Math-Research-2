#!/usr/bin/env python3
"""Prepare deterministic 0-1 products for GPU root-group allocation.

float64 only: prepared means numerically well reconstructed, never proved.
An unknown/mismatched input is retained as status=unknown and is not an
exclusion.  No GPU code is launched here.
"""
from __future__ import annotations

import hashlib, json
from pathlib import Path
import numpy as np

SEED=20261115
PAIR_TOL=2e-7
RECON_TOL=2e-7
ROOT_TOL=2e-7
OUT=Path(__file__).with_name('root_group_candidates_20261115.json')

def conv(*ps):
    a=np.array([1],dtype=np.int64)
    for p in ps: a=np.convolve(a,p)
    if np.any((a!=0)&(a!=1)): raise ValueError('support collision')
    return a

def sparse_poly(deg, support):
    a=np.zeros(deg+1,dtype=np.int64); a[list(sorted(set(support)))]=1
    assert a[0]==a[-1]==1
    return a

def groups_from_roots(c):
    roots=np.roots(c[::-1].astype(float))
    used=np.zeros(len(roots),bool); groups=[]; pair_err=0.0
    for i,z in enumerate(roots):
        if used[i]: continue
        if abs(z.imag)<=PAIR_TOL*max(1.0,abs(z)):
            groups.append(np.array([-z.real,1.0])); used[i]=True
        else:
            dist=np.where(used,np.inf,abs(roots-np.conj(z)))
            j=int(np.argmin(dist)); err=float(dist[j]/max(1.0,abs(z)))
            if j==i or err>PAIR_TOL: raise RuntimeError(f'unmatched root {z}, err={err}')
            pair_err=max(pair_err,err)
            groups.append(np.array([abs(z)**2,-2*z.real,1.0])); used[i]=used[j]=True
    groups.sort(key=lambda a:(len(a),tuple(np.round(a,12))))
    return roots,groups,pair_err

def root_residual(c,roots):
    out=0.0
    for z in roots:
        val=np.polynomial.polynomial.polyval(z,c)
        scale=sum(abs(float(a))*abs(z)**i for i,a in enumerate(c))
        out=max(out,float(abs(val)/max(scale,1.0)))
    return out

def make_record(name,kind,c,meta):
    digest=hashlib.sha256(json.dumps(c.tolist(),separators=(',',':')).encode()).hexdigest()
    rec={'name':name,'kind':kind,'degree':int(len(c)-1),'support':[int(i) for i,x in enumerate(c) if x],
         'coefficients_ascending':[int(x) for x in c],'sha256':digest,'meta':meta}
    try:
        roots,gs,perr=groups_from_roots(c)
        rebuilt=np.array([1.0])
        for g in gs: rebuilt=np.polynomial.polynomial.polymul(rebuilt,g)
        recon=float(np.max(abs(rebuilt-c))) if len(rebuilt)==len(c) else float('inf')
        rres=root_residual(c,roots)
        rec.update({'group_count':len(gs),'groups_ascending_coefficients':[[float(x) for x in g] for g in gs],
                    'pairing_relative_residual':perr,'root_relative_residual':rres,
                    'group_product_inf_residual':recon})
        rec['status']='prepared' if len(gs)<=25 and max(perr,rres,recon)<=max(PAIR_TOL,RECON_TOL,ROOT_TOL) else 'unknown'
    except Exception as exc:
        rec.update({'status':'unknown','error':str(exc),'group_count':None})
    return rec

def main():
    global OUT
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument("--out",type=Path,required=True)
    OUT=parser.parse_args().out
    assert Path("/Volumes/data").is_mount()
    rng=np.random.default_rng(SEED); items=[]
    # Calibration: collision-free digit product, total degree 24.
    items.append(('digit_direct_24','digit_direct',conv(np.ones(3,dtype=int),
        np.array([1,0,0,0,1,0,0,0,1]),np.r_[1,np.zeros(13,dtype=int),1]),{'radices':[3,3,2],'steps':[1,4,14]}))
    # Stretched direct sums; stretch_24 is a cyclotomic calibration.
    for name,H,L in [
        ('stretch_24',sparse_poly(4,[0,1,3,4]),5),
        ('stretch_35',sparse_poly(5,[0,1,3,5]),6),
        ('stretch_48',sparse_poly(6,[0,1,2,5,6]),7),
    ]:
        lifted=np.zeros((len(H)-1)*L+1,dtype=int); lifted[::L]=H
        items.append((name,'stretched_direct_sum',conv(H,lifted),{'H_support':[int(i) for i,x in enumerate(H) if x],'L':L}))
    # Fixed-seed general and reciprocal sparse masks.  Keep exact 0-1 C by construction.
    specs=[('random_28',28,8,False),('palindrome_30',30,8,True),('random_34',34,9,False),
           ('palindrome_38',38,10,True),('random_42',42,10,False),('random_46',46,11,False)]
    for name,d,count,pal in specs:
        if pal:
            half=sorted(rng.choice(np.arange(1,d//2),size=(count-2)//2,replace=False).tolist())
            sup=[0,d]+half+[d-i for i in half]
        else:
            sup=[0,d]+rng.choice(np.arange(1,d),size=count-2,replace=False).tolist()
        items.append((name,'palindromic_sparse_mask' if pal else 'seeded_sparse_mask',
                      sparse_poly(d,sup),{'seed':SEED,'requested_terms':count,'palindromic':pal}))
    records=[make_record(n,t,c,m) for n,t,c,m in items]
    OUT.write_text(json.dumps({'seed':SEED,'float_dtype':'float64','pair_tolerance':PAIR_TOL,
        'root_tolerance':ROOT_TOL,'reconstruction_tolerance':RECON_TOL,'candidates':records},ensure_ascii=False,indent=2)+'\n')
    for r in records:
        print(r['name'],r['degree'],r.get('group_count'),r['status'],r.get('group_product_inf_residual'))

if __name__=='__main__': main()
