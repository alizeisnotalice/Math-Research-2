#!/usr/bin/env python3
"""Independent exact verifier for a round005 matching certificate."""
import argparse, hashlib, json
from pathlib import Path

def parity(x): return x.bit_count() & 1
def q(x, profile):
    u, v = x & 7, x >> 3
    if profile == 'base': return parity(u & v)
    if profile == 'mm-nonlinear':
        pi = (0,1,3,2,6,7,5,4)
        return parity(u & pi[v]) ^ int(((v & 3) == 3) ^ ((v >> 2) & 1))
    raise AssertionError('unknown seed profile')
def rank(rows):
    piv = {}
    for x in rows:
        while x:
            b=x.bit_length()-1
            if b in piv: x ^= piv[b]
            else: piv[b]=x; break
    return len(piv)

p=argparse.ArgumentParser(); p.add_argument('certificate',type=Path); p.add_argument('--out',required=True,type=Path); a=p.parse_args()
raw=a.certificate.read_bytes(); d=json.loads(raw); matching=d['matching_by_delta']
profile=d['seed_profile']
if len(matching)!=63 or any(len(x)!=16 for x in matching): raise AssertionError('shape')
rows=[]; legal=[]
for delta, pairs in enumerate(matching,1):
    usedp=set(); usedn=set(); plus=minus=0
    for x,y,z,w in pairs:
        if x^y != delta or z^w != delta or x>=y or z>=w: raise AssertionError(('pair',delta,x,y,z,w))
        if q(x,profile)^q(y,profile) or not(q(z,profile)^q(w,profile)): raise AssertionError(('sign',delta,x,y,z,w))
        usedp.add((x,y)); usedn.add((z,w)); plus += 1; minus += 1; rows.append((1<<x)^(1<<y)^(1<<z)^(1<<w))
    if len(usedp)!=16 or len(usedn)!=16: raise AssertionError(('bijection',delta))
    legal.append({'delta':delta,'positive_pairs':plus,'negative_pairs':minus})
R=rank(rows)
if R != d['best_rank_R'] or 64-R != d['best_dimV'] or 63-R != d['best_graph_affine_dimension']: raise AssertionError(('rank',R))
stored_rows=[int(x,16) for x in d['incidence_rows_hex']]
if stored_rows != rows: raise AssertionError('stored incidence rows disagree')
basis=[int(x,16) for x in d['nullspace_basis_hex']]
if len(basis)!=64-R or d['nullspace_basis_dimension'] != len(basis) or rank(basis)!=len(basis): raise AssertionError('basis rank')
for h in basis:
    if any(parity(h&r) for r in rows): raise AssertionError('basis not null')
    for delta in range(1,64):
        corr=sum(-1 if q(x,profile)^q(x^delta,profile)^((h>>x)&1)^((h>>(x^delta))&1) else 1 for x in range(64))
        if corr: raise AssertionError(('corr',delta,corr))
    for t in range(64):
        w=sum(-1 if q(x,profile)^((h>>x)&1)^parity(x&t) else 1 for x in range(64))
        if abs(w)!=8: raise AssertionError(('walsh',t,w))
for mask in range(1 << len(basis)):
    h=0
    for j,b in enumerate(basis):
        if mask >> j & 1: h ^= b
    for t in range(64):
        w=sum(-1 if q(x,profile)^((h>>x)&1)^parity(x&t) else 1 for x in range(64))
        if abs(w)!=8: raise AssertionError(('combined-walsh',mask,t,w))
C=[int(x,16) for x in d['C_basis_hex']]
if len(C)!=7 or rank(C)!=7 or d['C_rank']!=7: raise AssertionError('C rank')
if any(parity(c&r) for c in C for r in rows): raise AssertionError('C containment')
if len({1 | (x << 1) for x in range(64)}) != 64: raise AssertionError('frequency columns')
a.out.write_text(json.dumps({'status':'passed','certificate_sha256':hashlib.sha256(raw).hexdigest(),'seed_profile':profile,'rank_R':R,'dimV':64-R,'graph_affine_dimension':63-R,'basis_vectors_exactly_verified':len(basis),'kernel_masks_fwht_exactly_verified':1 << len(basis),'all_delta_pair_counts':legal},indent=2)+'\n')
