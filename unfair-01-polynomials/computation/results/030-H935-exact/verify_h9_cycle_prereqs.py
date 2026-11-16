#!/usr/bin/env python3
"""Exact H9(3,5) prerequisites and the Phi75 zero-sum cycle certificate."""
from __future__ import annotations
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sympy as sp

HERE=Path(__file__).resolve().parent; x=sp.symbols('x')
H=1+x+x**2+x**8+x**9; G=x**5-x**3+1


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def factor_degrees(f,p): return sorted(q.degree() for q,m in sp.factor_list(sp.Poly(f,x,modulus=p))[1] for _ in range(m))
def mask_support(m): return [i for i in range(16) if m>>i&1]


def main():
    phi5,phi15,phi25,phi75=(sp.cyclotomic_poly(n,x) for n in (5,15,25,75))
    assert sp.expand(H-phi5*G)==0
    assert sp.expand(phi5.subs(x,x**3)-phi5*phi15)==0
    assert sp.expand(phi5.subs(x,x**5)-phi25)==0
    C=sp.Poly(sp.expand(H.subs(x,x**3)*H.subs(x,x**5)),x)
    assert C.degree()==72 and sum(C.nth(i) for i in range(73))==25 and all(C.nth(i) in (0,1) for i in range(73))
    # Sturm is exact: one negative real root, hence two nonreal conjugate pairs.
    PG=sp.Poly(G,x)
    sturm={'all_real':int(PG.count_roots(-sp.oo,sp.oo)),'negative_real':int(PG.count_roots(-sp.oo,0)),'positive_real':int(PG.count_roots(0,sp.oo))}
    assert sturm=={'all_real':1,'negative_real':1,'positive_real':0}
    galois_mod={str(p):factor_degrees(G,p) for p in (2,3,5)}
    assert galois_mod=={'2':[5],'3':[2,3],'5':[1,4]}
    # Cyclotomic pairs: zeta_n^j+zeta_n^-j = zeta_75^k+zeta_75^-k.
    groups=[]
    for n in (5,15,25):
        for j in range(1,(n-1)//2+1):
            if sp.gcd(j,n)==1: groups.append((n,j,75*j//n))
    assert groups==[(5,1,15),(5,2,30),(15,1,5),(15,2,10),(15,4,20),(15,7,35),(25,1,3),(25,2,6),(25,3,9),(25,4,12),(25,6,18),(25,7,21),(25,8,24),(25,9,27),(25,11,33),(25,12,36)]
    columns=[]
    for _,_,k in groups:
        rem=sp.Poly(sp.rem(x**k+x**(75-k),phi75,x),x)
        columns.append([int(rem.nth(i)) for i in range(40)])
    M=sp.Matrix(40,16,lambda r,c:columns[c][r]); assert M.rank()==12
    pivcols=list(M.rref()[1]); assert len(pivcols)==12
    pivrows=list(M[:,pivcols].T.rref()[1]); minor=M.extract(pivrows,pivcols).det(); assert minor!=0
    block_masks=[22,41,26176,39296]
    vectors=[sp.Matrix([(m>>i)&1 for i in range(16)]) for m in block_masks]
    assert all(M*v==sp.zeros(40,1) for v in vectors)
    supports=[mask_support(m) for m in block_masks]
    assert set().union(*map(set,supports))==set(range(16)) and sum(map(len,supports))==16
    assert all(not(set(a)&set(b)) for a,b in combinations(supports,2))
    zeros=[]
    for m in range(1<<16):
        if M*sp.Matrix([(m>>i)&1 for i in range(16)])==sp.zeros(40,1): zeros.append(m)
    subsets=sorted(sum(block_masks[i] for i in range(4) if q>>i&1) for q in range(16))
    assert zeros==subsets and len(zeros)==16
    index_by_pair={k:i for i,(_,_,k) in enumerate(groups)}
    lift_blocks=[]
    expected_lifts={(1,3):22,(2,3):41,(1,5):26176,(2,5):39296}
    for (j,p), want in expected_lifts.items():
        roots=[q for q in range(1,75) if (p*q-15*j)%75==0]
        pair_indices=sorted({index_by_pair[min(q,75-q)] for q in roots})
        got=sum(1<<i for i in pair_indices)
        assert got==want
        lift_blocks.append({'Phi5_parent_j':j,'stretch_power':p,'root_exponents_mod_75':roots,
                            'conjugate_pair_indices':pair_indices,'mask_u16':got})
    assert [z['mask_u16'] for z in lift_blocks]==block_masks
    out={'format':'H9-a3b5-prereq-cycle-trace-v1','H9_factorization':{'Phi5':str(phi5),'G5':str(G),'identity_verified':True},
         'target_binary_25_term_verified':True,'G5':{'sturm_exact':sturm,'modular_factor_degrees':galois_mod,
         'root_schema':'mother_0=unique negative real root; mother_1,mother_2=upper-half-plane pairs ordered by increasing real part (to be instantiated by root disks, never raw numerical order).'},
         'cycle_groups_n_j_k':groups,'Phi75_modulus_ascending':[int(sp.Poly(phi75,x).nth(i)) for i in range(41)],
         'pair_sum_remainder_matrix_40x16':[[int(M[r,c]) for c in range(16)] for r in range(40)],
         'matrix_certificate':{'rank':12,'nonzero_12_minor':{'rows':pivrows,'columns':pivcols,'determinant':int(minor)},
         'kernel_block_masks':block_masks,'kernel_block_supports':supports,'disjoint_supports_partition_all_16':True,
         'argument':'Rank 12 and four independent disjoint-support kernel indicators imply every 0/1 kernel vector is a union of these four blocks.'},
         'Phi5_lift_block_binding':{'identity':'For Phi5 parent pair j and stretch p, solve p*q = 15*j (mod 75), then pair q with -q.',
                                    'blocks':lift_blocks,
                                    'mother_real_quadratics':'Phi5 has two real quadratic parent factors x^2-t_j*x+1 with t_j=(-1 plus_or_minus sqrt(5))/2; each listed block is exactly one full p-fibre of one parent pair.'},
         'zero_pair_sum_masks_exhaustive_2pow16':zeros,'zero_pair_sum_count':len(zeros),
         'claim_boundary':'This proves the cyclotomic selected-root sum is zero exactly for these 16 masks. It does not prove coefficient nonnegativity.'}
    path=HERE/'H9-cycle-prerequisites.json';path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'path':str(path),'sha256':sha(path),'rank':12,'zero_masks':zeros},indent=2))

if __name__=='__main__':main()
