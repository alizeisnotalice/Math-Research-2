#!/usr/bin/env python3
"""Generate all 80 canonical necessary H9(3,5) real root layouts."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
PRE=HERE/'H9-cycle-prerequisites.json'
FULL=(1<<37)-1

# Stable 37-bit schema: cycle 0..15; then G5 mother 0 (real), 1, 2.
FAMILIES=((16,2,18,3),(21,3,24,5),(29,3,32,5)) # x3 start/count, x5 start/count


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def selected(z): return [i for i in range(37) if z>>i&1]
def fill(start,count): return ((1<<count)-1)<<start


def g5_word(eps,eta):
    z=0
    for i,(as_,ac,bs,bc) in enumerate(FAMILIES):
        assert eps[i]+eta[i] in (0,1,2)
        if eps[i]: z|=fill(as_,ac)
        if eta[i]: z|=fill(bs,bc)
    return z


def decode_g5(z):
    eps=[];eta=[]
    for as_,ac,bs,bc in FAMILIES:
        a=(z>>as_)&((1<<ac)-1);b=(z>>bs)&((1<<bc)-1)
        assert a in (0,(1<<ac)-1) and b in (0,(1<<bc)-1)
        eps.append(int(bool(a)));eta.append(int(bool(b)))
    return eps,eta


def record(kind,z):
    cycle=z&0xffff; g=z>>16; eps,eta=decode_g5(g<<16)
    assert cycle in ZERO and len(set(a+b for a,b in zip(eps,eta)))==1
    return {'class_kind':kind,'mask_u37':z,'mask_hex':f'0x{z:010x}','cycle_mask_u16':cycle,
            'epsilon_x3_by_G5_mother':eps,'eta_x5_by_G5_mother':eta,
            'common_e':eps[0]+eta[0],'selected_group_indices':selected(z),
            'complement_mask_u37':z^FULL,'complement_hex':f'0x{z^FULL:010x}'}


def main():
    p=json.loads(PRE.read_text()); global ZERO
    ZERO=set(p['zero_pair_sum_masks_exhaustive_2pow16']); assert len(ZERO)==16
    # e=0/e=2 are one complement family: represent by e=0 (all G5 roots in B).
    layouts=[]
    for c in sorted(ZERO):
        z=c; assert z<z^FULL
        layouts.append(record('e0_e2',z))
    # e=1: use A/B exchange to fix the negative-real G5 mother (block 0)
    # in the x3 fibre.  The other two conjugate mother pairs remain free.
    raw=[]
    for c in sorted(ZERO):
        for bits in range(4):
            eps=[1,bits&1,(bits>>1)&1]; eta=[0,1-eps[1],1-eps[2]]
            raw.append(c|g5_word(eps,eta))
    assert len(raw)==64 and len(set(raw))==64 and all(decode_g5(z)[0][0]==1 for z in raw)
    assert not {z^FULL for z in raw}&set(raw)
    layouts.extend(record('e1',z) for z in sorted(raw))
    assert len(layouts)==80 and len({r['mask_u37'] for r in layouts})==80
    assert {r['common_e'] for r in layouts if r['class_kind']=='e0_e2'}=={0}
    assert {r['common_e'] for r in layouts if r['class_kind']=='e1'}=={1}
    # Standard pair is present before quotient: H9(x^3) has cycle Phi5*Phi15 mask 63 and eps=(1,1,1).
    standard_x3=63|g5_word([1,1,1],[0,0,0]);standard_x5=standard_x3^FULL
    assert standard_x3 in {r['mask_u37'] for r in layouts}
    assert 0 in {r['mask_u37'] for r in layouts} # constant/full pair.
    schema=[]
    for n,j,k in p['cycle_groups_n_j_k']: schema.append(f'cycle:Phi{n}:pair_j{j}:zeta75_power{k}')
    schema += ['G5:mother0:x3_linear','G5:mother0:x3_quadratic','G5:mother0:x5_linear','G5:mother0:x5_quadratic0','G5:mother0:x5_quadratic1']
    for block in (1,2):
        schema += [f'G5:mother{block}:x3_quadratic{i}' for i in range(3)]
        schema += [f'G5:mother{block}:x5_quadratic{i}' for i in range(5)]
    assert len(schema)==37
    out={'format':'H9-a3b5-real-necessary-layout-v1','prerequisite_path':PRE.name,'prerequisite_sha256':sha(PRE),
         'group_bit_order':schema,'group_counts':{'cyclotomic_quadratic':16,'G5_quadratic':19,'G5_linear':2,'total_bits':37},
         'full_fibre_necessity_dependency':'Uses the already supplied/revised two-fibre root-product lemma: for every G5 mother root, x3/x5 fibre selections are whole fibres and eps_i+eta_i is the same e in {0,1,2}. This generator records its necessary consequences but does not reprove that lemma.',
         'cycle_trace_necessity_dependency':'Uses the independently verified Phi75 zero-pair-sum mask set in the prerequisite file.',
         'exchange_normalization':{'involution':'bitwise 37-bit complement exchanges A and B','e0_e2':'choose e=0 representatives','e1':'fix G5 negative-real mother x3 fibre (epsilon_0=1), leaving two conjugate mother-pair fibre bits free'},
         'coverage':{'e0_e2_raw':32,'e0_e2_canonical':16,'e1_raw':128,'e1_canonical':64,'canonical_total':80,
                     'constant_full_pair_present':True,'standard_pair_masks_u37':[standard_x3,standard_x5]},
         'canonical_layouts':layouts,
         'claim_boundary':'This is an exact necessary root-layout coverage ledger. It does not test coefficients, exclude layouts, or prove a counterexample.',
         'first_unproved_step':'For each nonstandard one of these 80 necessary classes, determine whether its associated real factor pair can have all coefficients in [0,1], and either exclude it rigorously or produce an exact certified counterexample.'}
    path=HERE/'H9-real-necessary-layouts.json';path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'path':str(path),'sha256':sha(path),'coverage':out['coverage']},indent=2))

if __name__=='__main__': main()
