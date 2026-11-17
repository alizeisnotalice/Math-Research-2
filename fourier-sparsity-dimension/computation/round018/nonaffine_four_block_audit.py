#!/usr/bin/env python3
"""Exact audit of the explicit non-affine 4-block common-bent code (R18)."""
import hashlib
import json
import platform
import sys
import time
import argparse
from pathlib import Path

N = 64
FULL = (1 << N) - 1

def fwht(a):
    a = list(a); h = 1
    while h < len(a):
        for i in range(0, len(a), 2*h):
            for j in range(i, i+h):
                a[j], a[j+h] = a[j]+a[j+h], a[j]-a[j+h]
        h *= 2
    return a

def bent(mask):
    return all(abs(x) == 8 for x in fwht([1 - 2*((mask >> x) & 1) for x in range(N)]))

def rank(vs):
    p = {}
    for v in vs:
        while v:
            j = v.bit_length()-1
            if j in p: v ^= p[j]
            else: p[j] = v; break
    return len(p)

def span(bs):
    ans = {0}
    for b in bs: ans |= {x^b for x in tuple(ans)}
    return ans

def three_flats():
    spaces = set()
    for raw in range(1 << 18):
        bs = [(raw >> (6*j)) & 63 for j in range(3)]
        if rank(bs) == 3: spaces.add(tuple(sorted(span(bs))))
    assert len(spaces) == 1395
    flats = {tuple(sorted(a^v for v in S)) for S in spaces for a in range(N)}
    assert len(flats) == 11160
    return [sum(1 << x for x in A) for A in flats]

def affine(mask):
    pts = [x for x in range(N) if mask >> x & 1]
    if not pts: return False
    H = {x^pts[0] for x in pts}
    return len(H) == 1 << rank(H) and span(H) == H

def block(a,b):
    return sum(((((x >> 0)&1) == a and ((x >> 5)&1) ==
                ((((x >> 1)&1)&((x >> 2)&1)) ^ (((x >> 3)&1)&((x >> 4)&1)) ^ b)) << x)
               for x in range(N))

def phase(x):
    q = (((x>>1)&1)&((x>>2)&1)) ^ (((x>>3)&1)&((x>>4)&1))
    return (((x>>1)&1)&((x>>2)&1)) ^ (((x>>0)&1)*(((x>>1)&1)^((x>>3)&1))) ^ (((x>>5)&1)^q)*(((x>>2)&1)^((x>>4)&1))

def maskof(fn): return sum(fn(x) << x for x in range(N))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', type=Path, help='Optional round-018 graph-cases.json to cross-check.')
    ap.add_argument('--output', type=Path,
                    default=Path(__file__).with_name('nonaffine-four-block-audit.json'))
    args = ap.parse_args()
    begun = time.monotonic()
    F = maskof(phase); assert bent(F)
    B = [block(a,b) for a in range(2) for b in range(2)]
    assert all(x.bit_count() == 16 and not affine(x) for x in B)
    assert sum(B) == FULL and rank(B) == 4
    flats = three_flats()
    input_sha = None
    if args.input is not None:
        raw = args.input.read_bytes(); supplied = json.loads(raw)
        input_sha = hashlib.sha256(raw).hexdigest()
        case = next(c for c in supplied['cases'] if c['r'] == 0)
        supplied_blocks = {sum(1 << x for x in A) for A in case['blocks']}
        assert int(case['f0_hex'], 16) == F and supplied_blocks == set(B)
        assert case['dimension'] == 9 and case['phase_count'] == 512
    # Independent exact condition: F plus all block flips and all affine phases.
    coords = [maskof(lambda x,j=j: (x>>j)&1) for j in range(6)]
    gens = B + [FULL] + coords
    # The block algebra already contains 1 and x0 (the sums over b and over
    # a respectively), so its intersection with RM(1,6) has dimension two.
    assert rank(gens) == 9
    code = span(gens)
    assert len(code) == 512 and all(bent(F^v) for v in code)
    partial = []
    legal8 = [A for A in flats if bent(F^A)]
    for idx,A in enumerate(B):
        W = fwht([((1-2*((F>>x)&1)) if A>>x&1 else 0) for x in range(N)])
        support = {z for z,v in enumerate(W) if v}
        a,b = divmod(idx,2)
        predicted = {z for z in range(N) if
                     (((z>>5)&1)==0 and ((z>>3)&1)==a and ((z>>4)&1)==b) or
                     (((z>>5)&1)==1 and ((z>>1)&1)==a and ((z>>2)&1)==b)}
        assert support == predicted and set(W) <= {-8,0,8}
        contained = [C for C in legal8 if not (C & ~A)]
        contained_all = [C for C in flats if not (C & ~A)]
        assert not contained
        assert not contained_all
        partial.append({'a':a,'b':b,'mask_hex':hex(A),'walsh_support':sorted(support),
                        'walsh_values_set':sorted(set(W)),
                        'contained_all_affine_q_flats':len(contained_all),
                        'contained_legal_affine_q_flats':len(contained)})
    assert set().union(*[{z for z,v in enumerate(fwht([((1-2*((F>>x)&1)) if A>>x&1 else 0) for x in range(N)])) if v} for A in B]) == set(range(N))
    # The four supports are a partition, so check disjointness separately.
    SS=[set(p['walsh_support']) for p in partial]
    assert all(not(SS[i]&SS[j]) for i in range(4) for j in range(i))
    result={'schema':'a28-r018-nonaffine-four-block-v1','arithmetic':'exact F2/integer FWHT',
      'phase_mask_hex':hex(F),'phase_bent':True,'block_count':4,'block_weight':16,
      'blocks_affine':False,'blocks_partition_domain':True,'block_span_rank':4,
      'RM1_plus_blocks_rank':9,'all_common_bent_phase_count':len(code),
      'all_common_bent_verified':True,'all_affine_q_flats':len(flats),
      'legal_affine_q_flats_for_F':len(legal8),'each_block_contains_any_affine_q_flat':False,
      'each_block_contains_legal_affine_q_flat':False,
      'dual_supports_disjoint_partition':True,'blocks':partial,
      'input_path':str(args.input) if args.input else None,'input_sha256':input_sha,
      'python':sys.version,'platform':platform.platform(),'elapsed_seconds':time.monotonic()-begun,
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['phase_mask_hex','all_common_bent_phase_count','legal_affine_q_flats_for_F','dual_supports_disjoint_partition','elapsed_seconds']},indent=2))

if __name__ == '__main__': main()
