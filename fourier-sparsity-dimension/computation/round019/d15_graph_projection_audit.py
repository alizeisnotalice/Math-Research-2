#!/usr/bin/env python3
"""Exact R19 scan: projection graphs among weight-16 directions in 20 V15 codes."""
import argparse
import hashlib
import json
import platform
import sys
import time
import random
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
    return all(abs(v) == 8 for v in fwht([1-2*((mask >> x)&1) for x in range(N)]))

def rank(vectors):
    pivots = {}
    for x in vectors:
        while x:
            p = x.bit_length()-1
            if p in pivots: x ^= pivots[p]
            else: pivots[p] = x; break
    return len(pivots)

def span(basis):
    ans = {0}
    for b in basis: ans |= {x^b for x in tuple(ans)}
    return ans

def kernel_two_spaces():
    spaces = set()
    for a in range(1,64):
        for b in range(a+1,64):
            if rank([a,b]) == 2: spaces.add(tuple(sorted(span([a,b]))))
    assert len(spaces) == 651
    return [sum(1 << x for x in S) for S in spaces]

def difference_mask(A):
    points = [x for x in range(N) if A >> x & 1]
    differences = {x ^ y for x in points for y in points}
    # The set is essential: ordinary integer summation over repeated
    # differences creates carries and is not an F2/subset bit mask.
    D = sum(1 << d for d in differences)
    direct_or = 0
    for d in differences:
        direct_or |= 1 << d
    assert D == direct_or and (D & 1) and D.bit_length() <= N
    return D

def graph_kernel_count(A, kernels):
    D = difference_mask(A)
    return sum(not (D & (K & ~1)) for K in kernels)

def direct_kernel_count(A, kernels):
    """Same predicate without a precomputed difference set (control only)."""
    points = [x for x in range(N) if A >> x & 1]
    out = 0
    for K in kernels:
        good = True
        for k in range(1, N):
            if K >> k & 1 and any((A >> (x ^ k)) & 1 for x in points):
                good = False
                break
        out += good
    return out

def coordinate_masks():
    return [sum(((x>>j)&1) << x for x in range(N)) for j in range(6)]

def seed_mask(row):
    return sum(((((x&7)&row['permutation'][x>>3]).bit_count()&1) ^
                ((row['h_bits']>>(x>>3))&1)) << x for x in range(N))

def code_basis(relative_generators):
    basis = [FULL] + coordinate_masks() + [int(x) for x in relative_generators]
    assert len(basis) == rank(basis) == 15
    return basis

def r18_blocks():
    ans=[]
    for a in range(2):
        for b in range(2):
            A=0
            for x in range(N):
                q=(((x>>1)&1)&((x>>2)&1)) ^ (((x>>3)&1)&((x>>4)&1))
                if ((x>>0)&1)==a and ((x>>5)&1)==(q^b): A|=1<<x
            ans.append(A)
    return ans

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--partitions',type=Path,required=True)
    ap.add_argument('--audit', type=Path, help='round-006 common-mm-audit.json; default is beside partitions')
    ap.add_argument('--output',type=Path,default=Path(__file__).with_name('d15-graph-projection-audit.json'))
    args=ap.parse_args(); begun=time.monotonic()
    raw=args.partitions.read_bytes(); parts=json.loads(raw)
    assert len(parts['rows']) == 20
    audit_path = args.audit or args.partitions.with_name('common-mm-audit.json')
    audit_raw = audit_path.read_bytes(); audit = json.loads(audit_raw)
    lookup = {(r['source_file'], r['row_index']): r for r in audit['rows']}
    assert len(lookup) == 20
    kernels=kernel_two_spaces()
    r18=[{'mask_hex':hex(A),'weight':A.bit_count(),'graph_kernels':graph_kernel_count(A,kernels)} for A in r18_blocks()]
    assert all(r['graph_kernels'] > 0 for r in r18)
    assert [r['graph_kernels'] for r in r18] == [direct_kernel_count(A,kernels) for A in r18_blocks()]
    # A deterministic negative control for the predicate itself, not a bent switch.
    negative_control = None
    for seed in range(4096):
        rng = random.Random(seed)
        A = sum(1 << x for x in rng.sample(range(N), 16))
        by_difference = graph_kernel_count(A, kernels)
        by_direct = direct_kernel_count(A, kernels)
        assert by_difference == by_direct
        if by_difference == 0:
            negative_control = {'seed':seed, 'mask_hex':hex(A), 'weight':A.bit_count(),
                                'graph_kernels':0, 'direct_graph_kernels':0}
            break
    assert negative_control is not None
    total_w16=total_unique=non_graph=0; seen=set(); per_code=[]; candidates=[]
    for ci,item in enumerate(parts['rows']):
        source=Path(item['source_file'])
        audited = lookup[(str(source), item['row_index'])]
        # The source path is provenance only: R6 already stores the seed fields.
        # This avoids depending on the original absolute checkout location.
        f=seed_mask(audited); assert bent(f)
        basis=code_basis(audited['relative_generators']); w16=[]
        for coeff in range(1<<15):
            v=0
            for j,b in enumerate(basis):
                if coeff>>j&1: v^=b
            if v.bit_count()==16: w16.append((v,coeff))
        local_bad=0; local_good=0
        for v,coeff in w16:
            gk=graph_kernel_count(v,kernels)
            if v not in seen:
                seen.add(v)
                total_unique += 1
            if gk: local_good += 1
            else:
                local_bad += 1
                if len(candidates)<16:
                    assert bent(f^v)
                    candidates.append({'code_index':ci,'source_file':str(source),'row_index':item['row_index'],
                      'seed_mask_hex':hex(f),'switch_mask_hex':hex(v),'code_coefficient_hex':hex(coeff),
                      'seed_walsh_abs_set':sorted(set(abs(x) for x in fwht([1-2*((f>>z)&1) for z in range(N)]))),
                      'flipped_walsh_abs_set':sorted(set(abs(x) for x in fwht([1-2*(((f^v)>>z)&1) for z in range(N)]))),
                      'difference_mask_hex':hex(difference_mask(v)),'all_651_kernels_fail':True})
        total_w16 += len(w16); non_graph += local_bad
        per_code.append({'code_index':ci,'source_file':str(source),'row_index':item['row_index'],
                         'seed_mask_hex':hex(f),'weight16_codewords':len(w16),
                         'graph_projection_words':local_good,'nongraph_words':local_bad})
    result={'schema':'a28-r019-d15-graph-projection-v2','scope':'all 20 archived dim-15 common-bent codes; directions of weight16 only',
      'arithmetic':'exact F2 difference sets and integer FWHT on reported candidates','tolerance':0,
      'supersedes':'v1 is invalid: it summed repeated difference bits as integers, allowing carries; discard its numerical results.',
      'partitions_path':str(args.partitions),'partitions_sha256':hashlib.sha256(raw).hexdigest(),
      'audit_path':str(audit_path),'audit_sha256':hashlib.sha256(audit_raw).hexdigest(),
      'two_dimensional_kernels':len(kernels),'r18_positive_control':r18,'negative_control':negative_control,
      'weight16_occurrences':total_w16,'unique_weight16_switches':total_unique,
      'nongraph_occurrences':non_graph,'nongraph_certificates':candidates,'per_code':per_code,
      'python':sys.version,'platform':platform.platform(),'elapsed_seconds':time.monotonic()-begun,
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['two_dimensional_kernels','weight16_occurrences','unique_weight16_switches','nongraph_occurrences','r18_positive_control']},indent=2))

if __name__=='__main__': main()
