#!/usr/bin/env python3
"""Independent integer replay of a d=15 prefix plan; no searcher imports."""
from __future__ import annotations

import argparse
import copy
import hashlib
import itertools
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

PAIR = ((0, 0), (1, 0), (0, 1))


def digest(x: object) -> str: return hashlib.sha256(json.dumps(x, separators=(',', ':')).encode()).hexdigest()


def decode(code: str):
    a, b = [1], [1]
    for c in code:
        x, y = PAIR[int(c)]; a.append(x); b.append(y)
    return a, b


def conv(a, b, end):
    out = []
    for k in range(end+1): out.append(sum(a[i]*b[k-i] for i in range(max(0,k-len(b)+1),min(k,len(a)-1)+1)))
    return out


def check_plan(plan: dict) -> dict:
    assert plan['schema'] == 'd15-necessary-prefix-plan-v1'
    p = plan['parameters']; assert p['d'] == 15 and p['m'] == 64 and p['n'] == 104 and p['minimum_starts'] == 8192
    valid, survivors, rejected = [], [], Counter()
    for chars in itertools.product('012', repeat=14):
        code=''.join(chars);a,b=decode(code)
        if set(conv(a,b,14)) - {0,1} or sum(a[i]*b[15-i] for i in range(1,15)): continue
        valid.append(code); product=conv(a,b,28); reasons=[]
        if any(x>1 for x in product):reasons.append('prefix_product_coefficient_gt_1')
        if any((product[15+i] if 15+i < len(product) else 0)>=1 and a[i]+b[i]>0 for i in range(1,15)):reasons.append('forced_positive_fractional_cross_at_d_plus_i')
        if reasons:rejected.update(reasons)
        else:survivors.append(code)
    assert plan['raw_ternary_count']==3**14 and plan['valid_low_prefix_count']==len(valid) and plan['survivor_count']==len(survivors)
    assert plan['valid_low_prefix_codes_sha256']==digest(valid) and plan['survivor_codes_sha256']==digest(survivors)
    assert plan['strict_filter_rejected_by_reason_counts']==dict(rejected) and plan['survivor_codes']==survivors
    starts=p['coverage_plan_starts'];assigned=plan['coverage_assignment_codes'];assert starts==len(assigned) and starts%8192==0
    expected=[survivors[i%len(survivors)] for i in range(starts)];assert assigned==expected
    counts=Counter(assigned);assert set(counts)==set(survivors) and min(counts.values())>=1
    assert plan['per_prefix_assignment_counts']==[{'code':x,'count':counts[x]} for x in survivors]
    return {'raw':3**14,'valid':len(valid),'survivors':len(survivors),'starts':starts,'rejected':dict(rejected)}


def reject(label, plan):
    try: check_plan(plan)
    except (AssertionError, KeyError, TypeError): return label
    raise AssertionError('production validator accepted '+label)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    assert Path('/Volumes/data').is_mount();started=time.monotonic();plan=json.loads(a.input.read_text());summary=check_plan(plan)
    p1=copy.deepcopy(plan);p1['survivor_codes'].pop()
    p2=copy.deepcopy(plan);p2['coverage_assignment_codes'][1]=p2['coverage_assignment_codes'][2]
    p3=copy.deepcopy(plan);p3['parameters']['coverage_plan_starts']-=8192
    p4=copy.deepcopy(plan);p4['parameters']['d']=14
    out={'status':'passed_independent_d15_integer_replay','summary':summary,
         'negative_tests_rejected':[reject('drop_survivor',p1),reject('wrong_rotation',p2),reject('false_starts',p3),reject('wrong_d',p4)],
         'input_sha256':hashlib.sha256(a.input.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'elapsed_seconds':time.monotonic()-started,'interpreter':sys.version,
         'claim_boundary':'Complete finite integer replay of the stated prefix plan only; no MPS/GPU and no global factorization/exclusion claim.'}
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))

if __name__=='__main__':main()
