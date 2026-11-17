#!/usr/bin/env python3
"""Compare independently recovered short bins with the earlier full n7 count."""
import argparse,hashlib,json
from fractions import Fraction
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('--short-gaps',type=Path,required=True)
p.add_argument('--full-count',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
s=json.loads(a.short_gaps.read_text()); g=s['global']
ref=json.loads(a.full_count.read_text())
assert s['status']=='passed' and ref['n']==7
A1,A2,B=int(ref['A1']),int(g['A2']),int(g['B'])
assert A2==int(ref['A2'])
h=[int(g['h'][f'h{i}']) for i in (1,2,3)]
assert 4*h[0]-2*B==A1
assert h[1]==2*h[2]
q2=2*h[0]-B
assert 2*q2==A1
tail=B-sum(h);tail_moment=A2-sum((i+1)*x for i,x in enumerate(h))
assert tail==int(g['T']['4']) and 4*tail<=tail_moment<=46*tail
f=lambda x,y:str(Fraction(x,y))
out={'status':'passed','scope':'finite n7 independent full-count crosscheck, not asymptotic',
     'A1':str(A1),'A2':str(A2),'B':str(B),'h':[str(x) for x in h],
     'A1_equals_4h1_minus_2B':True,'h2_equals_2h3':True,
     'immediate_merge_q2_count':str(q2),
     'immediate_merge_original_probability':f(q2,A2),
     'identity_immediate_merge_probability_equals_1_over_2R':f(A1,2*A2),
     'L_at_least_4_B_mass':f(tail,B),
     'L_at_least_4_original_mass':f(tail_moment,A2),
     'conditional_tail_mean_L':f(tail_moment,tail),
     'inputs':{'short_gaps_sha256':hashlib.sha256(a.short_gaps.read_bytes()).hexdigest(),
               'full_count_sha256':hashlib.sha256(a.full_count.read_bytes()).hexdigest()}}
a.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
