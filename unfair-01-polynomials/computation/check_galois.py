#!/usr/bin/env python3
import json
from algebraic_verify import RealField

f=RealField([1,-3,1],[0,1])
alpha=f.norm([0,1]); conjugate=f.norm([3,-1])
assert f.sign(alpha)==1
assert f.sign(f.norm([1,-1]))==1
assert f.sign(f.norm([2,-1]))==1
assert f.mul(alpha,conjugate)==f.one
assert f.norm([1,-3,1])==f.zero
print(json.dumps({'minpoly':[1,-3,1],'chosen_interval':[0,1],
                  'alpha_between_zero_and_one':True,'other_conjugate_above_one':True,
                  'norm_exactly_one':True,'is_original_problem_counterexample':False},indent=2))
