#!/usr/bin/env python3
"""Exact auxiliary counterexample; not a counterexample to uniform sums."""
import json, hashlib
from pathlib import Path
from collections import Counter
from fractions import Fraction
u=[0,1,2,3,4];v=[0,6,12,18,24]
c=Counter(i+j for i in u for j in v)
assert len(c)==25 and set(c.values())=={1}
k=[c[i] for i in range(max(c)+1)]
q=[1,2,3,4,5,4,3,2,1] # Phi5 squared
r=list(map(Fraction,k)); quotient=[Fraction(0)]*(len(r)-len(q)+1)
for j in range(len(quotient)-1,-1,-1):
 t=r[j+len(q)-1];quotient[j]=t
 for i,b in enumerate(q):r[j+i]-=t*b
assert not any(r)
assert all(x.denominator==1 for x in quotient)
print(json.dumps(dict(status="passed_exact_phi5_squared_binary_multiple",binary_support=sorted(c),quotient=[int(x) for x in quotient],source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope="Refutes only the auxiliary all-quartic square-free assertion; the exhibited factors are fair."),sort_keys=True))
