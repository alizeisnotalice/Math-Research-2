#!/usr/bin/env python3
"""Exact SOS coefficients and interval scalar constants; no domain discretization."""
from pathlib import Path
from fractions import Fraction as Q
import math,json,sys,hashlib,os
import flint
from flint import arb,ctx
ctx.prec=128
ROOT=Path(__file__).resolve().parents[1]

def square(a):
    out=[Q(0)]*(2*len(a)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(a):out[i+j]+=x*y
    return out

def run():
    assert os.path.ismount('/Volumes/data')
    a=square([-Q(1,3),Q(0),Q(1)])
    b=[Q(5,3)*x for x in square([-Q(18,25),Q(1)])]
    rhs=[a[i]+(b[i] if i<len(b) else 0)+(Q(28,1125) if i==0 else 0) for i in range(5)]
    lhs=[Q(1),-Q(12,5),Q(1),Q(0),Q(1)]
    assert rhs==lhs
    const=6*arb.pi()/5
    assert const>arb(18)/5
    def cost(q):return math.pi*(1-q**8)/(4*q*q*(-math.log(q)))
    lo,hi=.05,.95;g=(math.sqrt(5)-1)/2
    for _ in range(100):
        a=hi-g*(hi-lo);b=lo+g*(hi-lo)
        if cost(a)>cost(b):lo=a
        else:hi=b
    q=(lo+hi)/2
    proof=ROOT/'proofs/killed-green-dihedral-domain-obstruction.md'
    result={'status':'PASS_exact_SOS_and_Arb_scalar_checks','python':sys.version,'python_flint':flint.__version__,'workers':1,'random_seed':None,'arb_precision_bits':128,'SOS_lhs_coefficients_low_to_high':[str(x) for x in lhs],'SOS_rhs_coefficients_low_to_high':[str(x) for x in rhs],'positive_remainder':'28/1125','cost_floor_arb':str(const),'strictly_greater_than_18_over_5':True,'exploratory_only':{'function':'pi*(1-q^8)/(4*q^2*(-log(q)))','initial_interval':[.05,.95],'golden_iterations':100,'q':q,'necessary_cost_value':cost(q),'scope':'not certified global minimum; not feasible majorant upper bound'},'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [proof,Path(__file__).resolve()]}}
    out=ROOT/'computation/results/killed-green-dihedral-check.json';out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':run()
