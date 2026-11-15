#!/usr/bin/env python3
"""Small scalar checks; global theorem is analytic, optimizer is exploratory only."""
import sys,math,json,hashlib,platform
from fractions import Fraction as Q
from pathlib import Path
import flint
from flint import arb,ctx
ctx.prec=128
ROOT=Path(__file__).resolve().parents[1]
N_H=128
N_S=4096

def H(q):
    b=1.;power=q;out=0.
    for n in range(N_H):
        out+=b*power/(4*n+1)
        b*=(2*n+1)/(2*n+2);power*=q**4
    return out

def denominator(q):return H(q)**2*(-math.log(q))

def main():
    assert __import__('os').path.ismount('/Volumes/data')
    c=arb.pi()*(1+arb(2).sqrt())/2
    assert c>arb(379)/100
    rational=Q(31415,10000)*(1+Q(14142,10000))/2
    assert rational>Q(379,100)
    assert arb.pi()>arb(31415)/10000
    assert arb(2).sqrt()>arb(14142)/10000
    b=1.;S=0.
    for n in range(N_S+1):
        S+=b*b/(4*n+1)**2;b*=(2*n+1)/(2*n+2)
    lo,hi=.05,.95;g=(math.sqrt(5)-1)/2
    for _ in range(100):
        a=hi-g*(hi-lo);b=lo+g*(hi-lo)
        if denominator(a)<denominator(b):lo=a
        else:hi=b
    q=(lo+hi)/2
    I=math.gamma(.25)*math.sqrt(math.pi)/(4*math.gamma(.75))
    proof=ROOT/'proofs/killed-green-homothetic-diamond-obstruction.md'
    result={'status':'PASS_scalar_interval_and_rational_checks','python':sys.version,'flint':flint.__version__,'platform':platform.platform(),'arithmetic':{'arb_precision_bits':128,'final_decimal_claim':'strict lower 379/100','random_seed':None,'workers':1},'rigorous_scalar_checks':{'cost_constant_arb':str(c),'constant_gt_379_over_100':True,'rational_lower':str(rational),'pi_gt_31415_over_10000':True,'sqrt2_gt_14142_over_10000':True},'exploratory_only':{'H_terms':N_H,'S_inclusive_index':N_S,'S_partial_float':S,'S_tail_bound_formula':'1/(32*N_S^2); floating rounding not enclosed here','golden_section_iterations':100,'search_interval':[.05,.95],'q':q,'R':I/H(q),'vertex_necessary_cost':math.pi*S/(4*denominator(q)),'meaning':'Numerical minimum of vertex-based necessary cost; neither certified global minimum nor feasible upper bound.'},'sources_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [proof,Path(__file__).resolve()]}}
    out=ROOT/'computation/results/killed-green-diamond-check.json';out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
