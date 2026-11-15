#!/usr/bin/env python3
"""Exact scalar sanity checks; analytic proof is not replaced by these tests."""
import hashlib,json,math,platform,sys
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def sqrtr(q):
    a,b=math.isqrt(q.numerator),math.isqrt(q.denominator)
    assert a*a==q.numerator and b*b==q.denominator
    return Q(a,b)

def run():
    checks={}
    alpha=Q(3,2)
    for v in (Q(1,4),Q(1,16),Q(1,64)):
        invroot=1/sqrtr(v)
        # Integral over t in [v,1] of 2(t-v)t^(-1-alpha).
        direct=2*(2*(invroot-1)-v*Q(2,3)*(invroot**3-1))
        full=2*invroot/(alpha*(alpha-1))
        tail=2/(alpha-1)-2*v/alpha
        assert direct==full-tail
        checks[f'cusp_integral_v_{v}']={'direct_to_1':str(direct),'full':str(full),'tail_after_1':str(tail),'equal':True}
    # Fourier reciprocal for p=1/2 and a=1/4, normalized by C_p.
    a=Q(1,4)
    psi_a1=(1-a*a)/(1-sqrtr(a))
    psi_a0=a*sqrtr(a)
    defect=psi_a1-psi_a0-1
    assert defect==Q(3,4)
    assert 4*defect==3
    checks['CND_four_point']={'p':'1/2','a':'1/4','psi_a1':str(psi_a1),'psi_a0':str(psi_a0),'parallelogram_defect':str(defect),'quadratic_form':str(4*defect)}
    # Optimized hinge cost at p=1/2,1: candidate p=1/2 is excluded above.
    checks['hypothetical_cost_p_half']=str(Q(4,3)**4)
    assert Q(4,3)**4==Q(256,81)<4
    checks['cauchy_cost']=4
    for t in (Q(1,4),Q(1,2),Q(3,4),Q(5,4)):
        lhs=1/t-4*(1-t)
        rhs=(2*t-1)**2/t
        assert lhs==rhs>=0
        checks[f'hinge_gap_t_{t}']=str(lhs)
    # Exact partition algebra: in the cusp branch p+alpha>=2, p<1 implies alpha>1.
    examples=[(Q(1,4),Q(1,2)),(Q(1,2),Q(3,2)),(Q(3,4),Q(7,4))]
    checks['exponent_examples']=[{'p':str(p),'alpha':str(a),'2-p-alpha':str(2-p-a),'branch':'cutoff' if p+a<2 else 'axis_cusp'} for p,a in examples]
    for p,a in examples:
        assert p+a<2 or a>1
    sources=[ROOT/'proofs/kernel-interface-audit.md',ROOT/'proofs/kernel-family-floor.md',ROOT/'proofs/kernel-subcritical-cusp-obstruction.md',Path(__file__).resolve()]
    return {'status':'PASS_scalar_checks_only','scope':'Exact rational scalar identities; no grid C4, LP, old-kernel full check, or numerical proof of the general theorem.','python':sys.version,'platform':platform.platform(),'parameters':{'alpha_for_integral':'3/2','v':['1/4','1/16','1/64'],'normalization':'half-line second difference','arithmetic':'Python Fraction','random_seed':None,'tolerance':'exact'},'checks':checks,'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sources}}

if __name__=='__main__':
    assert __import__('os').path.ismount('/Volumes/data')
    out=ROOT/'computation/results/kernel-subcritical-cusp-check.json'
    out.parent.mkdir(parents=True,exist_ok=True)
    result=run();out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
