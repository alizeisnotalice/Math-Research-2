"""Exact continuous-radius incidence coefficients, floating LP route probe."""
from fractions import Fraction as F
from itertools import product
import json, random, platform, sys
from datetime import datetime, timezone


def masks(points,r):
    axes=[]
    for k in (0,1):
        ends=sorted({p[k]+s*r for p in points for s in (-1,1)})
        axes.append(ends+[(a+b)/2 for a,b in zip(ends,ends[1:])])
    return {tuple(i for i,p in enumerate(points) if max(abs(x-p[0]),abs(y-p[1]))<=r)
            for x,y in product(*axes)}-{()}


def coefficients(points,masses):
    critical=sorted({F(0)} | {abs(p[k]-q[k])/2 for p,q in product(points,repeat=2) for k in (0,1)})
    # Last finite cap exceeds sqrt(total mass)/2 and all arrangement changes.
    critical.append(max(critical[-1]+1,sum(masses)+1))
    coeff={}
    for a,b in zip(critical,critical[1:]):
        for s in masks(points,(a+b)/2):
            mass=sum(masses[i] for i in s)
            if 4*a*a<mass:
                coeff[s]=max(coeff.get(s,F(0)),min(mass,4*b*b))
    for r in critical[1:-1]:
        for s in masks(points,r):
            if 4*r*r<sum(masses[i] for i in s):
                coeff[s]=max(coeff.get(s,F(0)),4*r*r)
    return coeff


def main():
    import scipy
    from scipy.optimize import linprog
    rng=random.Random(20261115)
    count=0
    equality_certified=0
    found=None
    certificates=[]
    for trial in range(100):
        points=rng.sample(list(product(range(5),repeat=2)),5)
        points=[tuple(map(F,p)) for p in points]
        masses=[F(rng.randrange(1,40),4) for _ in points]
        c=coefficients(points,masses)
        rows=[[-int(i in s) for i in range(len(points))] for s in c]
        res=linprog([1.]*len(points),A_ub=rows,b_ub=[-float(v) for v in c.values()],bounds=(0,None),method='highs')
        assert res.success
        count+=1
        z=[F(-float(v)).limit_denominator(1000000) for v in res.ineqlin.marginals]
        loads=[sum(v for s,v in zip(c,z) if i in s) for i in range(len(points))]
        objective=sum(v*cost for v,cost in zip(z,c.values()))
        if (all(v>=0 for v in z) and all(v<=1 for v in loads)
            and objective==sum(masses)):
            equality_certified+=1
            certificates.append({'trial':trial,'points':points,'masses':masses,
                                 'total_mass':sum(masses),
                                 'positive_dual_support':[{'S':list(s),'c_S':c[s],'z':v}
                                                          for s,v in zip(c,z) if v>0],
                                 'source_loads':loads,'dual_objective':objective})
        if res.fun<float(sum(masses))-1e-7:
            w=[F(float(v)).limit_denominator(1000000) for v in res.x]
            if all(sum(w[i] for i in s)>=v for s,v in c.items()) and sum(w)<sum(masses):
                found={'trial':trial,'points':points,'masses':masses,'total_mass':sum(masses),
                       'feasible_reduced_weights':w,'reduced_total':sum(w),
                       'coefficients':{','.join(map(str,s)):v for s,v in c.items()}}
                break
    print(json.dumps({'python':platform.python_version(),'scipy':scipy.__version__,'seed':20261115,
                      'run_completed_utc':datetime.now(timezone.utc).isoformat(),
                      'run_purpose':'Actual rerun to archive per-input exact certificates; supersedes summary-only run',
                      'trials_run':count,'coefficient_arithmetic':'exact Fraction',
                      'P_equals_mass_exact_dual_certificates':equality_certified,
                      'lp_search':'HiGHS floating; found feasible weights rationally checked',
                      'coefficient_reproduction':'coefficients(points,masses): all exact continuous-radius arrangement layers and critical radii',
                      'strict_witness_note':'P=M is the supremum-capacity LP equality; actual strict receiver certificates only approach it as epsilon tends to zero',
                      'certificates':certificates,'found':found},default=str,indent=2))


def verify_saved(path):
    """Read-only verification using standard-library exact arithmetic; no LP solver."""
    with open(path) as handle:
        data=json.load(handle)
    certs=data['certificates']
    assert len(certs)==data['P_equals_mass_exact_dual_certificates']==100
    assert [row['trial'] for row in certs]==list(range(100))
    for cert in certs:
        points=[tuple(map(F,p)) for p in cert['points']]
        masses=list(map(F,cert['masses']))
        assert len(set(points))==len(points)==len(masses)
        assert all(m>0 for m in masses)
        c=coefficients(points,masses)
        active=cert['positive_dual_support']
        assert len({tuple(row['S']) for row in active})==len(active)
        for row in active:
            assert F(row['z'])>0 and c[tuple(row['S'])]==F(row['c_S'])
        loads=[sum(F(row['z']) for row in active if i in row['S']) for i in range(len(points))]
        value=sum(F(row['z'])*F(row['c_S']) for row in active)
        assert loads==list(map(F,cert['source_loads'])) and all(v<=1 for v in loads)
        assert value==F(cert['dual_objective'])==F(cert['total_mass'])==sum(masses)
        assert all(v<=sum(masses[i] for i in s) for s,v in c.items())
    print(json.dumps({'certificates_verified':len(certs),'arithmetic':'Fraction',
                      'LP_solver_used':False,'capacity_coefficients_recomputed':True,
                      'verified_utc':datetime.now(timezone.utc).isoformat()},indent=2))


if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='--verify':
        verify_saved(sys.argv[2])
    else:
        main()
