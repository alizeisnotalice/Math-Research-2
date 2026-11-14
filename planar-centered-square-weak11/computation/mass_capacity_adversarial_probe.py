"""Bounded adversarial probe of P(mu)=mass; exact continuous capacities.

Single-thread HiGHS; process CPU search budget <=105 s (120 s total target).
Reuses frozen coefficient enumeration read-only; never edits R5/R9 files.
"""
from fractions import Fraction as F
from itertools import product
import argparse,json,random,time,platform,hashlib
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import linprog
from general_dual_mass_probe import coefficients


def disjoint_partition(labels,n):
    full=(1<<n)-1
    masks=[sum(1<<i for i in s) for s in labels]
    memo={0:[]}
    def visit(remaining):
        if remaining in memo:return memo[remaining]
        i=(remaining&-remaining).bit_length()-1
        for mask in masks:
            if mask&(1<<i) and mask&remaining==mask:
                tail=visit(remaining^mask)
                if tail is not None:
                    memo[remaining]=[mask]+tail
                    return memo[remaining]
        memo[remaining]=None
        return None
    return visit(full)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    parser.add_argument('--cpu-budget',type=float,default=105)
    parser.add_argument('--targeted',action='store_true')
    args=parser.parse_args()
    if args.targeted:
        start=time.process_time()
        points=[(-F(1,50),F(0)),(F(1,50),F(0))]+[(F(x,10),F(y,10)) for x,y in product((-1,1),repeat=2)]+[(F(1),F(1))]
        masses=[F(1,10),F(1,10)]+[F(1,1000)]*5
        c=coefficients(points,masses)
        parts=[tuple(range(6)),(6,)]
        assert c[(0,)]<masses[0] and c[(1,)]<masses[1]
        assert all(c[s]==sum(masses[i] for i in s) for s in parts)
        overlap_points=[(F(x),F(0)) for x in ('0','9/10','47/50','46/25')]
        overlap_masses=[F(1,2)]*4
        overlap_c=coefficients(overlap_points,overlap_masses)
        assert overlap_c[(0,1,2)]==overlap_c[(1,2,3)]==F(3,2)
        assert (0,1,2,3) not in overlap_c
        result={'kind':'deliberate seven-source example; no general equality claim',
                'points':points,'masses':masses,'M':sum(masses),'diameter_squared':F(121,100),
                'singleton_capacities':[c[(i,)] for i in range(7)],
                'exact_saturated_partition':parts,
                'all_capacities':{','.join(map(str,s)):v for s,v in c.items()},
                'exact_equality_certified':True,'cpu_seconds':time.process_time()-start,
                'overlapping_saturated_union_failure':{'points':overlap_points,'masses':overlap_masses,
                     'saturated_S':[0,1,2],'saturated_T':[1,2,3],'each_capacity':F(3,2),
                     'union_has_no_valid_strict_square':True,'union_diameter_squared':F(2116,625),
                     'union_mass':F(2),'all_capacities':{','.join(map(str,s)):v for s,v in overlap_c.items()}},
                'python':platform.python_version(),'coefficient_source_sha256':hashlib.sha256(Path(__file__).with_name('general_dual_mass_probe.py').read_bytes()).hexdigest()}
        Path(args.out).write_text(json.dumps(result,default=str,indent=2))
        print(json.dumps({k:v for k,v in result.items() if k!='all_capacities'},default=str,indent=2))
        return
    rng=random.Random(20261116)
    start=time.process_time(); wall=time.time()
    assert 0<args.cpu_budget<=105
    counts={'tried':0,'nontrivial_singletons':0,'multiple_unsaturated_singletons':0,'no_disjoint_saturated_partition':0,'exact_equalities':0,'counterexamples':0}
    examples=[]; partition_examples=[]; best_ratio=1.; best=None
    for trial in range(10000):
        if time.process_time()-start>args.cpu_budget:break
        n=rng.choice([5,6,7,8,9])
        # Include uneven coordinate gaps and clustered interior sites.
        if trial%3==0:
            universe=list(product([0,1,2,6,7],repeat=2))
        elif trial%3==1:
            universe=list(product(range(7),repeat=2))
        else:
            universe=[(-4,-4),(-4,4),(4,-4),(4,4),(-1,0),(0,0),(1,0),(0,-1),(0,1),(2,1),(-2,-1)]
        points=[tuple(map(F,p)) for p in rng.sample(universe,n)]
        diameter=max(max(p[k] for p in points)-min(p[k] for p in points) for k in (0,1))
        target_mass=diameter**2*F(rng.randrange(3,19),20)
        weights=[rng.choice([1,2,3,5,10,30,100,300]) for _ in points]
        masses=[target_mass*F(w,sum(weights)) for w in weights]
        c=coefficients(points,masses)
        tight=[s for s,v in c.items() if v==sum(masses[i] for i in s)]
        singleton_tight=[(i,) in tight for i in range(n)]
        counts['tried']+=1
        if all(singleton_tight):
            counts['exact_equalities']+=1
            continue
        counts['nontrivial_singletons']+=1
        multiple=sum(not v for v in singleton_tight)>=2
        if multiple:counts['multiple_unsaturated_singletons']+=1
        partition=disjoint_partition(tight,n)
        if partition is not None:
            counts['exact_equalities']+=1
            if multiple and len(partition_examples)<3:
                labels_partition=[tuple(i for i in range(n) if mask&(1<<i)) for mask in partition]
                partition_examples.append({'trial':trial,'points':points,'masses':masses,'M':target_mass,
                                           'diameter_squared':diameter**2,'singleton_saturated':singleton_tight,
                                           'partition':[{'S':s,'c_S':c[s],'mass_S':sum(masses[i] for i in s)} for s in labels_partition]})
            continue
        counts['no_disjoint_saturated_partition']+=1
        labels=list(c)
        G=np.array([[-int(i in s) for i in range(n)] for s in labels],float)
        sol=linprog(np.ones(n),A_ub=G,b_ub=[-float(c[s]) for s in labels],bounds=(0,None),method='highs',options={'threads':1})
        assert sol.success
        z=[F(-float(v)).limit_denominator(10000000) for v in sol.ineqlin.marginals]
        loads=[sum(v for s,v in zip(labels,z) if i in s) for i in range(n)]
        dual=sum(c[s]*v for s,v in zip(labels,z))
        exact_equal=all(v>=0 for v in z) and all(v<=1 for v in loads) and dual==target_mass
        record={'trial':trial,'points':points,'masses':masses,'M':target_mass,'diameter_squared':diameter**2,
                'singleton_saturated':singleton_tight,'all_capacities':{','.join(map(str,s)):v for s,v in c.items()},
                'saturated_labels':tight,'exact_equality_certified':exact_equal,
                'positive_dual':[{'S':s,'c_S':c[s],'z':v} for s,v in zip(labels,z) if v],
                'loads':loads,'dual_objective':dual,'float_objective':float(sol.fun)}
        if exact_equal:
            counts['exact_equalities']+=1
            if len(examples)<3:examples.append(record)
        else:
            w=[F(float(v)).limit_denominator(10000000) for v in sol.x]
            # Repair floating residue only upward, then verify all exact capacities.
            deficit=max([F(0)]+[c[s]-sum(w[i] for i in s) for s in labels])
            w=[v+deficit for v in w]
            ratio=float(sum(w)/target_mass)
            if ratio<best_ratio:
                best_ratio=ratio; best=record
                best['candidate_reduced_weights']=w
                best['reduced_total']=sum(w)
            if all(v>=0 for v in w) and all(sum(w[i] for i in s)>=c[s] for s in labels) and sum(w)<target_mass:
                counts['counterexamples']+=1
                best['exact_counterexample']=True
                break
    result={'seed':20261116,'counts':counts,'search_cpu_seconds':time.process_time()-start,'wall_seconds':time.time()-wall,
            'budget_cpu_seconds':args.cpu_budget,'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
            'coefficient_source_sha256':hashlib.sha256(Path(__file__).with_name('general_dual_mass_probe.py').read_bytes()).hexdigest(),
            'capacities':'all exact incidences, complete continuous-radius layers and endpoints',
            'best_ratio':best_ratio,'best':best,'nonpartition_equality_examples':examples,
            'multiple_exception_partition_examples':partition_examples}
    Path(args.out).write_text(json.dumps(result,default=str,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ('best','nonpartition_equality_examples','multiple_exception_partition_examples')},indent=2))


if __name__=='__main__':main()
