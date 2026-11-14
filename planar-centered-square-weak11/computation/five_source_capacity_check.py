"""Small exact checks of the five-source structural proof; no random search/LP."""
from fractions import Fraction as F
from itertools import combinations,permutations,product
from pathlib import Path
import argparse,hashlib,json,platform,time
from general_dual_mass_probe import coefficients


def full_quadrants(p,others):
    return all(any(s*(q[0]-p[0])>=0 and t*(q[1]-p[1])>=0 for q in others)
               for s,t in product((-1,1),repeat=2))


def diagonal_pattern(guards):
    for s,r,t in permutations(guards):
        if s[0]<0 and s[1]==0 and r[0]==1 and r[1]>1 and t[0]>=1 and t[1]<=0:
            return 'staircase'
        if s[0]==0 and s[1]<0 and r[0]>1 and r[1]==1 and t[0]<=0 and t[1]>=1:
            return 'swapped_staircase'
    return None


def aligned_pattern(guards):
    left=[p for p in guards if p[0]<0 and p[1]==0]
    right=[p for p in guards if p[0]>1 and p[1]==0]
    if left and right:return 'four_collinear'
    if left:
        assert len(left)==1
        others=[p for p in guards if p!=left[0]]
        if all(p[0]>=1 for p in others) and others[0][1]*others[1][1]<0:
            return 'left_anchor_split'
    if right:
        assert len(right)==1
        others=[p for p in guards if p!=right[0]]
        if all(p[0]<=0 for p in others) and others[0][1]*others[1][1]<0:
            return 'right_anchor_split'
    return None


def fixture(points,masses,partition):
    points=[tuple(map(F,p)) for p in points]; masses=list(map(F,masses))
    c=coefficients(points,masses)
    assert c[(0,)]<masses[0] and c[(1,)]<masses[1]
    assert sorted(i for s in partition for i in s)==list(range(5))
    assert all(c[tuple(s)]==sum(masses[i] for i in s) for s in partition)
    D=max(max(p[k] for p in points)-min(p[k] for p in points) for k in (0,1))
    return {'points':points,'masses':masses,'M':sum(masses),'diameter_squared':D*D,
            'all_capacities':{','.join(map(str,s)):v for s,v in c.items()},
            'singleton_capacities':[c[(i,)] for i in range(5)],'saturated_partition':partition,
            'P_equals_M_verified':True}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);args=parser.parse_args()
    start=time.process_time()
    counts={}
    for name,q,check in [('diagonal',(1,1),diagonal_pattern),('aligned',(1,0),aligned_pattern)]:
        p=(0,0); universe=[v for v in product(range(-2,4),repeat=2) if v not in (p,q)]
        examined=0; matches={}
        for guards in combinations(universe,3):
            examined+=1
            if full_quadrants(p,[q,*guards]) and full_quadrants(q,[p,*guards]):
                pattern=check(guards)
                assert pattern is not None,(name,guards)
                matches[pattern]=matches.get(pattern,0)+1
        counts[name]={'guard_triples_examined':examined,'full_quadrant_cases':matches}
    examples={
        'axis_aligned_two_exceptions':fixture([(0,0),(1,0),(-1,0),(4,1),(4,-1)],
                                           [5,17,'1/10','1/10','1/10'],[(0,2),(1,3,4)]),
        'diagonal_two_exceptions':fixture([(0,0),(1,1),(-1,0),(1,2),(2,-1)],
                                         [10,10,'1/10','1/10','1/10'],[(0,1,2,3,4)])}
    a=examples['axis_aligned_two_exceptions']; d=examples['diagonal_two_exceptions']
    assert a['M']<a['diameter_squared'] and d['M']>d['diameter_squared']
    triangle_points=[(F(0),F(0)),(F(2),F(0)),(F(1),F(2))]
    triangle_masses=[F(3)]*3
    triangle_c=coefficients(triangle_points,triangle_masses)
    assert all(triangle_c[s]==6 for s in [(0,1),(0,2),(1,2)])
    unbalanced={'points':triangle_points,'masses':triangle_masses,
                'odd_cycle_labels':[(0,1),(0,2),(1,2)],
                'capacities':[triangle_c[s] for s in [(0,1),(0,2),(1,2)]],
                'meaning':'Unbalanced incidence submatrix; not a failure of fractional perfect cover.'}
    result={'scope':'finite sanity checks of a separate analytic five-source proof; no sampled proof for arbitrary coordinates',
            'python':platform.python_version(),'arithmetic':'Fraction and exact integer order tests',
            'random_search':False,'LP_solver_used':False,'grid_counts':counts,'exact_examples':examples,
            'unbalanced_hypergraph_example':unbalanced,
            'cpu_seconds':time.process_time()-start,
            'capacity_enumerator_sha256':hashlib.sha256(Path(__file__).with_name('general_dual_mass_probe.py').read_bytes()).hexdigest()}
    Path(args.out).write_text(json.dumps(result,default=str,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='exact_examples'},default=str,indent=2))


if __name__=='__main__':main()
