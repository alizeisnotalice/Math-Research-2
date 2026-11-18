import os
"""Protected-peak upper bounds for standard and T-shaped rectangular tilings."""
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import permutations, product
import json
from math import comb, factorial, gcd, prod
from pathlib import Path
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
PROJECT=Path(os.environ.get('A022_RESEARCH_ROOT', Path(__file__).resolve().parents[3]))
R49=PROJECT/'certificates/r049/odd-quadrants/four-block-bound.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verified_primes(path):
    certificate = json.loads(path.read_text())
    known = set()

    def visit(p):
        if p in known:
            return
        node = certificate['nodes'][str(p)]
        assert int(node['p']) == p
        if p == 2:
            assert node == {'p': '2', 'base': 'two'}
        else:
            factors = {int(q): k for q, k in node['p_minus_1_factorization'].items()}
            assert all(2 <= q < p and k >= 1 for q, k in factors.items())
            assert prod(q**k for q, k in factors.items()) == p - 1
            for q in factors:
                visit(q)
            a = int(node['witness'])
            assert 1 < a < p and pow(a, p - 1, p) == 1
            assert all(gcd(pow(a, (p - 1) // q, p) - 1, p) == 1 for q in factors)
        known.add(p)

    primes = tuple(map(int, certificate['roots']))
    assert len(primes) == len(set(primes)) == 5
    for p in primes:
        assert p.bit_length() == 61
        visit(p)
    assert prod(primes) > factorial(61)
    return primes, len(known)


def independent_sizes(profile):
    """Forward cardinality DP, independent of the C++ suffix DP."""
    last = {0: [1]}
    for width in map(int, profile):
        after = {}
        for s in range(1 << width):
            if s & (s << 1):
                continue
            shift = bin(s).count('1')
            current = [0] * (max(map(len, last.values())) + shift)
            for t, counts in last.items():
                if s & t:
                    continue
                for k, count in enumerate(counts):
                    current[k + shift] += count
            after[s] = current
        last = after
    histogram = [0] * max(map(len, last.values()))
    for counts in last.values():
        for k, count in enumerate(counts):
            histogram[k] += count
    while not histogram[-1]:
        histogram.pop()
    return histogram


def crt(residues, primes):
    result, modulus = 0, 1
    for r, p in zip(residues, primes):
        result += ((r - result) * pow(modulus, -1, p) % p) * modulus
        modulus *= p
    assert all(result % p == r for r, p in zip(residues, primes))
    return result


def dec(value):
    r=Fraction(value)
    with localcontext() as ctx:
        ctx.prec=57
        return str(Decimal(r.numerator)/Decimal(r.denominator))


def plan(shape):
    w,h,edges=shape
    argv=[str(HERE/'rect_protected'),str(w),str(h),edges,'0','--plan','--max-table-mib','128']
    result=subprocess.run(argv,capture_output=True,text=True,timeout=10,check=True)
    assert result.stderr==''
    data=json.loads(result.stdout)
    a,b=max(0,w-int('L' in edges)-int('R' in edges)),max(0,h-int('T' in edges)-int('B' in edges))
    x,y=sorted((a,b))
    expected=independent_sizes(str(x)*y)
    assert data['protected_vertices']==a*b and data['area']==w*h
    assert data['states_by_layer']==data['values_by_layer']==expected
    assert data['max_table_bytes']==max(8*(s+t) for s,t in zip([0]+expected,expected))
    return data


def worker(shape,idx,prime,predicted,out):
    w,h,em=shape
    tag=f'{w}x{h}-{em}-p{idx}'
    cmd=['/usr/bin/time','-l',str(HERE/'rect_protected'),str(w),str(h),em,str(idx),'--max-table-mib','128']
    result=subprocess.run(cmd,capture_output=True,text=True,timeout=50,check=True)
    q=json.loads(result.stdout)
    assert (q['width'],q['height'],q['unsafe_edges'])==shape and q['mode']=='count'
    assert q['processed_by_layer']==q['states_by_layer']==predicted['states_by_layer']
    assert int(q['prime'])==prime and q['prime_index']==idx
    assert q['max_table_bytes']==predicted['max_table_bytes']
    m=re.search(r'^\s*(\d+)\s+maximum resident set size$',result.stderr,re.M)
    assert m is not None
    p=out/f'{tag}.json';qf=out/f'{tag}.resource'
    p.write_text(result.stdout);qf.write_text(result.stderr)
    receipt={'shape':list(shape),'prime_index':idx,'stdout_sha256':digest(p),'resource_sha256':digest(qf),
             'max_rss_bytes':int(m[1]),'exit_status':result.returncode}
    return (shape,idx),q,receipt


def readback(shape,raw,primes):
    keys=('peak_frequency_residues','unprotected_max_frequency_residues','bad_second_frequency_residues')
    w,h,e=shape
    data=[raw[shape,i] for i in range(len(primes))]
    f,t,v=[[crt(series,primes) for series in zip(*(r[k] for r in data))] for k in keys]
    assert sum(f)==factorial(w*h)
    assert sum(t)==(w*h-data[0]['protected_vertices'])*factorial(w*h-1)
    assert all(0<=z<=x-y<=x for x,y,z in zip(f,t,v))
    return [f,t,v]


def merge(areas,histograms):
    B=range(len(areas))
    @lru_cache(maxsize=None)
    def interleaving(*prefix):
        occupancy=list(areas)
        for i in prefix: occupancy[i]-=1
        if min(occupancy)<0: return 0
        return factorial(sum(occupancy))//prod(factorial(a) for a in occupancy)
    raw=first_bad=second_bad=third_bad=0
    for orders in product(range(3),repeat=len(areas)):
        if sum(orders)>2: continue
        f,t,v=[[histograms[i][j][orders[i]] if orders[i]<len(histograms[i][j]) else 0 for i in B] for j in range(3)]
        raw+=interleaving()*prod(f)
        if sum(orders)!=2:continue
        s=[f[i]-t[i] for i in B]
        for i in B:
            separate=prod(f[z] for z in B if z!=i)
            first_bad+=interleaving(i)*t[i]*separate
            second_bad+=interleaving(i,i)*v[i]*separate
            for j in B:
                if j==i:continue
                separate2=prod(f[z] for z in B if z!=i and z!=j)
                second_bad+=interleaving(i,j)*s[i]*t[j]*separate2
                third_bad+=(interleaving(i,j,i)*v[i]*s[j]+interleaving(i,j,j)*s[i]*v[j])*separate2
                for k in B:
                    if k in (i,j):continue
                    third_bad+=interleaving(i,j,k)*s[i]*s[j]*t[k]*prod(f[z] for z in B if z not in (i,j,k))
    up2=raw-first_bad-second_bad
    up3=up2-third_bad
    assert 0<up3<=up2<=raw<=factorial(sum(areas))
    return [raw,first_bad,second_bad,third_bad,up2,up3]


def brute(shape):
    w,h,e=shape
    P={y*w+x for y in range(int('T' in e),h-int('B' in e)) for x in range(int('L' in e),w-int('R' in e))}
    nbr=[{v for v in range(w*h) if abs(u%w-v%w)+abs(u//w-v//w)==1} for u in range(w*h)]
    f,t,v=([0]*(len(P)+1) for _ in range(3))
    for word in permutations(range(w*h)):
        idx={x:i for i,x in enumerate(word)}
        k=sum(all(idx[u]<idx[z] for z in nbr[u]) for u in P)
        f[k]+=1
        t[k]+=word[0] not in P
        v[k]+=word[0] in P and word[1] not in P and word[1] not in nbr[word[0]]
    return f,t,v


def brute_merge_gate():
    # Three abstract path blocks, each with vertex zero protected.  The
    # paths need not form a grid tiling for the shuffle count to be tested.
    sizes=[3,3,2]
    charts=[]
    local=[0,1,2,0,1,2,0,1]
    region=[0,0,0,1,1,1,2,2]
    # Count the histograms independently: protected vertex 0 is a peak iff
    # it outranks 1. For a three-path, local maximum 2 allows one such order.
    for a in sizes:
        f,t,v=([0,0] for _ in range(3))
        for perm in permutations(range(a)):
            z=int(perm.index(0)<perm.index(1))
            f[z]+=1;t[z]+=perm[0]!=0;v[z]+=perm[0]==0 and perm[1]>1
        charts.append([f,t,v])
    direct=[0]*6
    for word in permutations(range(sum(sizes))):
        seen={v:i for i,v in enumerate(word)}
        Z=sum(seen[r]<seen[r+1] for r in (0,3,6))
        if Z>2:continue
        direct[0]+=1
        a,b,c=word[:3]
        sa,sb,sc=(local[u]==0 for u in (a,b,c))
        def adjacent(x,y):return region[x]==region[y] and abs(local[x]-local[y])==1
        if Z==2 and not sa:
            direct[1]+=1
        elif Z==2 and not sb and not adjacent(a,b):
            direct[2]+=1
        else:
            direct[4]+=1
            if Z==2 and sa and sb and region[a]!=region[b] and not sc and not adjacent(a,c) and not adjacent(b,c):
                direct[3]+=1
            else:
                direct[5]+=1
    assert direct==merge(sizes,charts) and direct[3]>0
    return direct


def main():
    ap=ArgumentParser()
    ap.add_argument('out',type=Path)
    ap.add_argument('--workers',type=int,choices=(1,2,3),default=2)
    args=ap.parse_args()
    args.out.mkdir(exist_ok=False)
    tilings={'quadrants':[(6,6,'BR')]*4,'three_block_T':[(12,4,'B'),(6,8,'TR'),(6,8,'TR')]}
    gate=[(2,3,'BT'),(3,3,'TR')]
    old=[(5,5,'BR'),(5,6,'BR'),(6,5,'BR'),(6,6,'BR')]
    keys=sorted(set(gate+old+[s for layout in tilings.values() for s in layout]))
    plan_all={s:plan(s) for s in keys}
    # Plans only. These comparator layouts are never sent to the worker pool.
    deferred=[plan(s) for s in ((12,6,'B'),(12,5,'B'),(6,7,'TR'))]
    assert all(p['fits_budget'] for p in plan_all.values())
    all_primes,nodes=verified_primes(HERE/'pratt-61bit.json')
    primes=all_primes[:4]
    assert prod(primes)>factorial(max(w*h for w,h,_ in keys))
    jobs=[(s,i,p,plan_all[s],args.out) for s in keys for i,p in enumerate(primes)]
    with ThreadPoolExecutor(max_workers=args.workers) as group:
        observed=list(group.map(lambda z:worker(*z),jobs))
    residue={k:v for k,v,_ in observed}
    hist={k:readback(k,residue,primes) for k in keys}
    for s in gate:
        for x,y in zip(brute(s),hist[s]):assert x==y+[0]*(len(x)-len(y)),s
    previous=json.loads(R49.read_text())
    for actual,expect in zip((hist[s] for s in old),previous['sides']):
        for a,b in zip(actual,expect):assert a+[0]*(len(b)-len(a))==b
    old_refine=merge([w*h for w,h,_ in old],[hist[s] for s in old])
    assert old_refine[4]==int(previous['both_protected_peaks_upper'])
    algebra_gate=brute_merge_gate()
    results={}
    for name,layout in tilings.items():
        sizes=[w*h for w,h,_ in layout]
        row=merge(sizes,[hist[s] for s in layout]);F=factorial(sum(sizes))
        result={'shape_representatives':[list(s) for s in layout],'protected_count':sum(plan_all[s]['protected_vertices'] for s in layout),
                'protected_count_le_two':str(row[0]),'excluded_global_max':str(row[1]),'excluded_second_global':str(row[2]),
                'excluded_third_global':str(row[3]),'top_two_upper':str(row[4]),'top_three_upper':str(row[5]),
                'top_two_probability':dec(Fraction(row[4],F)),'top_three_probability':dec(Fraction(row[5],F))}
        results[name]=result
        print(name,result['top_two_probability'],result['top_three_probability'])
    ledger={'scope':'Exact finite upper counts from necessary protected-peak conditions','python':sys.version,'pratt_nodes_verified':nodes,
            'computed_primes':list(map(str,primes)),'counter_sha256':digest(HERE/'rect_protected.cpp'),'driver_sha256':digest(Path(__file__)),
            'computed_results':results,'deferred_plan_only_regions':deferred,'n11_legacy_ledger_sha256':digest(R49),
            'n11_quadrants_top_two_checked':str(old_refine[4]),'n11_quadrants_extra_third_filter':str(old_refine[5]),
            'symbolic_three_block_full_permutation_check':algebra_gate,'waited_jobs':len(observed),
            'max_worker_rss':max(x['max_rss_bytes'] for _,_,x in observed),'receipts':[x for _,_,x in observed],
            'regional_histograms':{f'{w}x{h}-{e}':[list(map(str,a)) for a in hist[w,h,e]] for w,h,e in keys},
            'regional_plans':{f'{w}x{h}-{e}':plan_all[w,h,e] for w,h,e in keys}}
    (args.out/'summary.json').write_text(json.dumps(ledger,indent=2)+'\n')


if __name__=='__main__':main()
