"""Recover the two-machine n9 integer ledger and disjoint class weights."""
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys

base=Path(__file__).resolve().parent
src=base.parents[1]/'r046/layered'
sys.path.insert(0,str(src))
from run_crt import independent_sizes,safe_roots,verified_primes

n=9
primes,_=verified_primes(src/'pratt-61bit.json')
P=math.prod(primes)
factors=[(P//p)*pow(P//p,-1,p) for p in primes]
plan=json.loads((src/'reference-plan.json').read_text())
pairs={d['profile']:tuple(s['profile'] for s in d['sides']) for d in plan['profile_pairs']}
batches=json.loads((base/'residue-ledger.json').read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def complement(h):return ''.join(str(n-int(c)) for c in h)
def mirror_root(v):return (n-1-v//n)*n+v%n
def transpose(x):return x%n*n+x//n
def reflect_x(x):return x//n*n+n-1-x%n

vectors={}
receipts=[]
for batch,expected_jobs in zip(batches.values(),(6,42,42),strict=True):
    summary=batch['summary']
    assert summary['counter_source_sha256']==sha(src/'root_mod.cpp')
    records=batch['receipts']
    assert len(records)==expected_jobs and summary['completed_jobs']==expected_jobs
    residues={}
    for rec in records:
        h,i=rec['profile'],rec['prime_index']
        data=batch['stdout'][f'{h}-p{i}']
        assert hashlib.sha256(data.encode()).hexdigest()==rec['stdout_sha256']
        x=json.loads(data)
        assert x['processed_by_layer']==independent_sizes(h)
        assert int(x['prime'])==primes[i] and rec['exit_code']==0
        assert (h,i) not in residues
        residues[h,i]=list(map(int,x['root_residues']))
        receipts.append(dict(profile=h,prime=i,rss=rec['max_rss_bytes'],sha256=rec['stdout_sha256']))
    for h,counts in summary['root_counts'].items():
        check=[sum(a*b for a,b in zip(column,factors))%P
               for column in zip(*(residues[h,i] for i in range(3)))]
        assert check==list(map(int,counts)) and h not in vectors
        assert sum(check)<=math.factorial(sum(map(int,h)))
        assert {i for i,x in enumerate(check) if x}==safe_roots(h)
        vectors[h]={i:x for i,x in enumerate(check) if x}

profiles=sorted({g for h in pairs for g in (h,h[::-1])})
assert len(profiles)==28 and len(vectors)==30
assert len({tuple(int(c)%2 for c in h) for h in profiles})==1
for h,k in itertools.combinations(profiles,2):
    assert any(abs(int(x)-int(y))>=2 for x,y in zip(h,k))
    assert all((int(x)-int(y))%2==0 for x,y in zip(h,k))
def vector(h):
    if h in vectors:return vectors[h]
    return {mirror_root(v):x for v,x in vectors[h[::-1]].items()}
counts=Counter()
event_counts={}
near_count=0
for h in profiles:
    factor=math.comb(n*n,sum(map(int,h)))
    total=0
    for u,a in vector(h).items():
        for w,b in vector(complement(h)).items():
            v=reflect_x(w)
            d=abs(u%n-v%n)+abs(u//n-v//n)
            assert d>n/3
            value=factor*a*b
            counts[tuple(sorted((u,v)))]+=value
            total+=value
            if d<=n/2:near_count+=value
    event_counts[h]=total
    assert total>0
total=sum(event_counts.values())
assert total==sum(counts.values())
classes=[]
for f in (lambda x:x,reflect_x,transpose,lambda x:transpose(reflect_x(x))):
    other=Counter({tuple(sorted(f(v) for v in r)):x for r,x in counts.items()})
    assert sum(other.values())==total and len(other)==len(counts)
    classes.append(other)
maxima={}
for other in classes:
    for r,a in other.items():
        maxima[r]=max(a,maxima.get(r,0))
mass=sum(maxima.values())
assert total<=mass<=4*total
summary=dict(profiles=profiles,independent_regions=len(vectors),residue_jobs=len(receipts),
             cutoff='rho>n/3',one_class_count=total,one_class_count_distance_le_n_over_2=near_count,
             rootwise_subunion_count=mass,root_pairs=len(maxima),
             one_class_unconditional=str(Fraction(total,math.factorial(81))),
             rootwise_unconditional=str(Fraction(mass,math.factorial(81))),
             event_counts=event_counts,rootwise_counts=[dict(peaks=r,count=v) for r,v in sorted(maxima.items())],
             residue_input_receipts=receipts)
with Path(sys.argv[1]).open('x') as f:
    json.dump(summary,f,indent=2)
    f.write('\n')
print(total,mass,len(maxima))
