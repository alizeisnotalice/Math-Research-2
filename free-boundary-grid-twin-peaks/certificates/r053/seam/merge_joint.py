"""Reconstruct joint peak-mask histograms and contract T5 upper conditions."""
from argparse import ArgumentParser
from collections import Counter
from decimal import Decimal,localcontext
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import permutations,product
import json
from math import comb,factorial,gcd,prod
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent
PROJ=Path('/Volumes/data/CodexProjects/T2/free-boundary-grid-twin-peaks')
R52=PROJ/'certificates/r052/t5'

# Integer certificate and row-transfer helpers from the round-51 runner.
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



def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def dec(x):
    with localcontext() as context:
        context.prec=57
        r=Fraction(x)
        return str(Decimal(r.numerator)/Decimal(r.denominator))


def modular_inversion(raw):
    p=int(raw['prime'])
    masks=raw['mask_count']
    rows=len(raw['states_by_layer'])
    A=[[list(r) for r in raw['moments'][k*masks:(k+1)*masks]] for k in range(rows)]
    assert len(raw['moments'])==rows*masks and all(len(row)==4 for row in raw['moments'])
    assert all(0<=x<p for row in raw['moments'] for x in row)
    for k in range(rows):
        for j in range(masks.bit_length()-1):
            b=1<<j
            for m in range(masks):
                if not m&b:
                    A[k][m]=[(a-z)%p for a,z in zip(A[k][m],A[k][m|b],strict=True)]
    return [[[sum((-1)**(k-i)*comb(k,i)*A[k][m][q] for k in range(i,rows))%p
              for i in range(rows)] for q in range(4)] for m in range(masks)]


def exact_histograms(raws,primes):
    matrices=[modular_inversion(raw) for raw in raws]
    ans={}
    masks=raws[0]['mask_count']
    for m in range(masks):
        f,t,v,a=[[crt(z,primes) for z in zip(*(X[m][q] for X in matrices),strict=True)] for q in range(4)]
        h=[x-y-z-w for x,y,z,w in zip(f,t,v,a,strict=True)]
        assert all(x>=0 for x in h)
        if m&(m<<1):
            assert sum(f)==0
        ans[m]=[f,t,v,h]
    assert sum(sum(z[0]) for z in ans.values())==factorial(raws[0]['area'])
    return ans


def histogram_local(words,protected,marked,adj):
    chart={m:[[0]*(len(protected)+1) for _ in range(4)] for m in range(1<<len(marked))}
    for word in words:
        loc={v:i for i,v in enumerate(word)}
        peaks={v for v in protected|set(marked) if all(loc[v]<loc[u] for u in adj[v])}
        k=len(protected&peaks)
        mask=sum(1<<j for j,v in enumerate(marked) if v in peaks)
        f,t,v,h=chart[mask]
        f[k]+=1
        if word[0] not in protected:t[k]+=1
        elif word[1] in protected:h[k]+=1
        elif word[1] not in adj[word[0]]:v[k]+=1
    return chart


def rectangle_gate(w,h):
    cmd=[str(HERE/'joint_seam'),str(w),str(h),'0','--max-table-mib','1','--verify-ranks']
    p=subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=10)
    assert p.stderr==''
    raw=json.loads(p.stdout)
    actual=exact_histograms([raw],[int(raw['prime'])])
    H=tuple(range(w*h))
    P={y*w+x for y in range(1,h) for x in range(w-1)}
    J=[y*w+w-1 for y in range(1,h)]
    adj={v:{u for u in H if abs(v%w-u%w)+abs(v//w-u//w)==1} for v in H}
    direct=histogram_local(permutations(H),P,J,adj)
    assert all(actual[m][q]+[0]*max(0,len(direct[m][q])-len(actual[m][q]))==
               direct[m][q]+[0]*max(0,len(actual[m][q])-len(direct[m][q]))
               for m in actual for q in range(4))
    return {'rectangle':[w,h],'exhaustive_words':factorial(w*h),'mask_count':len(actual),
            'rank_and_histograms_match':True}


def get(seq,i):return seq[i] if i<len(seq) else 0

def old_chart(three):
    f,t,v=[list(map(int,z)) for z in three]
    H=[max(0,a-b-c-(t[i-1] if i else 0)) for i,(a,b,c) in enumerate(zip(f,t,v,strict=True))]
    return [f,t,v,H]


def compact_sum(hist):
    values=list(hist.values())
    return [[sum(v[k][i] for v in values) for i in range(len(values[0][k]))] for k in range(4)]


@lru_cache(maxsize=None)
def suffix(areas,prefix):
    used=Counter(prefix)
    left=[a-used[i] for i,a in enumerate(areas)]
    if min(left)<0:return 0
    return factorial(sum(left))//prod(factorial(x) for x in left)


def count_cells(areas,charts,overlap):
    B=range(len(areas))
    stat=[0]*7 # raw; rank1; rank2; rank3 different/same; rank4; rank5
    for order in product(range(3),repeat=len(areas)):
        Z=sum(order)
        if Z+overlap>2:continue
        f,t,v,h=[[get(charts[i][q],order[i]) for i in B] for q in range(4)]
        if not all(f):continue
        stat[0]+=suffix(areas,())*prod(f)
        if Z<2:continue
        s=[f[i]-t[i] for i in B]
        def c(prefix,override):
            return suffix(areas,tuple(prefix))*prod(override.get(i,f[i]) for i in B)
        for i in B:
            stat[1]+=c((i,),{i:t[i]})
            stat[2]+=c((i,i),{i:v[i]})
            for j in B:
                if j==i:continue
                stat[2]+=c((i,j),{i:s[i],j:t[j]})
                stat[3]+=c((i,j,i),{i:v[i],j:s[j]})+c((i,j,j),{i:s[i],j:v[j]})
                stat[4]+=c((i,i,j),{i:h[i],j:t[j]})
                for k in B:
                    if k in (i,j):continue
                    stat[3]+=c((i,j,k),{i:s[i],j:s[j],k:t[k]})
                if order[i]!=1 or order[j]!=1:continue
                # First 3 protected labels occupy i twice and j once.
                stat[5]+=3*c((i,i,j,j),{i:h[i],j:v[j]})
                for k in B:
                    if k in (i,j):continue
                    stat[5]+=3*c((i,i,j,k),{i:h[i],j:s[j],k:t[k]})
                    if i<j:
                        stat[6]+=6*c((i,i,j,j,k),{i:h[i],j:h[j],k:t[k]})
    assert stat[0]>=sum(stat[1:]),stat
    return stat


def count_all(areas,top,L,R,joint):
    stat=[0]*7
    checked=0
    for x,X in L.items():
        if not sum(X[0]):continue
        for y,Y in R.items():
            if not sum(Y[0]):continue
            z=(x&y).bit_count() if joint else 0
            add=count_cells(areas,[top,X,Y],z)
            checked+=1
            stat=[a+b for a,b in zip(stat,add,strict=True)]
    return stat,checked


def final(stat):return stat[0]-sum(stat[1:])


def assembly_gate():
    regions=((0,1),(2,3,4),(5,6,7))
    P={0,2,5}
    edges={(0,1),(2,3),(3,4),(5,6),(6,7)}
    cross={(1,3),(1,6),(4,7)}
    owner={v:i for i,reg in enumerate(regions) for v in reg}
    def adjacency(e):return {v:{b if a==v else a for a,b in e if v in (a,b)} for v in range(8)}
    local=adjacency(edges)
    full=adjacency(edges|cross)
    inputs=[]
    for R,J in zip(regions,((),(4,),(7,)),strict=True):
        inputs.append(histogram_local(permutations(R),P&set(R),list(J),local))
    truth=[0,0]
    original=0
    for word in permutations(range(8)):
        rank={v:i for i,v in enumerate(word)}
        peaks={v for v in range(8) if all(rank[v]<rank[u] for u in full[v])}
        z=len(peaks&P)
        overlap=int(rank[4]<rank[3] and rank[7]<rank[6])
        if len(peaks)==2:
            original+=1
            assert z+overlap<=2
        if z>2:continue
        reject=False
        if z==2:
            first=next((i for i,v in enumerate(word) if v not in P),8)
            u=word[first] if first<8 else None
            if first==0:reject=True
            if first==1:reject=word[0] not in full[u]
            if first==2:
                a,b=word[:2]
                reject=(owner[a]!=owner[b] and not full[u].intersection((a,b))) or (owner[a]==owner[b]!=owner[u])
            if first in (3,4):
                occ=Counter(owner[v] for v in word[:first])
                if len(occ)==2 and max(occ.values())<=2:
                    reject=occ[owner[u]]==0 or (occ[owner[u]]==1 and not full[u].intersection(word[:first]))
        if not reject:
            truth[0]+=1
            truth[1]+=z+overlap<=2
        if len(peaks)==2:assert not reject
    sum0,cell=count_all((2,3,3),inputs[0][0],inputs[1],inputs[2],False)
    sum1,_=count_all((2,3,3),inputs[0][0],inputs[1],inputs[2],True)
    assert truth==[final(sum0),final(sum1)] and truth[0]>truth[1]>=original
    return {'all_words':factorial(8),'true_two_peak_count':original,
            'before_joint_filter':truth[0],'after_joint_filter':truth[1],'mask_cells':cell}


def main():
    p=ArgumentParser();p.add_argument('batch',type=Path);p.add_argument('out',type=Path)
    p.add_argument('--prior',type=Path,default=R52);a=p.parse_args()
    reg_gates=[rectangle_gate(*z) for z in ((2,3),(3,3))]
    assemble=assembly_gate()
    cert=HERE/'pratt-61bit.json';all_primes,nodes=verified_primes(cert)
    primes=all_primes[:3]
    assert prod(primes)>factorial(42)
    receipts=read(a.batch/'receipts.json')
    assert receipts['counter_source_sha256']==sha(HERE/'joint_seam.cpp')
    assert receipts['plan_sha256']==sha(HERE/'n12-plan.json')
    reference=independent_sizes('666666')
    raws=[]
    for i,record in enumerate(receipts['jobs']):
        fp=a.batch/f'prime-{i}.json';fr=a.batch/f'prime-{i}.resource'
        assert record['prime_index']==i and record['waited'] and record['exit_status']==0
        assert record['stdout_sha256']==sha(fp) and record['resource_sha256']==sha(fr)
        raw=read(fp)
        assert int(raw['prime'])==primes[i] and raw['states_by_layer']==raw['processed_by_layer']==reference
        raws.append(raw)
    hist=exact_histograms(raws,primes)
    archive=read(a.prior/'r052-denominator-refined.json')
    T=old_chart(archive['t5_histograms']['12x5-B'])
    B=old_chart(archive['t5_histograms']['6x7-TR'])
    B_exact=compact_sum(hist)
    for q in range(3):
        size=max(len(B[q]),len(B_exact[q]))
        assert all(get(B[q],i)==get(B_exact[q],i) for i in range(size)),q
    assert all(x>=get(B[3],i) for i,x in enumerate(B_exact[3]))
    area=(60,42,42)
    r52=final(count_cells(area,[T,B,B],0))
    assert r52==int(archive['denominator'])
    new_no_joint,cell0=count_all(area,T,hist,hist,False)
    assert final(new_no_joint)==final(count_cells(area,[T,B_exact,B_exact],0))
    new_joint,cell1=count_all(area,T,hist,hist,True)
    assert final(new_joint)<final(new_no_joint)<=r52 and cell0==cell1
    assert 0<final(new_joint)
    upper=final(new_joint)
    old_combined=read(a.prior/'combined.json')
    tails=[]
    for item in old_combined['tail_lower_bounds']:
        count=int(item['combined_numerator']);ratio=Fraction(count,upper)
        tails.append({'strict_threshold':item['strict_threshold'],'r52_constructed_count':str(count),
                      'conditional_lower_exact':str(ratio),'conditional_lower_decimal':dec(ratio)})
    ledger={'scope':'Strict finite T5 upper bound using joint marked-seam peaks, not an exact A2 count',
            'r52_denominator_ledger_sha256':sha(a.prior/'r052-denominator-refined.json'),
            'r52_numerator_ledger_sha256':sha(a.prior/'combined.json'),'source_sha256':sha(Path(__file__)),
            'verified_pratt_nodes':nodes,'moduli_used':list(map(str,primes)),'exact_small_rectangles':reg_gates,
            'exact_small_three_block_gate':assemble,'worker_receipts_sha256':sha(a.batch/'receipts.json'),
            'r52_exact_upper':str(r52),'new_upper_with_exact_bottom_H_only':str(final(new_no_joint)),
            'new_upper_with_joint_seam_filter':str(upper),'new_upper_over_144factorial_exact':str(Fraction(upper,factorial(144))),
            'new_upper_over_144factorial_decimal':dec(Fraction(upper,factorial(144))),
            'old_over_new_upper':dec(Fraction(r52,upper)),'new_tail_bounds':tails,
            'mask_cells':cell1,'counts_before_joint_and_losses_R1_R2_R3different_R3same_R4_R5':list(map(str,new_no_joint)),
            'counts_after_joint_and_losses_R1_R2_R3different_R3same_R4_R5':list(map(str,new_joint)),
            'bottom_exact_histogram_F_T_V_H_per_mask':{str(m):[list(map(str,rows)) for rows in z] for m,z in hist.items()}}
    with a.out.open('x') as f:
        json.dump(ledger,f,indent=2);f.write('\n')
    print('new normalized upper',ledger['new_upper_over_144factorial_decimal'],'old/new',ledger['old_over_new_upper'])
    for t in tails:print(t['strict_threshold'],t['conditional_lower_decimal'])

if __name__=='__main__':main()
