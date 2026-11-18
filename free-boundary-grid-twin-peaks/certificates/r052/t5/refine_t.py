"""CRT reconstruction and disjoint fourth/fifth-label exclusions."""
from collections import Counter
from functools import cache
from fractions import Fraction
from hashlib import sha256
from itertools import permutations, product
from math import factorial, prod
from pathlib import Path
import argparse
import json
import sys

r051 = Path(__file__).resolve().parents[2]/"r051/n12"
if r051.exists():
    sys.path.insert(0, str(r051))

from count_rectangles import dec, merge, readback, verified_primes

ROOT = Path(__file__).resolve().parent


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def two_safe_lower(chart, k):
    f,t,v=chart
    return max(0, f[k]-t[k]-v[k]-(t[k-1] if k else 0))


def exclusions(areas, charts, exact=None):
    blocks=range(len(areas))
    singles=[c[0][1]-c[1][1] for c in charts]
    unsafe=[c[2][1] for c in charts]
    unseen=[c[1][0] for c in charts]
    zero=[c[0][0] for c in charts]
    doubles=[two_safe_lower(c,1) for c in charts] if exact is None else [h[1] for h in exact]

    @cache
    def shuffles(used):
        left=tuple(a-x for a,x in zip(areas,used,strict=True))
        return factorial(sum(left))//prod(factorial(x) for x in left) if min(left)>=0 else 0

    # First two labels can both be protected in a single block.
    # A third label in a different block has no earlier cross-neighbor.
    rank3=0
    for order in product(range(3),repeat=len(areas)):
        if sum(order)!=2:
            continue
        for i in blocks:
            if not order[i]:
                continue
            di=two_safe_lower(charts[i],order[i]) if exact is None else exact[i][order[i]]
            for k in blocks:
                if i==k:
                    continue
                remaining=shuffles(tuple(2*(a==i)+(a==k) for a in blocks))
                rank3 += remaining * di * charts[k][1][order[k]] * prod(
                    charts[a][0][order[a]] for a in blocks if a not in (i,k))
    result=[rank3]
    for length, patterns in ((3,((2,1),(1,2))),(4,((2,2),))):
        count=0
        for i in blocks:
            for j in range(i+1,len(areas)):
                for x,y in patterns:
                    used={i:x,j:y}
                    multiplicity=factorial(length)//(factorial(x)*factorial(y))
                    before=[singles[a] if used[a]==1 else doubles[a] for a in (i,j)]
                    for k in blocks:
                        remaining=shuffles(tuple(used.get(a,0)+(a==k) for a in blocks))
                        if k not in used:
                            count += multiplicity * remaining * prod(before) * unseen[k] * prod(
                                zero[a] for a in blocks if a not in (i,j,k))
                        elif used[k]==1:
                            other=j if k==i else i
                            count += multiplicity * remaining * doubles[other] * unsafe[k] * prod(
                                zero[a] for a in blocks if a not in used)
        result.append(count)
    return result,doubles


def result(areas, charts):
    top_three=merge(areas,charts)[-1]
    extra,doubles=exclusions(areas,charts)
    upper=top_three-sum(extra)
    nwords=factorial(sum(areas))
    assert 0<upper<=top_three
    return {'areas':areas,'locally_protected_prefix_lower_at_one_peak':list(map(str,doubles)),
            'extra_excluded_rank3_rank4_rank5':list(map(str,extra)),'top_three_upper':str(top_three),
            'top_five_upper':str(upper),'top_three_probability':dec(Fraction(top_three,nwords)),
            'top_five_probability':dec(Fraction(upper,nwords)),'ratio_after_top_five':dec(Fraction(upper,top_three))}


def gate():
    regions=[(0,1,2),(3,4,5),(6,7)]
    owner={v:i for i,R in enumerate(regions) for v in R}
    safe={0,1,3,4,6}
    local_edges={(u,v) for R in regions for u,v in zip(R,R[1:])}

    def adjacency(edges):
        return {v:{y if x==v else x for x,y in edges if v in (x,y)} for v in range(8)}

    local=adjacency(local_edges)
    full=adjacency(local_edges|{(2,5),(5,7)})
    charts=[]
    exact=[]
    for R in regions:
        f,t,v,h,swap=([0]*3 for _ in range(5))
        for word in permutations(R):
            rank={v:i for i,v in enumerate(word)}
            z=sum(all(rank[u]<rank[x] for x in local[u]) for u in R if u in safe)
            a,b=word[:2]
            f[z]+=1
            t[z]+=a not in safe
            v[z]+=a in safe and b not in safe and b not in local[a]
            h[z]+=a in safe and b in safe
            swap[z]+=a not in safe and b in safe and b in local[a]
        assert all(h[k]==f[k]-t[k]-v[k]-(swap[k-1] if k else 0) for k in range(3))
        charts.append([f,t,v])
        exact.append(h)
    truth=[0,0,0]
    for word in permutations(range(8)):
        rank={v:i for i,v in enumerate(word)}
        peaks={v for v in range(8) if all(rank[v]<rank[u] for u in full[v])}
        if len(peaks & safe)!=2:
            continue
        a,b,c=word[:3]
        if a in safe and b in safe and c not in safe and owner[a]==owner[b]!=owner[c]:
            assert c in peaks and len(peaks)>=3
            truth[0]+=1
        for index,length in enumerate((3,4),1):
            prefix=word[:length]
            if not all(v in safe for v in prefix) or word[length] in safe:
                continue
            counts=Counter(owner[v] for v in prefix)
            if len(counts)!=2 or max(counts.values())>2:
                continue
            v=word[length]
            if counts[owner[v]]==0 or (counts[owner[v]]==1 and not full[v].intersection(prefix)):
                assert v in peaks and len(peaks)>=3
                truth[index]+=1
    areas=list(map(len,regions))
    assert exclusions(areas,charts,exact)[0]==truth and min(truth)>0
    bound=exclusions(areas,charts)[0]
    assert all(a<=b for a,b in zip(bound,truth,strict=True))
    return {'full_words':factorial(8),'exact_rank3_rank4_rank5':truth,'recovered_lower':bound}


def shards(old, root):
    primes,nodes=verified_primes(root/'pratt-61bit.json')
    expected={(12,5,'B'):old['deferred_plan_only_regions'][1],(6,7,'TR'):old['deferred_plan_only_regions'][2]}
    assert prod(primes)>factorial(60) and prod(primes[:4])>factorial(42)
    raw={}
    receipts={}
    for site in ('local','macbook'):
        directory=root/site
        report=json.loads((directory/'finished.json').read_text())
        assert report['waited_all']
        assert report['source_sha256']==digest(root/'rect_protected.cpp')
        assert report['certificate_sha256']==digest(root/'pratt-61bit.json')
        for saved in report['jobs']:
            w,h,e,i=saved['job']
            shape=w,h,e
            stem=f'{w}x{h}-{e}-p{i}'
            result=directory/f'{stem}.json'
            resource=directory/f'{stem}.resource'
            assert saved['waited_reaped'] and saved['returncode']==0
            assert (digest(result),digest(resource))==(saved['stdout_sha256'],saved['resource_sha256'])
            z=json.loads(result.read_text())
            assert (z['width'],z['height'],z['unsafe_edges'],z['prime_index'])==(w,h,e,i)
            assert int(z['prime'])==primes[i]
            assert z['processed_by_layer']==z['states_by_layer']==expected[shape]['states_by_layer']
            assert z['max_table_bytes']==expected[shape]['max_table_bytes']==saved['planned_array_bytes']
            assert (shape,i) not in raw
            raw[shape,i]=z
        receipts[site]=report
    assert len(raw)==9
    top=readback((12,5,'B'),raw,primes)
    bottom=readback((6,7,'TR'),raw,primes[:4])
    return top,bottom,receipts,nodes


def main():
    p=argparse.ArgumentParser()
    p.add_argument('r051_denominator',type=Path)
    p.add_argument('r051_numerator',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--shards',type=Path,default=ROOT)
    args=p.parse_args()
    old=json.loads(args.r051_denominator.read_text())
    numer=json.loads(args.r051_numerator.read_text())
    top,bottom,hosts,nodes=shards(old,args.shards)
    hist={name:[[int(x) for x in line] for line in rows] for name,rows in old['regional_histograms'].items()}
    inputs={'quadrants':([36]*4,[hist['6x6-BR']]*4),
            'T4':([48]*3,[hist['12x4-B'],hist['6x8-TR'],hist['6x8-TR']]),
            'T5':([60,42,42],[top,bottom,bottom])}
    layouts={name:result(*configuration) for name,configuration in inputs.items()}
    assert all(layouts[key]['top_three_upper']==old['computed_results'][name]['top_three_upper']
               for key,name in [('quadrants','quadrants'),('T4','three_block_T')])
    key=min(layouts,key=lambda name:int(layouts[name]['top_five_upper']))
    denominator=int(layouts[key]['top_five_upper'])
    lower=[]
    oldden=int(old['computed_results']['three_block_T']['top_three_upper'])
    for row in numer['thresholds']:
        value=int(row['sum_parts']['total'])
        lower.append({'strict_threshold':row['strict_threshold'],'count':str(value),
                      'lower_exact':str(Fraction(value,denominator)),
                      'lower_decimal_rounded':dec(Fraction(value,denominator))})
    protocol={'gate':gate(),'layouts':layouts,'selected_layout':key,
              'denominator':str(denominator),'old_denominator_over_new':dec(Fraction(oldden,denominator)),
              'r051_denominator_sha256':digest(args.r051_denominator),'r051_numerator_sha256':digest(args.r051_numerator),
              'reused_r051_numerator':lower,'t5_histograms':{'12x5-B':top,'6x7-TR':bottom},
              'pratt_nodes_checked':nodes,'actual_batch_receipts':hosts}
    with args.output.open('x') as stream:
        json.dump(protocol,stream,indent=2)
        stream.write('\n')


if __name__=='__main__':
    main()
