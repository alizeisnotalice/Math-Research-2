from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path
import platform
import sys

P = Path(__file__).resolve().parent
A2 = {3:143808,4:779299174080,5:7063261481362958280000}

def check(n):
    data = P/f'n{n}-contact.tsv'
    projected=Counter()
    groups=defaultdict(Counter)
    by_b=defaultdict(Counter)
    for row in data.read_text().splitlines():
        i,b,h,d,D,count=map(int,row.split())
        assert 1<=b<=8 and b+2<=i<n*n and d>=2 and D>=1
        assert 1<=h<=n*n-i-D+1 and count>0
        groups[i,b,d,D][h]+=count
        projected[i,b,h,d]+=count
        by_b[b][i,h,d]+=count
    comparisons=0
    for (i,b,d,D), counts in groups.items():
        m=n*n-i
        total=sum(counts.values())
        remaining=total
        pop=comb(m,D)
        for t in range(m-D+2):
            if t:
                remaining-=counts[t]
            assert remaining*pop <= comb(m-t,D)*total
            comparisons+=1
        mean=sum(v*comb(h+b-1,b) for h,v in counts.items())
        assert mean*comb(D+b,b) <= total*comb(m+b,b)
    if n in (3,4):
        old = json.loads((P.parents[1]/'r048/merge-repair/verification.json').read_text())['exact']
        expected={(i,b,h,d):v for i,b,h,d,v in old[n-3]['histogram'] if b<=8}
        assert dict(projected)==expected
    total=sum(projected.values())
    assert total<=A2[n]
    if n==5:
        old=json.loads((P.parents[2]/'computation/results/singleton-n5-r003.json').read_text())
        assert sum(by_b[1].values())==old['beta_positive_count']
        if len(sys.argv)>1:
            comparison=[list(map(int,x.split())) for x in Path(sys.argv[1]).read_text().splitlines()]
            assert [[*row,c] for row,c in sorted(projected.items())]==comparison
    table=[]
    for b,counts in sorted(by_b.items()):
        amount=sum(counts.values())
        rank_moment=sum(v*Fraction(comb(h+b-1,b),comb(i-1,b))
                        for (i,h,d),v in counts.items())/amount
        table.append(dict(b=b,words=amount,conditional_b_given_K2=str(Fraction(amount,A2[n])),
                          conditional_h1_given_b=str(Fraction(sum(v for (i,h,d),v in counts.items() if h==1),amount)),
                          h_rising_over_i_falling_mean=str(rank_moment)))
    return dict(n=n,source_sha256=hashlib.sha256((P/'low_branch.cpp').read_bytes()).hexdigest(),
                data_sha256=hashlib.sha256(data.read_bytes()).hexdigest(),
                cells=len(data.read_text().splitlines()),projected_cells=len(projected),
                groups=len(groups),tail_comparisons=comparisons,covered_words=total,
                covered_fraction=str(Fraction(total,A2[n])),by_b=table)

if __name__=='__main__':
    print(json.dumps(dict(python=platform.python_version(),checks=[check(n) for n in (3,4,5)]),indent=2))
