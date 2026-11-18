"""python3 -B check_support_capacity.py; rational linear-cost checks, no dependencies."""
from fractions import Fraction
from itertools import product
import json
import platform
from check_unrestricted import alt, columns, min_cover, rank, rows


def test(catalog,a,b,alternating):
    tight=total=0
    for A,B in product(catalog,repeat=2):
        matrices=(A,B,A^B)
        allcols=[columns(M,a,b) for M in matrices]
        ranks=list(map(rank,allcols))
        image=rank(v for group in allcols for v in group)
        if alternating:
            assert all(r%2==0 for r in ranks)
            cost=[r//2 for r in ranks]
            target=image
        else:
            cost=ranks
            target=image+rank(v for M in matrices for v in rows(M,a,b))
        minimum=min_cover(cost)
        assert minimum>=target
        tight+=minimum==target
        total+=1
    return {'all':total,'equality':tight}


def full_rows():
    output=[]
    for n in range(2,13):
        N=1<<n
        cut=max(Fraction(n+2,N),*(Fraction(c+1,N-(N>>c)) for c in range(1,n+1)))
        assert cut==Fraction(n+2,N)
        output.append({'n':n,'b_plus_squared':str(cut**2*N)})
    return output


def addresses():
    output=[]
    for m in range(1,11):
        N=1<<m
        weight=max(Fraction(N+m-q-j+1,N*N-q*(1<<j))
                  for q in range(N+1) for j in range(m+1) if q<N or j<m)
        assert weight==Fraction(2,N)
        output.append({'N':N,'b_plus_squared':str(weight**2*N*N)})
    return output


if __name__=='__main__':
    result={'python':platform.python_version(),'tolerance':0,'seed':None,
            'exhaustive_rectangles':{},'exhaustive_alternating':{}}
    for a,b in ((1,4),(2,2),(2,3),(3,2)):
        result['exhaustive_rectangles'][str((a,b))]=test(range(1<<(a*b)),a,b,False)
    for n in (3,4):
        mats=[alt(mask,n) for mask in range(1<<(n*(n-1)//2))]
        result['exhaustive_alternating'][str(n)]=test(mats,n,n,True)
    result['full_supports']=full_rows()
    result['standard_addresses']=addresses()
    print(json.dumps(result,indent=2))
