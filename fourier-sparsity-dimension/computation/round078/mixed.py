from functools import cache
from collections import Counter
from fractions import Fraction
import json
from audit import (all_subspaces,basis_label_subset,character,indicator_bands,
                   quotient_bins,r49_cells,rank_space,weak_folding_cells)

def check(make_cells):
    address,cells=make_cells(1)
    coeff,N,bands=indicator_bands(cells,address)
    n=address+N
    A=(1<<address)-1
    @cache
    def active(entries):
        return sum(sum(v*character(x,z) for x,v in entries)!=0 for z in range(1<<n))
    checked=Counter()
    for level in all_subspaces(n):
        for H in level:
            bins=quotient_bins(coeff,H)
            weights={q:active(tuple(entries)) for q,entries in bins.items()}
            labels=sorted(weights,key=lambda q:-weights[q])
            chosen=basis_label_subset(labels)
            tau=Fraction(sum(weights.values())-sum(weights[i] for i in chosen),1<<n)
            h0=rank_space({v for v in H if v<=A})
            h1=rank_space(H)-h0
            if make_cells==r49_cells:
                assert (4*tau+4)*N*(h0+4)**2>=max(0,N-2*h1)**3
            else:
                assert N*(4*tau+4)>=max(0,N-2*h1-1)*(1<<(address-h0))
            outside=[i for i in range(N) if all(not ((h>>(address+i))&1)
                       or (h>>address).bit_count()>2 for h in H)]
            target=sum(Fraction(m*(m-1-(m.bit_length()-1)),N) for i in outside
                       for m in [N//len(H&bands[i])])
            assert tau>=target
            checked[rank_space(H)]+=1
    return dict(checked)
print(json.dumps({f.__name__:check(f) for f in (r49_cells,weak_folding_cells)},indent=2))
