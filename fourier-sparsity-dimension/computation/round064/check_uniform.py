"""All 4-variable truth tables and linear flats; verify the deterministic rank-density transfer."""
import json
import platform
from collections import Counter
from functools import lru_cache

def spaces(n):
    known={1}
    for a in range(1,1<<n):
        old=list(known)
        for H in old:
            if not H>>a&1:
                shifted=sum(1<<(v^a) for v in range(1<<n) if H>>v&1)
                known.add(H|shifted)
    return sorted(known)

def dim(H): return H.bit_count().bit_length()-1

def span_mask(mask,n=4):
    H=1
    for a in range(1,1<<n):
        if mask>>a&1 and not H>>a&1:
            H |= sum(1<<(x^a) for x in range(1<<n) if H>>x&1)
    return H

def run():
    n=4; N=16
    subs=spaces(n)
    signs=[sum(1<<x for x in range(N) if (x&a).bit_count()%2) for a in range(N)]
    freq=Counter()
    for f in range(1<<N):
        mask=sum(1<<a for a,m in enumerate(signs) if N!=2*(f^m).bit_count())
        freq[mask]+=1
    opt_count=0
    transferred=0
    nonzero_core=0
    records=0
    maxratio={}
    for mask in freq:
        s=mask.bit_count();d=dim(span_mask(mask))
        if s not in maxratio or d>maxratio[s]:
            maxratio[s]=d
    for mask,mult in freq.items():
        s=mask.bit_count()
        space=span_mask(mask)
        d=dim(space)
        if not d:
            continue
        local=[T for T in subs if T&~space==0]
        merit=lambda H:d*(mask&H).bit_count()-2*s*dim(H)
        value=max(map(merit,local))
        record=all(d*d*ss>=dd*dd*s for ss,dd in maxratio.items() if ss<=s)
        if record:
            records+=mult
        for H in local:
            if merit(H)!=value:
                continue
            opt_count+=1
            nonzero_core+=dim(H)>0
            for T in local:
                if T&H==H:
                    assert d*(mask&(T^H)).bit_count()<=2*s*(dim(T)-dim(H))
            # The identity is an admissible injective kernel compression whenever
            # h <= k, and all these finite spectra have h <= k.
            h=dim(H)
            if h and d*((1<<h)-1)>2*s*h:
                continue
            for T in local:
                # Contributions in the two disjoint parts, including rank identities.
                j=dim(T&H)
                u=dim(T)-j
                assert d*(mask&T&~H).bit_count()<=2*s*u
                assert d*(mask&T&H&~1).bit_count()<=2*s*j
                assert d*(mask&T&~1).bit_count()<=2*s*dim(T)
                transferred+=1
    return dict(python=platform.python_version(),ambient=n,linear_spaces=len(subs),truth_tables=1<<N,
                support_masks=len(freq),all_optimizer_checks=opt_count,positive_dimension_optimizers=nonzero_core,
                finite_class_record_truth_tables=records,all_transferred_flat_checks=transferred,
                note="Finite identities only. Records are relative to n<=4, and no asymptotic records are assumed.")
if __name__=="__main__":
    print(json.dumps(run(),indent=2))
