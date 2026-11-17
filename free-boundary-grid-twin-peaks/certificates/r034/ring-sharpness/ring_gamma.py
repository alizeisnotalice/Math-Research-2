#!/usr/bin/env python3
"""Exact Gamma and its all-nonroot extension for a corner L-infinity word."""
from fractions import Fraction
import sys
import json

# Fractions at n=128 have numerators with more than Python's default
# decimal-conversion guard permits.  The calculation remains exact.
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

def word(n):
    out=[0]
    for k in range(1,n):
        out += [y*n+k for y in range(k+1)]       # right side, bottom-up
        out += [k*n+x for x in range(k-1,-1,-1)] # top side, right-left
    assert len(out)==n*n==len(set(out))
    return out

def gamma(n):
    w=word(n);rank={v:i+1 for i,v in enumerate(w)};N=n*n
    nb=lambda v:[u for u in (v-1 if v%n else -1,v+1 if v%n<n-1 else -1,v-n if v>=n else -1,v+n if v<N-n else -1) if u>=0]
    early={v:[u for u in nb(v) if rank[u]<rank[v]] for v in w}
    l={v:min([rank[u] for u in early[v]],default=N+1) for v in w}
    b={p:N+1 for p in w}
    unique_children={p:[] for p in w}
    for child in w[1:]:
        if len(early[child])==1:
            p=early[child][0]
            b[p]=min(b[p],rank[child])
            unique_children[p].append(child)
    terms={p:Fraction(1,b[p]-l[p]-1) for p in w[1:] if l[p]<=N and b[p]>l[p]+1}
    # The proof's Gamma excludes neighbours of the root.  Retain terms for
    # every non-root vertex as a nonnegative, deliberately larger quantity.
    root_neighbors=set(nb(w[0]))
    eligible={p:q for p,q in terms.items() if p not in root_neighbors}
    return (sum(terms.values(),Fraction()), sum(eligible.values(),Fraction()),
            terms,eligible,w,rank,early,b,l,unique_children)

def expected_unique_children(n):
    """The two singleton-parent children on each nonzero L-infinity ring."""
    ans=[]
    for k in range(1,n):
        ans += [k, k*n+k]  # (k,0), (k,k)
    return ans

rows=[]
for n in (8,16,32,64,128):
    tilde,g,t,eligible,w,rank,early,b,l,uc=gamma(n)
    observed=sorted(c for cs in uc.values() for c in cs)
    assert observed==sorted(expected_unique_children(n)), (n, observed, expected_unique_children(n))
    # (1,0) is the parent of both (1,1) and (2,0); all later special
    # parents are (k,0) and (k,k-1), as used in the written bound.
    assert sorted(p for p,cs in uc.items() if cs)==sorted([0,1]+list(range(2,n-1))+[(k-1)*n+k for k in range(2,n)])
    harmonic=sum((Fraction(1,j) for j in range(1,n)), Fraction())
    # Exact verification of the elementary bound derived in REPORT.md.
    upper=Fraction(1)+Fraction(11,3)*harmonic
    assert tilde <= upper
    assert g <= tilde
    rows.append({'n':n,
                 'eligible_Gamma_float':float(g),
                 'all_nonroot_tilde_Gamma_float':float(tilde),
                 'eligible_terms':len(eligible), 'all_nonroot_terms':len(t),
                 'unique_early_children':len(observed),
                 'upper_float':float(upper),
                 'exact_fraction_checks_passed':True})
print(json.dumps({'construction':'corner L-infinity rings, right side bottom-up then top side right-left',
                  'Gamma_definition':'sum excludes the two root-neighbours; tilde sums every non-root vertex',
                  'rows':rows},indent=2))
