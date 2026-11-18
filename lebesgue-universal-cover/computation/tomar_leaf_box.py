"""Decode one 11D preorder leaf using exact rational endpoints."""
from fractions import Fraction
import numpy as np

def leaf_box(prefix,codes,part,idx,meta):
    lows=list(map(Fraction,meta["lo"]))
    highs=list(map(Fraction,meta["hi"]))
    for k,bit in enumerate(np.unravel_index(part,(2,)*11)):
        mid=(lows[k]+highs[k])/2
        if bit: lows[k]=mid
        else: highs[k]=mid
    i=0
    while i!=idx:
        k=int(codes[i]);assert k<11
        mid=(lows[k]+highs[k])/2
        # The end of the left subtree is the first return to the slot count
        # before its parent was parsed. If it precedes idx, take the right child.
        returned=np.flatnonzero(prefix[i+2:idx+1]==prefix[i])
        if len(returned):
            lows[k]=mid
            i+=2+int(returned[0])
        else:
            highs[k]=mid
            i+=1
    center=[(a+b)/2 for a,b in zip(lows,highs)]
    widths=[(b-a)/2 for a,b in zip(lows,highs)]
    c=np.array(center,dtype=float); h=np.array(widths,dtype=float)
    assert list(map(Fraction,c))==center and list(map(Fraction,h))==widths
    return c,h
