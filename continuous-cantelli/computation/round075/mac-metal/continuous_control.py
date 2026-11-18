"""Continuous digit-Cantor exact low-mode control and sharpened error guards."""
import json
from fractions import Fraction as F
from pathlib import Path
import sys

import mpmath as mp

p,q=7,11
B=5*p*q
left=[1+5*q*j for j in range(p)]
right=[3+5*p*j for j in range(q)]
E=sorted(left+right)
D=[i for i in range(B) if i not in E]
assert set(left).isdisjoint(right)
assert {0,B-1}.issubset(D)
assert not any((x+1)%B in E for x in E)
assert all(n%p and n%q for n in range(1,7))
# Largest gaps share width 1/B. At the symmetry-breaking n=p,q,
# compare the primary integral to its exact moment bound on the tail.
guards={}
for n in (p,q):
    cos_bound=1-F(8*n*n,B*B)   # cos(pi*n/B) >= 1-8n^2/B^2 using pi<4.
    loss_bound=1-F(n*n,8*B)    # At alpha=1 and height<=1.
    relative_tail=F(len(E),B*B-len(D))
    assert cos_bound*loss_bound>relative_tail
    guards[str(n)]={'scaledPrimaryLower':str(cos_bound*loss_bound),
                     'scaledTailUpper':str(relative_tail),
                     'strictMargin':str(cos_bound*loss_bound-relative_tail)}
with mp.workdps(65):
    rows=[]
    g=mp.mpf(1)/B
    for n in (*range(1,7),p,q):
        amplitude=mp.mpf(n*n)/2
        event=sum(mp.exp(2j*mp.pi*n*d/B) for d in E)
        gap=2*mp.quad(lambda t:-mp.expm1(-amplitude*t)*mp.cos(2*mp.pi*n*(t-g/2)),[0,g/2])
        c=-mp.exp(1j*mp.pi*n/B)*event*gap
        rows.append({'n':n,'omittedDigitPhaseMagnitude':mp.nstr(abs(event),34),
                     'firstLevelHarmonicMagnitude':mp.nstr(abs(c),35)})
        if n<=6: assert abs(event)<mp.mpf('1e-62')
        else: assert abs(c)>mp.mpf('1e-8')
out={'base':B,'omittedDigits':E,'retainedCount':len(D),'exactZeroModes':list(range(1,7)),
     'detectingModes':[p,q],'alphaOneHeightAtMostOneGuards':guards,'numericChecks':rows}
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
print('base',B,'prime cycles',p,q,'rational detection guards',guards)
