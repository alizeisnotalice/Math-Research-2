#!/usr/bin/env python3
"""Finite exact cross-check; general conclusions require the binomial proof."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import platform

pairs = events = 0
for q in range(2, 101):
    for a in range(1, q):
        b = q-a
        total = comb(q, a)
        top2_one_orientation = comb(q-2, a-1)
        high_a, high_b = comb(q-1, a-1), comb(q-1, b-1)
        assert high_a + high_b == total and F(high_a,total) == F(a,q)
        assert F(high_a,top2_one_orientation) == F(q-1,b)
        assert F(total,2*top2_one_orientation) == F(q*(q-1),2*a*b)
        assert F(a,q)*F(q,a)+F(b,q)*F(q,b) == 2
        for delta in (F(1,10),F(1,4),F(1,3),F(49,100)):
            mass = F(a,q)*(F(a,q)<=delta)+F(b,q)*(F(b,q)<=delta)
            assert mass <= delta
            events += 1
        for k in (1,2,3,5,10,20):
            assert a*(a<=k)+b*(b<=k) <= 2*k
        pairs += 1
out = dict(kind='finite-binomial-ordered-first-merge-fibre-audit',
           python=platform.python_version(),q_range=[2,100],fibres=pairs,
           endpoint_events=events,
           checks=['orientation a/q','ordered density (q-1)/b',
                   'symmetric density q(q-1)/(2ab)',
                   'conditional reciprocal branch fraction mean 2',
                   'high-root endpoint probability <= delta for delta<1/2',
                   'high-root small-branch distance moment <=2k'],
           scope='exact finite index checks; general proof is binomial identity, not finite extrapolation')
Path(__file__).with_name('finite-fibre-audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
