# R78 far-face actual input: non-removable joint labels

Exactly one specified input was computed: S20 has 400 atoms
(21i/400,21j/400), 0<=i,j<20, each mass 1/400; the two outside atoms
are (-1/2,0) and (0,-1/2), both mass 37/100. No old atom is collapsed.
The region is Gamma=[-1/2,-1/400]^2. Halo witnesses are strict open squares
with captured mass >4r^2. This component does not change either previously
frozen R78 or R77 component.

## Caps and complete actual menu

The frozen R77 unequal-axis all-n cap lemma reduces all square caps to
20 left, 20 bottom and 400 joint prefix constraints, with the individual
outside atom caps. The computation checks all 440 rows exactly and also
the three pure-outside caps. All pass. The minimum joint outside capacity
is 1197/1600 >37/50. A new full-bbox enumeration is unnecessary and was
not performed. The original S20 positions are used in every row.

The strict halo uses r<1/2 by these caps. Its actual nonempty old capture
within Gamma is an origin-anchored k by l prefix, 1<=k,l<=10. We retain
all 100 prefixes times outside masks 0,1,2,3, and all three pure-outside
labels, for 403 full labels. The complete residual uses the same three
pure-outside labels. The separately owned R78 theory proof supplies the
complete-menu argument. This computation retains pure-old labels too.

As in the first frozen R78 component, exact 72-bit root enclosures produce
inner and outer center rectangles. The lower area is full-inner minus
residual-outer; the upper reverses both directions. Both exact integer
vertical-event interval union and compressed Cartesian midpoint coverage
evaluate every endpoint. They agree. The frozen R75 union core is imported
read-only with its SHA256 checked.

| Actual area quantity | Decimal display |
|---|---:|
| Full Gamma deletion | .016617159266386683 |
| Pure-old labels after complete residual | .00050625 |
| Exactly-one-outside group after complete residual | .002569342143534639 |
| Both-outside group after complete residual | .015337925518279819 |
| k,l>=2 joint group after complete residual | .015012389420974992 |
| Full deletion after removing all k,l>=2 joint labels | .004715488935277193 |
| Net area lost by removing those joint labels | .01190167033110949 |

These displays are rounded; the result JSON contains the exact rational
endpoints. The last row is certified directly as the multirow/multicolumn
joint union minus the union of all kept full labels and the complete residual.
Thus it already accounts for overlapping labels. The other group areas
are individual unions and cannot be added as independent credits.

Of the 100 labels for outside masks 0,1,2,3, respectively 1,21,21,76 have
nonempty open center rectangles in Gamma. All 81 possible k,l>=2 joint
labels are removed in the comparison, including those that are empty.
The main and auxiliary corner bounds both hold for this input: the full
certified upper endpoint is below 1/8 and 80003/640000. The point of the
example is that joint labels with both k,l>=2 cannot simply be discarded.

## Positive-area exact witness and exclusion certificate

Let E=[-61/1000,-59/1000]^2, of exact area 1/250000. The common radius
r=89/200 captures the actual 8x8 old prefix and both outside atoms for every
center in E. Their mass is 9/10, maximum distance over E is 441/1000,
the open geometric margin is 1/250, and the mass margin is 1079/10000.
The next old row and column are not captured by this radius.

For every one of the 322 kept labels (all labels except k,l>=2 joint), the
code computes its minimum possible farthest source distance over E. For
bbox [a,b]x[c,d] and E=[u,v]x[s,t], this is exactly
max{(b-a)/2,(d-c)/2,b-v,u-a,d-t,s-c,0}. Every kept label satisfies the strict
reverse inequality mass <4delta_E^2, with minimum slack 2721/250000.
This includes all thin joint labels (k=1 or l=1), all one-outside labels,
all pure-old labels and all pure-outside residual labels. Thus E is entirely
outside the kept halo and the complete residual halo, while entirely inside
the full halo. Its positive area certifies a strict net loss independently
of the larger directed numerical area enclosure.

## Reproduction and source pins

From the topic directory, existing Python 3.13.18:

```
/Users/zzh/.venvs/main/bin/python3 computation/round078-farface/certify_farface.py --input computation/inputs/r078-farface/farface-actual-input.json --out NEW-OUTPUT.json --bits 72
```

Use a new output path to preserve the frozen production result. The actual
run used the Python math helper; its receipt records exact command, runtime,
source hash, PID 5632 and exit 0. The manifest records observed process
absence. Exact deterministic arithmetic uses no random seed or numerical
tolerance. No parameter search, MATLAB, new agent, commit or push occurred.

The source adapts the already frozen R78 exact-prefix implementation into
a new file; it does not modify its source. Frozen dependencies:

- R75 union core `computation/round075/certify_three_atom_collapse.py`,
  SHA256 `a93cf5a038ae2579326e8c48170bfa102f9735cda982878fa75b95a5e91a4456`.
- R77 cap theorem `proofs/r077-mixed-only-obstruction-and-unequal-axis-caps.md`,
  SHA256 `3c71513ed72567445d4900b96f819f49d59481a7a6e9a1d1b4c1550cb4a211aa`,
  section 3. This is a mathematical dependency, not imported code.

This is one finite actual-input regional certificate. It gives no general
C2 improvement and no counterexample to the main corner bound.
