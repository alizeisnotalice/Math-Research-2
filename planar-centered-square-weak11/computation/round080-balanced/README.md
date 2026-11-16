# R80 actual balanced four-corner block with two outside atoms

Exactly one specified six-source input was computed. The old block is the
four corners of [0,999/1000]^2, each mass 1/4. The outside atoms are
q1=(-71/200,0), mass 1/4, and q2=(0,-171/200), mass 499/1000. All original
old positions are retained. This is an actual finite measure, not a collapsed
proxy and not an S20 uniform-grid input. No parameter search was run.

## Full qualification

The program enumerates all 63 nonempty actual source subsets, records their
tight bboxes, and checks mass <= max(1,bbox-span^2). For each bbox it also
captures all actual sources in that closed bbox and checks the resulting
full captured mass. These give 23 distinct bbox closures. Every subset and
closure cap passes exactly. The 12 eligible subsets satisfy the strict
test mass >bbox-span^2; all other subsets have empty open center receivers.

The unique maximum eligible subset is the old block, mask 15, mass 1.
It is an eligible four-source hyperedge because its span squared is
998001/1000000 <1. It connects all old sources even though its separated
two-source old pairs need not be eligible. The eligible origin-q1 and
origin-q2 pairs connect both outside atoms. The program independently
forms the connectivity relation from every eligible subset and obtains
one component. Total mass is 1749/1000, whole-source span squared is
859329/250000, and the nonterminal gap
span(all)^2-(35/33)mass(all)=395579/250000 is strictly positive.

## Complete strict halo and directed areas

For a source subset of mass M and bbox [a,b]x[c,d], its strict open center
rectangle is (b-sqrt(M)/2,a+sqrt(M)/2) x
(d-sqrt(M)/2,c+sqrt(M)/2). This rectangle is empty when M<=span^2.
All 63 full labels and all three outside-only residual labels are retained
in the catalog. The three full outside-only labels equal the complete
residual halo and cancel algebraically before interval bounding, so actual
deletion is computed from the 60 old-involving labels minus the same three
residual labels. This avoids artificial positive residual-root strips.

Every square root has a directed exact rational enclosure of width at most
1/(1000*2^72). Lower area uses old-full inner minus residual outer rectangles;
upper area reverses both directions. A vertical x-event/y-interval union
difference and a compressed Cartesian midpoint-cell coverage algorithm
give identical integer area endpoints. Both routines are imported read-only
from the frozen R75 core, whose SHA256 is checked before use.

Let e=1/1000, w=999/1000 and C=[-1/2,w+1/2]^2. All old-involving receiver
outer rectangles lie inside C. The four corner regions are
SW=[-1/2,-e]^2, SE=[w+e,w+1/2]x[-1/2,-e], and their northern reflections
about the old block. Noncorner means C minus these four corner squares.
Boundaries have area zero and the four corners have disjoint interiors.
Each directed total endpoint equals the sum of the four corresponding
corner endpoints and the noncorner endpoint, exactly in integer arithmetic.

| Actual deletion region | Area, decimal display |
|---|---:|
| Entire plane | 1.0913027092428533 |
| SW corner | .12991656304344726 |
| SE corner | .062001, exact 62001/1000000 |
| NW corner | .062001, exact 62001/1000000 |
| NE corner | .062001, exact 62001/1000000 |
| Noncorner | .7753831461994061 |

The exact rational endpoints in the result JSON establish SW>1299/10000>
1/8 and whole-plane deletion <273/250<7/2. The actual internal mass in the
SW source window (-1,1/2-e)^2 is 1/4: only the origin old atom lies there.
Thus this legal, unique-maximum, connected, nonterminal actual input
strictly violates the local auxiliary estimate Delta_SW<=m_P/2.
It remains far below the global deletion threshold 3.5.

The result also records the numerical per-corner allowance
1/8+3e^2/4=500003/4000000 and its SW comparison. That field is an additional
quota diagnostic, not a claim that the general problem imposes this quota
separately on every corner. In fact the total four-corner deletion is about
.3159195630, below the aggregate geometric allowance 1/2+3e^2=.500003.
The failed separate-corner auxiliary estimate does not refute that aggregate
route or the main global bound. No general C2 improvement is claimed.

## Actual run and reproduction

One run used existing Python 3.13.18 on macOS arm64. PID 19040 exited 0,
with helper wall time .083 seconds and script time .023 seconds. The final
manifest records observed process absence. The receipt preserves the exact
command, interpreter version and executed source hash. Seed and numerical
solver tolerance are inapplicable to this deterministic exact computation.

From the topic directory, using a new output path:

```
/Users/zzh/.venvs/main/bin/python3 computation/round080-balanced/certify_balanced.py --input computation/inputs/r080-balanced/six-atom-input.json --out NEW-OUTPUT.json
```

Frozen dependency: `computation/round075/certify_three_atom_collapse.py`,
SHA256 `a93cf5a038ae2579326e8c48170bfa102f9735cda982878fa75b95a5e91a4456`.
This component was produced independently of the already frozen R79 payload.
Its final mathematical status awaits independent exact replay and the
separate local-geometry proof. There was no MATLAB, new agent, commit or
push; R80 archival integration remains the root agent's responsibility.
