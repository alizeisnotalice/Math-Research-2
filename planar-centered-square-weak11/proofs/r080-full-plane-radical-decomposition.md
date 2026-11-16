# R80: exact full-plane deletion, including the cross-old narrow rectangle

2026-11-16. This note concerns the single actual six-source input in the
R80 balanced-block certificate. It does not change a general weak bound.
The original four old corners have equal mass 1/4 and side w=999/1000.
The two outside atoms are (-d,0) of mass 1/4 and (0,-D) of mass 499/1000,
where d=71/200 and D=171/200. Set

    e=1/1000, a=sqrt(1/8), U=sqrt(999/4000), V=sqrt(749/4000),
    s=sqrt(499/4000).

All coordinates and thresholds refer to one common radius and strict halo
threshold 1. The independent finite certificate checks all 63 nonempty
source subsets and all 23 source-bbox closures. Its nine eligible labels
involving old mass are precisely masks

    1, 2, 4, 8, 15, 17, 33, 35, 49,

with bit indices ordered as origin, (w,0), (0,w), (w,w), (-d,0), (0,-D).
The remaining three eligible labels are exactly the complete residual.
This finite exhaustive statement is supplied by the raw-source certificate,
not inferred by only drawing the most visible rectangles.

## 1. The complete geometry

The origin-containing three-source cluster has four old-bearing receiving
rectangles in the original center coordinates:

    O=[-1/4,1/4]^2,
    A=[-a,a-d] x [-a,a],
    B=[-V,V] x [-V,V-D],
    T=[-U,U-d] x [-U,U-D].

The complete residual consists of

    R1=[-d-1/4,-d+1/4] x [-1/4,1/4],
    R2=[-s,s] x [-D-s,-D+s],
    R12=[-V,V-d] x [-V,V-D].

Closed endpoints here are for Lebesgue area only. The underlying strict
receivers are open, and no source-boundary atom is removed from a capacity
test. All radical comparisons in the following partition are strict.

The old-bearing labels outside this local cluster are the other three
old singleton squares (total area 3/4), the all-old central square
[w-1/2,1/2]^2 (area e^2), and the label of mask 35:

    N=[w-U,U] x [-U,U-D].

The last label contains the old bottom edge and the far outside atom.
Its mass is 999/1000, strictly larger than w^2=998001/1000000; it is
eligible and has positive area

    |N|=(2U-w)(2U-D).

It is disjoint, up to boundaries, from all the other listed old-bearing
receivers and from the residual. Indeed w-U>V while U<w-1/4; its y interval
is negative, whereas the central square and the upper old squares have
positive y. All other old singleton squares and the central square are
likewise pairwise disjoint, lie outside the residual, and do not meet the
origin cluster. For example the central square starts at .499 on both
axes, while the cluster has x or y at most V<.433 or a<.354, and the
right/bottom old square starts at x=.749.

An initial hand-derived candidate omitted N. That candidate was never
certified: cross-review identified mask 35 before this note and the
symbolic run were frozen. The full producer and independent area replay
included this label from the start. The formula below includes its positive
contribution explicitly.

## 2. Two disjoint pieces of the origin cluster

The inequalities D-U>a and D-s>U show that O union A is separated from
B union T, and that R2 lies strictly below every old-bearing rectangle.
R1 meets only O union A; R12 is contained in B intersection T and does
not meet the upper piece. Hence there is one residual subtraction per piece.

The overlap of O with A has area (a-d+1/4)/2. The intersection of R1 with
O union A has this same area: the part of R1 in O already lies in A.
Consequently

    |(O union A) minus R1|
      =1/4 + (2a-d)2a - (a-d+1/4)
      =D-2Da.

The part of B outside T is its right strip, of width V-U+d and height
2V-D. Since R12 is contained in T, the lower piece has area

    |(B union T) minus R12|
      =(2U-d)(2U-D)+(V-U+d)(2V-D)-(2V-d)(2V-D)
      =(2U-d)(2U-D)+(2d-U-V)(2V-D).

All these regions are actual original receivers, not moved-source proxies.
Adding the disjoint old singleton, central, cross-old and cluster pieces
gives the exact full-plane identity

    |H(mu) minus H(nu)|
      =3/4+e^2+D-2Da
        +(2U-d)(2U-D)+(2d-U-V)(2V-D)
        +(2U-w)(2U-D)
      =3779121/1000000 - (171/100)a - (5273/1000)U
        +(91/40)V - 2UV.

Its decimal display is approximately 1.0913027092428533. This full deletion
is strictly below 273/250 and hence below 7/2. The separately proved actual
SW corner deletion exceeds 1/8. The two assertions coexist: a failed local
half-mass quota is not a counterexample to the total deletion budget.

## 3. Verification and reproduction

The formal local-corner proof is `r080-balanced-corner-halfmass-obstruction.md`.
The full producer result is
`computation/results/round080-balanced/six-atom-result.json`, SHA256
`f22e7216fe6a39e36e991fd94a24b00047320376de257fc258e3ebd553f0c9a3`.
The separate independent review reconstructs its complete labels and areas
from the raw six sources, using a horizontal union algorithm.

Root's `computation/round080-root/check_full_plane.wls` checks the displayed
algebraic identities and strict geometric orderings with exact radicals.
It compares both the full expression and the author's SW expression to
the frozen producer's directed rational endpoints, by exact algebraic
signs. It does not enumerate another halo or verify all caps anew. Its
actual version, process IDs, exit status and outputs are stored separately.

From the topic directory, reproduce into a fresh output path:

```
wolframscript -file computation/round080-root/check_full_plane.wls computation/results/round080-balanced/six-atom-result.json NEW-OUTPUT.json
```

The analytic partition and exhaustive source menu carry the geometry;
symbolic agreement alone would not establish that a receiver had not been
omitted. This is why the mask-35 check precedes the final formula.
