# R69: a positive net mixed increment larger than one fifth

2026-11-15. Finite positive atomic measures, threshold 1, strict centered
axis-parallel square witnesses. Write I=H(S), R=H(A\S),
P=H(A)\(I union R), and D=H(A)\R. Since I is contained in H(A),

    |D| = |I| - |I intersect R| + |P|.                 (1)

The following potential shortcut is **false**, even for a connected,
nonterminal input with a unique maximum eligible block, two thick axes, and
minimum support-tree cost greater than 5L/2:

    |P| - |I intersect R| <= L/5.                     (2)

Had (2) held generally, the compact-support estimate B<3.3 applied to S
would have given |D|<3.5L. This note gives a strict counterexample to (2).
It does not refute the 3.5 deletion or W-prime selection candidates. Indeed,
the same example is paid by the new heavy-core criterion at the end.

## 1. A seven-source input in the support-tree remaining domain

Let S consist of the five positive sources

| Source | Position | Mass |
| --- | --- | --- |
| p0 | (0,0) | 499/1000 |
| p1 | (1/2,0) | 499/1000 |
| v0 | (9/10,0) | 1/1500 |
| v1 | (0,9/10) | 1/1500 |
| v2 | (9/10,9/10) | 1/1500 |

The two residual sources are

    qL=(-4/5,0),  qR=(13/10,0),
    m(qL)=m(qR)=249/500.

Thus m(S)=1 and both widths of S are 9/10. In particular S is eligible.
It is the unique maximum eligible source set, with L=1. To see this without
relying on enumeration, consider any subset containing residual sources.
If it contains both qL and qR, its horizontal span is 21/10 and its mass is
at most the full mass 499/250, so it is not eligible. If it contains just
one residual source and both p0,p1, its horizontal span is at least 13/10
and its mass is at most 749/500, again strictly below its squared span.
If it omits either p0 or p1, its mass is at most

    1 - 499/1000 + 249/500 = 999/1000 < 1.

Every proper subset of S also has mass less than 1. This proves the asserted
unique maximum, and not merely maximality under inclusion.

The full eligible hypergraph is connected: S is one eligible set, and each
pair (p0,qL), (p1,qR) is eligible, since its mass is 997/1000 and its squared
span is 16/25. The full-input nonterminal margin is

    (21/10)^2 - (35/33)(499/250) = 7567/3300 > 0.

The support-tree cost at L=1 is exactly 27/10. In the complete support
graph, the two cheapest edges are the bottom-row consecutive edges, of
costs 2/5 and 1/2. Kruskal then adds two edges of cost 9/10 to connect the
two upper corners. All remaining edges cost at least 9/10; the nontrivial
diagonal costs are 47/50, 19/20 and 99/100. Thus the minimum is 27/10>5/2.
The four corner sources also make the side-1 support-square union the entire
expanded bbox, of area 361/100>3.5. These are failures of sufficient support
bounds, not lower bounds for W-prime.

## 2. Rational outer covers prove that the old credit is zero

Use the following closed-square notation:

    Q(p,r) = p + [-r,r]^2.

Every internal eligible label contains two, one, or no heavy sources p0,p1.
If it contains both, its total mass is at most 1, and its receiver lies in

    E = [0,1/2] x [-1/2,1/2].

If it contains just one heavy source, its mass is at most 501/1000. Since
4(177/500)^2>501/1000, its receiver lies in Q(pi,177/500). If it contains
neither, its mass is at most 1/500. As 4(23/1000)^2>1/500, its receiver
lies in the union of Q(vj,23/1000), j=0,1,2. Therefore the explicit set

    Iplus = E union Q(p0,177/500) union Q(p1,177/500)
              union union_j Q(vj,23/1000)              (3)

contains the complete internal halo I. The argument includes every
internal label and does not assume that only singleton witnesses occur.

The two residual atoms cannot be captured together by a strict residual
witness, since their squared separation exceeds their total mass. Thus R
is exactly the union of their individual open halos. Using

    4(353/1000)^2 > 249/500,

we have

    R subset Rplus = Q(qL,353/1000) union Q(qR,353/1000). (4)

All rectangles in (3) and (4) are disjoint across the two lists. For example,
the left residual cover has right edge -447/1000, while the leftmost
internal cover has left edge -354/1000. The right residual cover has left
edge 947/1000; the internal covers near the bottom right end at most at
923/1000, and the upper tiny-source covers are vertically separated from
both residual squares. Consequently

    I intersect R = empty.                            (5)

In particular the subtraction in (2) cannot repair the mixed increment of
this example.

## 3. Two fixed-witness rectangles give an exact positive patch

Set r=499/1000 and define the open rectangles

    JL=(-499/1000,-301/1000) x (-r,r),
    JR=( 801/1000, 999/1000) x (-r,r).

At every center of JL, the square with half-side r contains p0 and qL;
at every center of JR it contains p1 and qR. Both witnesses are strict:

    4r^2 = 249001/250000 = 0.996004 < 997/1000.

Hence JL union JR is contained in H(A). Define the measurable patch

    F = (JL union JR) \ (Iplus union Rplus).

Then F is contained in P. Boundaries have zero area and do not affect the
following exact calculation. Each J rectangle has width 198/1000 and
height 998/1000. Within each J, the neighboring residual cover removes a
rectangle of width 52/1000 and height 706/1000, and the neighboring heavy
internal cover removes one of width 53/1000 and height 708/1000. The tiny
cover Q(v0,23/1000) lies inside JR and removes one further square of side
46/1000. The two upper tiny covers and E do not meet the J rectangles.
The five removed pieces just listed are pairwise disjoint. Thus

    |F| = 2(198/1000)(998/1000)
             -2(52/1000)(706/1000)
             -2(53/1000)(708/1000)
             -(46/1000)^2
         = 12231/50000
         = 0.24462
         > 1/5.                                      (6)

Combining (5) and (6),

    |P|-|I intersect R| >= 12231/50000 > L/5,
    |D|-|I| >= 12231/50000.

The strict excess is 2231/50000. This is an analytic patch certificate,
not an inferred difference between approximate full-halo areas.

## 4. The new mass-sensitive criterion pays this same input

This example belongs to the R68 support-tree remaining domain, but it does
not survive the R69 heavy-core sufficient criterion. In the notation of
that separate general lemma, choose

    H={p0,p1}, M=499/500, V={v0,v1,v2}, E_mass=1/500,
    beta=0, r_core=43/50.

Both residual sources have distance greater than sqrt(E_mass) from V, so
beta=0. The admissibility condition 4r_core^2<=3M holds. The support
container C=[-1/2,7/5]^2 and the core rectangle

    K=[-9/25,43/50] x [-43/50,43/50]

give

    |C|=361/100,
    |C intersect K|=1037/625,
    |Z_H intersect K|=61/50,
    Gamma=7927/2500.

Here K is not wholly contained in C: its lower vertical side extends
beyond C. The intersection in the area formula is essential. The criterion
then gives

    |W'| <= |W| <= Gamma+9(E_mass+beta)
                      = 1993/625 = 3.1888 < 3.5.

This application uses the separately proved heavy-core lemma; the scalar
check alone would not establish that lemma. The conclusion also makes
clear why refuting the uniform one-fifth increment shortcut does not
refute the original 3.5 goal.

## 5. Actual finite verification and reproduction

The minimal rational input is
`computation/inputs/r069-net-mixed-seven-source.json`. The standalone
`computation/round069/check_net_mixed_obstruction.py` enumerated all 127
nonempty subsets, with 30 closed bboxes and 17 eligible closed bboxes. It
checked the unique maximum, actual eligible connectivity, nonterminal
margin, exact Prim tree cost, every outer-cover separation and fixed-witness
inequality, and the patch area by a separate rational coordinate-cell union
difference. It also evaluated the heavy-core scalar quantities above.

The actual Python 3.15.0 run returned PASS, script wall time 0.035263 seconds,
exit code zero, and the owned child PID was absent after waiting. Source,
input and output hashes and the process receipt are saved in
`computation/results/round069-net-mixed/`. No full-halo evaluation, LP,
parameter scan, random search or MATLAB session was used for this example.

From the integrated topic directory:

```sh
python3 computation/round069/check_net_mixed_obstruction.py \
  --input computation/inputs/r069-net-mixed-seven-source.json \
  --output /Volumes/data/CodexTemp/T2/r069-net-mixed-replay.json
```

Independent review is tracked separately in the R69 round report. The
general W-prime selection, actual deletion selection, general 3.5 bound,
and matching upper/lower determination of the exact C2 remain open.
