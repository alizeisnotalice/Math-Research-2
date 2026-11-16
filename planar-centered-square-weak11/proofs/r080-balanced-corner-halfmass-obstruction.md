# R80: an actual balanced-block corner half-mass obstruction

2026-11-16. Centered axis-parallel squares, one common radius, strict halo
threshold 1. All sources in this note are atoms of one actual six-source
input. No internal mass is moved or collapsed.

## 1. Input and claim

Set w=999/1000, e=1/1000 and

    sigma=(delta_(0,0)+delta_(w,0)+delta_(0,w)+delta_(w,w))/4,
    nu=(1/4)delta_(-71/200,0)+(499/1000)delta_(0,-171/200),
    mu=sigma+nu,
    Gamma=[-1/2,-e]², P=(-1,1/2-e)².

The block S consisting of the four old corners is the unique maximum-mass
eligible block, with mass 1. The six-source input satisfies all global closed
square caps, is connected by eligible blocks, and is nonterminal in the
project's span-squared versus (35/33)total-mass criterion. These qualifications
are proved in Section 2 and exactly checked from the raw atoms.

Only the old origin lies in P, so m_P=sigma(P)=1/4. Nevertheless the actual
corner deletion Delta=H(mu)\H(nu), with the complete same residual, satisfies

    |Delta intersect Gamma|
       =244823/500000-(171/200)a-(857/1000)(U-V)
       >129915191/1000000000 >1/8=m_P/2,               (1)

where a=sqrt(1/8), U=sqrt(999/4000), V=sqrt(749/4000).
The strict surplus above 1/8 in the rational certificate is
4915191/1000000000. The exact area is approximately 0.1299165630434473;
the radical formula and rational inequalities, not this decimal, certify it.

This refutes an actual single-corner m_P/2 bound for general balanced
four-corner blocks. This old block is not the prescribed uniform S_n grid
(in particular its side is not the n=2 grid side 3/4). There is no inference
against the uniform-grid single-corner question, a four-corner total payment,
the target total deletion <=7/2, or the general weak constant.

## 2. Exact global qualification from six raw atoms

Use source coordinates in units 1/1000 and masses in units 1/1000:

    index       position       mass
      0          (0,0)           250
      1          (999,0)         250
      2          (0,999)         250
      3          (999,999)       250
      4          (-355,0)        250
      5          (0,-855)        499.

There are three distinct coordinates on each axis. Enumerating the six
closed intervals whose endpoints are those coordinates on each axis gives
36 coordinate boxes and 23 distinct nonempty captured bbox closures. For a
closure T, let M_T be its total mass and s_T its maximum coordinate span.
The exact checks are

    M_T<=max(1,s_T²) for all 23 closures;
    max{M_T:M_T>s_T²}=1,
    the unique closure attaining this eligible maximum is T=S.

Every square's captured set is contained in its source-bbox closure, whose
span is no greater than the square's side. Conversely the source bbox fits
in a square of side s_T, so all-square caps imply the closure inequalities.
Thus the finite table is an exact certificate for all square positions and
sizes, including source-boundary ties, rather than a sample of windows.
For maximum blocks, taking the bbox closure of an eligible subset preserves
eligibility and cannot decrease mass. All masses are positive, so a proper
subset cannot tie the mass of its closure. The unique maximum-closure
certificate therefore proves the claimed unique maximum eligible S among
all source subsets as well.

The old block S has span squared 998001/1000000<1 and is eligible. The
two additional eligible labels {0,4} and {0,5} have respectively

    mass 1/2 > (71/200)²,
    mass 749/1000 > (171/200)².

They connect both outside atoms to the old block, proving the input has one
eligible-block component. The full source span is 927/500 and total mass is
1749/1000. Hence its nonterminal gap is exactly

    (927/500)²-(35/33)(1749/1000)=395579/250000>0.     (2)

The scalar checker saves all 23 masks, bboxes, masses, eligibility decisions
and cap slacks; no numerical tolerance is used for these conclusions.

## 3. Full local menu and the complete residual

The caps imply that no strict full or residual witness has r>=1/2. For a
receiver in Gamma and r<1/2, the open square lies inside P. Its upper
coordinate endpoints are strictly below 1/2-e=.499 and its lower endpoints
strictly above -1. The three other old corners at coordinate w=.999 cannot
be captured. Consequently the actual full corner halo is exactly that of
the three actual atoms at origin, q1, q2, with the stated masses. This is
localization of the original input, not a replacement of its old measure.

Work in positive receiver coordinates X=-x,Y=-y, cropped by
C=[e,1/2]², and put

    d=71/200, D=171/200, v=1/4,
    a=sqrt(1/8), s=sqrt(499/4000),
    V=sqrt(749/4000), U=sqrt(999/4000).

For any nonempty actual label T, its full receiver rectangle is given by
farthest source distance <sqrt(mass(T))/2. All local full masses are below
1, so these radii are below 1/2. Membership gives one common radius between
the farthest distance and the mass radius. Any additional captured atom
only helps; a label need not be the exact captured subset. Conversely the
exact captured subset of a strict witness is one such label. This proves
the full seven-label menu and the complete residual three-label menu.

The four old-containing receivers, in positive coordinates before the crop,
are

    A0=(0,v)²,
    A1=(d-a,a)x(0,a),
    A2=(0,V)x(D-V,V),
    T =(d-U,U)x(D-U,U).

The negative-axis truncations at zero describe their portions in the positive
orthant. The residual receivers there are

    R1=(d-v,d+v)x(0,v),
    R2=(0,s)x(D-s,D+s),
    R12=(d-V,V)x(D-V,V).

All statements about these displayed rectangles are understood after
intersection with C; open rectangle endpoints and closed crop endpoints
have no effect on their Lebesgue areas. Halo thresholds remain strict and
source-boundary mass is not discarded. A zero-width receiver would be empty;
the inequalities below show the receivers used in the area computation
have positive width and height.

The relevant strict orderings are

    e<d-a<v<a<d<V<U<1/2,
    D-s>1/2,  D-U>a,  D-V>a,  D<2V.

In particular R2 misses C, while A2 intersect C=R12 intersect C, since
d-V<0<e. Thus A2 contributes no deletion. The joint T lies above both A0
and A1, intersects neither R1 nor R2, and contains R12 intersect C. The
only residual meeting A0 union A1 is R1. These facts partition the actual
deletion into two pieces with disjoint interiors:

    [(A0 union A1)\R1] intersect C,
    [T\R12] intersect C.                              (3)

Every residual label has been used. In particular the outside pair is
eligible and nonempty; it cannot be omitted or replaced by its singletons.

## 4. Exact area and strict rational lower bound

The first piece of (3) includes a thin genuine old-only strip. Since A1
covers the part of A0 with X>d-a, that extra strip has area
(d-a-e)(v-e). The A1 rectangle itself has cropped area (2a-d)(a-e), and
its overlap with R1 has area (a-d+v)(v-e). Therefore

    area_first=(2a-d)(a-e)-(a-d+v)(v-e)+(d-a-e)(v-e)
              =182323/500000-(171/200)a.              (4)

The last positive strip would be missed by computing only the mixed labels.
The second piece of (3) has area

    area_joint=(U-e)(2U-D)-(V-e)(2V-D)
               =1/8-(857/1000)(U-V).                 (5)

Adding (4)-(5) proves the equality in (1). All overlaps are subtracted using
the same residual: R1 once in the first piece and R12 once in the second.

The following elementary directed bounds suffice for both the orderings
and the strict lower bound:

    176776/500000<a<176777/500000,
    499749/1000000<U<1999/4000,
    432723/1000000<V<432724/1000000.

Each is certified by squaring positive rationals against the corresponding
radicand. The further comparison s<D-1/2=d follows from
s²=499/4000<d². These inequalities imply every ordering in Section 3.
In particular U+a<.49975+.353554<.855=D. Using the upper bounds on a,U
and lower bound on V in (4)-(5) gives

    area_first>6235733/100000000,
    area_joint>67557861/1000000000,
    area_first+area_joint>129915191/1000000000>1/8.

The calculation concerns the original full halo and full residual on the
original Gamma. It needs neither a comparison to a synthetic grid proxy nor
an atomic approximation to a continuous outside density.

## 5. Actual checks and scope

The minimal raw input is
`computation/inputs/r080-kernel/balanced-six-source.json` and the exact checker
is `computation/round080-kernel/check_balanced_corner.py`. The checker
enumerates all coordinate boxes and all subset-derived bbox closures from
the six raw atoms, verifies their agreement and the complete exact cap,
maximum and qualification table, then checks the radical bounds, geometric
orderings, algebraic rational constants and strict lower bound. Its result
and actual Python/version/input/script/timing receipt are saved at
`computation/results/round080-kernel/checks.json`.

The area identity is the analytic rectangle partition (3)-(5). This author's
small check does not rerun a rectangle-union area engine, compute the other
three corners or total deletion, or expand to another input. Independent
full-area computation and review, if supplied, retain their own source-bound
receipts. R79 and earlier frozen files remain unchanged; this component uses
no MATLAB, new agent, commit or push.
