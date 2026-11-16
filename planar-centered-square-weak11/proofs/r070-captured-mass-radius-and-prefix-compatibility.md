# R70: captured mass, receiver radius, and all-core compatibility

## 1. Finite setting and status

Let A be a finite set of distinct points of R² with positive masses. Squares
are centered, axis parallel, and tested at the strict threshold 1. Write
m(T) for the mass of a nonempty subset and d(T) for the larger of its two
bounding-box spans. A subset is eligible when d(T)² < m(T). Let

    L = max{m(T): T eligible},   s = sqrt(L),

and fix an eligible S of mass L. Thus both spans of S are strictly below s.
For any nonempty T and positive v define

    E_T(v) = {x: max_{p in T} ||x-p||_infinity < sqrt(v)/2}.

This is the open receiver rectangle, possibly empty. The halo of an input
is the union of E_T(m(T)) over its nonempty subsets. Let I=H(S),
R=H(A\S), and D=H(A)\R. Use the R67 definitions of W and W′: their
individual outer rectangles are E_T(m(T)) for eligible, bounding-box-closed
T meeting S; their cores use respectively the same full bounding box and
the actual residual subset. In particular D ⊂ W′ ⊂ W.

The results below are **necessary receiver constraints**, valid for every
finite A and every maximizing S. The scalar elimination in Section 3 is
exact for the stated relaxed scalar problem. It does not assert the
existence of a common actual outside input. Neither the original nor the
all-prefix outer region is asserted to have area at most 7L/2.

## 2. Captured mass forces distance from S

Take any mixed eligible T, put U=T∩S, a=m(U)>0, b=m(T\S)>0,
and m=a+b. Maximality gives m≤L, so 0<a<L. For x∈E_T(m), write

    c_S(x) = max_{p in S} ||x-p||_infinity.

Then

    c_S(x) > sqrt(L+b) - sqrt(a+b)/2.                    (1)

Indeed S∪T has mass L+b>L. If (1) failed, every cross distance between a
point of S and a point of T would be strictly less than sqrt(L+b), because
the T point is strictly within sqrt(m)/2 of x. The within-S spans are below
sqrt(L), and the within-T spans are below sqrt(m); both are also below
sqrt(L+b). Thus S∪T would be eligible with mass greater than L, a
contradiction. T need not be closed, and no eligibility of a smaller
internal subset is assumed. The strict receiver inequality is what makes
the contradiction valid even at equality in the proposed reverse of (1).

For fixed a, set f_a(m)=sqrt(L+m-a)-sqrt(m)/2 on a<m≤L. Its derivative
has the sign of 3m-(L-a). Consequently its infimum is

    rho_L(a) = sqrt(3(L-a))/2,       0<a≤L/4;
               sqrt(L)-sqrt(a)/2,  L/4≤a<L.             (2)

For a<L/4 the minimum is attained at m=(L-a)/3, equivalently
b=(L-4a)/3>0. At a=L/4 and above, the endpoint b=0 is excluded; the
displayed value is an infimum only. In every case (1) implies the strict
inequality c_S(x)>rho_L(a). The two formulas agree at a=L/4.

For comparison, with z=c_S(x)/s this forces

* z>1/2;
* a/L>4(1-z)² when 1/2<z≤3/4;
* a/L>1-4z²/3 when 3/4<z<sqrt(3)/2.

At z≥sqrt(3)/2 this inversion supplies no additional positive lower bound
beyond a>0. These are still necessary conditions, not a prescription for
placing outside mass.

## 3. Exact scalar elimination retaining the internal receiver radius

For a fixed nonempty proper U⊂S put a=m(U) and
r_U(x)=max_{p in U}||x-p||_infinity. Define

    q_U(x) = max{a, 4r_U(x)², (L-a)/3}.                 (3)

There exists a real scalar m satisfying

    a<m≤L,    4r_U(x)²<m,    c_S(x)>f_a(m)              (4)

if and only if

    q_U(x)<L,    c_S(x)>f_a(q_U(x)).                    (5)

To prove this, f_a decreases up to (L-a)/3 and increases thereafter.
The infimum on the interval in (4) is therefore f_a(q_U). If q_U is an
excluded lower endpoint, the strict final inequality in (5) permits a
sufficiently small positive increase of m. If q_U is the interior critical
point, choose that point. Conversely (4) bounds c_S strictly above this
infimum. The condition q_U<L is essential: q_U=L could only arise from
4r_U²≥L, when the strict radius requirement leaves no feasible m. It is
not a permissible equality case in the exact set definition.

Let M(S) be the union of (5) over all nonempty proper subsets U that are
closed relative to S under their own bounding box. Then

    D ⊂ W′ ⊂ W ⊂ I ∪ M(S).                            (6)

For a pure internal label its receiver is contained in I. For a mixed
closed label T, the subset U=T∩S is nonempty, proper, and closed relative
to S; any S point in bbox(U) also lies in bbox(T). Its actual m=m(T)
satisfies (4) by (1). This proves the last inclusion for every outer
rectangle used by W, and hence proves (6).

If x∉I, every U satisfies a≤4r_U(x)². Thus the first entry a in (3) can
be removed for this part of the plane. This simplification does not remove
the spatial radius r_U or the shared-outside-input obligation.

For rational a,L,c,r and q as in (3), let

    z = L + 3q/4 - a - c².

The strict inequality c>f_a(q) is exactly

    z<0,  or  [z≥0 and c²q>z²].                        (7)

Indeed it is equivalent to c sqrt(q)>z; c≥0 and q>0 justify the case
split before squaring. This permits exact rational directed-cell
classification. Boundary cells must remain unresolved or be included
only on the outer side; finite precision is not a general equality
decision procedure. R70 production saves all unresolved leaves.

## 4. Compatibility with every uncaptured internal core

Keep a mixed T, U=T∩S and m=m(T). For every nonempty
H⊂S\U with M=m(H) and M+m>L, one has

    c_H(x) > sqrt(M+m)-sqrt(m)/2,
    c_H(x)=max_{p in H}||x-p||_infinity.                 (8)

The proof is the same union argument, now for H∪T, whose mass is M+m.
The within-H spans are below s because H⊂S; they are strictly below
sqrt(M+m) under the hypothesis. The within-T spans and all cross spans
also are strictly below sqrt(M+m) if the reverse of (8) holds. Thus H∪T
would be an eligible subset of mass greater than L. No eligibility of H
is required. The trigger is strictly M+m>L: equality does not impose
this constraint by maximality alone.

At a fixed receiver x, list the distinct distances of S\U from x as
R_1<...<R_k, and let

    H_j={p in S\U: ||x-p||_infinity≤R_j},  M_j=m(H_j).

All sources at a tied distance are included. Then the entire collection
of constraints (8) is equivalent to the at most |S\U| prefix tests

    M_j+m≤L,  or
    R_j>sqrt(M_j+m)-sqrt(m)/2.                          (9)

One direction follows by taking H=H_j. Conversely, given H, take the
prefix at R_j=c_H(x). It has the same maximum distance and at least the
mass of H. If H triggers, the prefix triggers, and the right side of (8)
is increasing in M. Thus the prefix inequality implies the one for H.
This also proves why tied distances must be taken together.

Writing r=sqrt(m)/2, the active branch in (9) is equivalently

    M_j < (R_j+r)²-4r²
        = (R_j-r)(R_j+3r).                             (10)

This equivalence follows by squaring the nonnegative quantities
R_j+r and sqrt(M_j+m). It retains a strict inequality. In particular,
an active prefix with R_j≤r is impossible. The alternative M_j≤L-m
must not be dropped from (9).

Define M_prefix(S) by requiring the existence of one common scalar m
for (4) and all tests (9), for at least one proper closed U. Then

    W ⊂ I ∪ M_prefix(S) ⊂ I ∪ M(S).                    (11)

The same actual T,m verifies all tests simultaneously. Optimizing m
separately for different prefixes would give a further, unjustified
relaxation of this definition. The finite-prefix theorem reduces the
number of tests at fixed x,U,m; it is not an area estimate and does not
establish realizability by an actual common outside input.

## 5. What is paid, and what remains

The separate R70 balanced-four-corner theorem proves a positive class:
for equal quarter masses at the corners of a square of side t sqrt(L),
sqrt(3)/2≤t≤9/10, every outside input retaining S as a maximizing block
satisfies |W|≤(32+sqrt(2))L/10<3.35L. That statement uses the necessary
captured-mass constraints with the particular internal geometry; it is
not inferred from finite experiments.

The companion `r070-four-corner-radius-envelope-obstruction.md` proves
that |I∪M(S)|≤7L/2 fails in a legal unique-max, fat, connected,
nonterminal input beyond the support-tree criterion. More strongly,
`r070-grid-clamped-envelope-near-four.md` proves, after normalization,
that both |I∪M(S_n)| and |I∪M_prefix(S_n)| tend to 4 for uniform finite
grids in the same legal domain. The latter result imposes every prefix
for the SAME scalar m_n at each receiver; it is not a separate
minimization over m for each prefix. Its uniform discrepancy proof
controls all active prefixes simultaneously. Hence even this strengthened
source-only envelope cannot have a universal area bound below 4L. The
three production pilot outputs evaluate M(S), not M_prefix(S); their
finite results are distinct from that analytic limit theorem.

The general selection problem for |W′| or |D|≤7L/2 remains unresolved,
as does the original matching determination of C₂. Actual outside source
positions, their common mass budget across receivers, and residual halo
payment remain missing from the virtual scalar conditions.
