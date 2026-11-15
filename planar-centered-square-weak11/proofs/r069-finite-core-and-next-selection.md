# R69: a finite core budget and a refuted mixed-increment alternative

2026-11-15. The finite positive atomic, threshold-one, strict centered-square
setting is fixed throughout. The heavy-core estimates cited below are proved
in `r069-heavy-core-mass-sensitive-payment.md`. This note records their exact
consequence for the remaining selection problem. Section 3 proves that a
conditional core-versus-mixed-increment alternative also fails. The core
theorem itself remains valid; no replacement general dichotomy is asserted.

## 1. The proved core test is finite and independent of halo deletion

Let S be a maximum-mass eligible block, L=m(S), s=sqrt(L), and let
C=bbox(S)+[-s/2,s/2]^2. For each nonempty subset H of S put

    M=m(H), E=L-M, V=S\H,
    B_H={q in A\S:dist_infinity(q,V)<sqrt(E)}, beta_H=m(B_H),
    K_H=intersection_(p in H) [p+[-sqrt(3M)/2,sqrt(3M)/2]^2].

If V is empty, take B_H empty. For t>=0 write Z_t(U) for the union of
side-t open squares centered at U, with Z_0(U) empty. Define

    Phi_S(H) = |C|-|C intersect K_H|+|Z_s(H) intersect K_H|
                   +min{9(E+beta_H), |Z_sqrt(E)(V)|},
    Phi_*(S) = min_(nonempty H subset S) Phi_S(H).       (1)

The main heavy-core theorem and its two refinements prove

    |W'_A(S)| <= |W_A(S)| <= Phi_*(S).                  (2)

The optimal allowed core radius is already sqrt(3M)/2: for a fixed H the
geometric term equals |C|-|(C\Z_s(H)) intersect K_H(r)| and is nonincreasing
in r, while the low-mass term does not depend on r. Thus (1) has at most
2^(rank S)-1 choices and no continuous radius optimization. Its inputs are
source masses, positions, nearby-source membership and finite rectangle
unions; neither the complete residual halo nor D, W or Wprime occurs in
the definition. For rational source data the membership test for B_H is the
strict rational comparison dist_infinity(q,V)^2<E. The remaining root
endpoints admit exact algebraic or directed interval evaluation.
Directed intervals at a fixed precision can remain undecided at equality;
no claim is made that ordinary floating arithmetic or fixed-bit bounds always
decide the finite comparison. Exact algebraic sign handling is a distinct
option when such a boundary must be resolved.

A rational radius below the optimum remains a valid positive certificate.
Conversely, finding one unsuccessful core does not establish Phi_*(S)>3.5L.
All nonempty cores must be controlled to prove that failure condition.

The halfplane-credit class for one residual atom is another proved sufficient
test, independent of whether it is realized by a core in (1). The general
weak-type problem does not require any one particular sufficient test to
work on every input.

## 2. The original uniform mixed allowance is false

Write I=H(S), R=H(A\S), P=H(A)\(I union R), and define

    J_A(S)=|P|-|I intersect R|.

The exact identity

    Delta_A(S)=|I|+J_A(S)                              (3)

does not itself bound either term. R69's seven-source example has a unique
maximum block with two thick axes, rank five, and minimum support-tree cost
27L/10. It satisfies I intersect R empty and

    J_A(S) >= 12231L/50000 > L/5.

Therefore the unconditional estimate J_A(S)<=L/5 is false even in the
R68 support-tree remaining domain. The same example has a two-source core
with a proved bound 1993L/625<3.5L. The core mechanism therefore handles
that particular obstruction. The separate balanced example below shows why
this success does not yield a general core-versus-mixed-increment alternative.

## 3. The conditional alternative is also false

A stronger proposed route was to impose the previous connected/nonterminal,
fat, rank and support-tree-failure domain and then ask:

> If Phi_*(S)>7L/2 for every maximum block S, must at least one maximum
> block satisfy J_A(S)<=L/5?

The answer is **no**. A balanced eight-source input has a unique maximum
block S, Phi_*(S)=361/100>7/2, and J_A(S)>1/5. The following finite formulas
exhibit both failures; the separate R69 balanced-core attachment gives the
full domain and rational patch certificate.

Take the four corners of {0,9/10}^2 as S, each of mass 1/4, and add

    q0=(-3/4,-11/50), q1=(28/25,-3/4),
    q2=(33/20,28/25), q3=(-11/50,33/20),
    m(qj)=749/1000.

These four outside sources are rotations of each other about (9/20,9/20).
The complete eligible catalogue consists of the eight singletons, S itself,
and the four pairs joining each outside source to its closest S corner.
Indeed a subset with one outside source and two S corners has mass at most
1249/1000 for the closest two, but span at least 28/25; adding a third S
corner raises the span to at least 33/20. A subset with two outside sources
has span at least 187/100 and mass at most 2498/1000. Three or four outside
sources force span 12/5, whose square exceeds the full mass 3996/1000.
Far outside/inside pairs are also ineligible, while each closest pair has
mass 999/1000 and span 3/4. Thus S is the unique maximum with L=1. It is
fat, of rank four and tree cost 27/10; the full input is connected and
nonterminal. The full S receiver is the small square (2/5,1/2)^2.

For all fifteen nonempty cores H, let Gamma denote the optimal-radius
geometric term in (1). Put r=sqrt(6)/4 in the two-source cases. Clipping to
C=[-1/2,7/5]^2 gives the exact table

| Core type | Gamma |
| --- | --- |
| one corner | 361/100 |
| two adjacent corners | 241/100+(19/10)r |
| two opposite corners | 59/25+2r |
| three corners | 1419/400 |
| all four corners | 361/100 |

For adjacent corners, one side of K extends beyond C; only one outer strip
of width r-1/2 contributes to the hole. Using an unclipped K would give the
wrong equality. The elementary inequality r>3/5 shows that both two-source
values exceed 1419/400. Thus Gamma>=1419/400 for every core.

For every proper H the leftover mass E is at least 1/4 and V is nonempty.
Consequently both 9(E+beta_H) and |Z_sqrt(E)(V)| are at least E. No near-
outside enumeration is needed for this lower bound. It follows that

    Phi_S(H)>=1419/400+1/4=1519/400>361/100

for every proper H. With H=S the low branch vanishes and Z_s(S)=C, so

    Phi_*(S)=Phi_S(S)=361/100.                         (4)

The complete internal halo consists of four disjoint singleton squares and
the central full-S square, of total area 101/100. The residual halo consists
of four singleton squares. The two halos are disjoint. A fixed half-side
999/2000 gives, at q0 and its nearest corner, an inner receiver rectangle of
width 249/1000 and height 779/1000. Its rotations give four mutually
disjoint rectangles. The residual half-side is less than 433/1000. Inside
each inner rectangle its only possible residual cover has width 365/2000
and height 285/400; the internal halo is disjoint from the entire rectangle.
Thus the rotated rational patch has area

    4[(249/1000)(779/1000)-(365/2000)(285/400)]
       =255759/1000000>1/5.                           (5)

The witness is strict because 4(999/2000)^2=998001/1000000<999/1000.
Equations (4)--(5), together with I intersect R empty, strictly refute the
conditional alternative, including its all-maximum and all-core quantifiers.
There is only one maximum block, and all fifteen cores have been accounted
for analytically; a fixed-bit near-equality issue does not occur here.

This is still not a 3.5 counterexample. The four mixed receiver rectangles
have area strictly below

    (1-3/4)(1-11/50)=39/200.

All W labels are internal or belong to these four pairs. Therefore

    |W'|<=|W|<=|I|+sum_four |E_pair|
                    <101/100+4(39/200)=179/100<7/2.   (6)

Here there are exactly four classified mixed pairs; no sum over an
uncontrolled repeated family is being used. The core envelope loses too
much area on this balanced configuration, while a small explicit label
cover pays it directly.

## 4. What remains after these refutations

The valid finite core theorem is a sufficient structural class. The
unconditional J<=L/5 rule and the conditional assertion that it becomes
valid after every core fails are both frozen as false. Neither should be
left in the next-round state as an unresolved conjecture. Formula (6)
shows that the original Wprime and true-deletion selection problems survive
these refutations. A subsequent joint bound must retain more of the actual
internal-halo slack or of the mixed labels' captured-source geometry than
the uniform one-fifth allowance and the core envelope alone.

R69's eight requested full state evaluations are six historical maximum
blocks from five inputs, plus two fixed mechanism controls. They evaluate
I, R, P, W, Wprime and D, with the residual-free central rectangle K_R and
its complement. They do not minimize (1) over every H, and they do not
include the new balanced eight-source analytic obstruction. Their passing
net bounds are finite evidence only. In particular K_R, defined from the
full maximum block and excluding all residual halo, generally differs from
the core-dependent K_H in (1); when H=S they agree. A pure mixed region may
occur inside K_R even though residual halo and old-halo surplus cannot.

No replacement universal dichotomy is asserted. A general 3.5 selection
bound, if established by another joint payment mechanism, would improve the
current upper estimate but would still not match the known lower bound or
solve the exact C2 objective. The original Wprime selection, true deletion
selection, and matching legal extremizing-sequence problem remain open.
