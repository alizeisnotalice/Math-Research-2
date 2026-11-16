# R71: shared outside mass, exact residual exclusion, and geometric dual certificates

2026-11-16, Asia/Shanghai. Author proof. Independent mathematical review
and final source binding are recorded separately.

## 1. Setting and evidence boundary

Fix distinct internal points S with positive masses sigma_p, total mass 1,
and coordinate spans strictly below 1. Squares are centered, axis parallel,
with one common radius and strict threshold 1. For a finite positive measure
mu write

    H(mu) = {x: mu(Q_open(x,r)) > 4r^2 for some r>0}.

The object below is actual deletion H(sigma+nu) minus H(nu), not the
R70 source-only outer envelope. It retains the same outside measure nu for
all receivers. No area bound for deletion or new bound for C2 is proved.

For finite sources, an eligible block has squared maximum coordinate span
strictly less than its mass. S is required to remain a maximum eligible
block of mass 1. We permit zero masses at candidate outside positions;
zero positions are ignored when interpreting the actual positive support.

The inherited exact closure-cap equivalence is restated and proved below
so boundary conventions are explicit. Source: R64, Section 2 of
`r064-source-state-mass-polytope-variation.md`; final inherited hash is
recorded in this round's dependency list. R70 has already excluded any
uniform area bound below 4 for its source-only all-prefix envelope. No
such budget is proposed again here.

## 2. All bounding-box caps are exact, including zero positions

Fix a finite list Q={q_1,...,q_N} of distinct positions disjoint from S,
and put nu_y=sum_l y_l delta_(q_l), y_l>=0. For every nonempty complete
bounding-box closure B of the fixed position list S union Q, let d_B be
the maximum span of its tight bounding box. Impose

    sum_(q_l in B) y_l <= b_B,
    b_B=max(1,d_B^2)-sigma(B).                            (1)

These finitely many linear inequalities are exactly the condition that
S remains a maximum eligible block of mass 1. To see necessity, any
violating closure has mass greater than 1 and greater than d_B^2, so is
an eligible heavier block after zero positions are omitted. Conversely,
a positive eligible subset of mass greater than 1 can be completed to
its tight bounding-box closure: its span does not change and its mass
can only increase, violating (1). S itself is eligible and has mass 1.
If a closure contains zero-mass positions at extremes, removing them
only decreases its actual span; this does not spoil the first argument.

Every b_B is nonnegative because sigma(B)<=1. The zero vector y=0 is
feasible. Singleton closures at the distinct outside positions give
0<=y_l<=1, so the feasible mass set is compact. Positions in the closed
maximal-block collar are forced to zero when a side-1 square contains
them together with S. Keeping those columns is harmless; pretending
they are positive support is not.

## 3. Residual exclusion is a complete closed-prefix condition

For a receiver x, sort the distinct distances ||q_l-x||_infinity as
0<=t_1<...<t_M. Every tie is included in its complete prefix

    R_(x,t)={q_l: ||q_l-x||_infinity<=t}.

Then

    x not in H(nu_y)
    iff sum_(q_l in R_(x,t)) y_l <=4t^2 at every t=t_i.   (2)

For necessity, if the closed-prefix mass were greater than 4t^2, choose
r>t sufficiently close that this same mass is greater than 4r^2. Its
points are then strictly inside Q_open(x,r), a residual witness. For
sufficiency, an open square at radius r captures a closed complete prefix
at its largest captured distance t<r, or no positive mass. Its mass is
at most 4t^2<=4r^2 by (2). Extra zero-mass positions and distance ties do
not change either argument. A distance-zero row is essential: it forces
all residual mass located exactly at x to vanish. Equalities in (2) are
allowed and must not be converted to strict inequalities.

Thus residual exclusion at finitely many receivers is a finite family
of closed linear halfspaces, not a sampled approximation in radius.
The largest-distance row is retained; large residual radii cannot be
ignored merely because the desired mixed witness has radius below 1/2.

## 4. One common mass vector and a strict-margin linear program

Choose receivers x_1,...,x_k, k>=1, and fixed open witness radii
0<r_j<1/2. Put

    W_j=Q_open(x_j,r_j),
    d_j=4r_j^2-sigma(W_j),
    B_(j,l)=1_(q_l in W_j).

Let A y<=b be (1), and R y<=c all residual rows (2) for all k receivers.
The exact fixed-menu feasibility problem is

    y>=0, A y<=b, R y<=c, B y>d,                         (3)

with componentwise strict inequality. It says exactly that one actual
outside input makes every prescribed radius a strict full witness,
while every receiver remains outside the complete residual halo and S
remains maximum. If x_j is also outside H(sigma), that receiver belongs
to the purely mixed part; this additional test depends only on S.

Maximize a single scalar epsilon over

    y>=0, A y<=b, R y<=c, B y>=d+epsilon*1.              (4)

Here epsilon is unrestricted in sign. Its maximum epsilon_* is finite
and attained: the allowed y-set is nonempty compact and the objective
is equivalently max_y min_j((By)_j-d_j), a continuous function. Then
(3) is feasible if and only if epsilon_*>0. The minimum of finitely many
strict positive witness slacks supplies one common positive epsilon.
An optimum equal to zero is strict infeasibility, not a nearly feasible
strict witness. Bounding epsilon below by zero before solving is not
necessary and can obscure the distinction between negative and zero
optima. Rational data give a rational LP and rational exact certificates.

The weak dual is explicit. Let alpha,beta,lambda be componentwise
nonnegative and satisfy

    sum_j lambda_j=1,
    A^T alpha + R^T beta >= B^T lambda.                  (5)

Then every feasible point of (4) obeys

    epsilon <= b.alpha+c.beta-d.lambda.                 (6)

Indeed, multiply the witness rows by lambda, the upper rows by alpha
and beta, and use y>=0 and (5). Thus a right side <=0 proves strict
infeasibility, including a zero-gap certificate. A right side <0 gives
a strict negative upper bound. Finite-dimensional LP duality also
provides such an attained dual optimum here; no infinite-dimensional
strong duality is invoked. Every solver-proposed certificate must be
rationalized and (5)-(6) checked exactly. A floating infeasibility status
or small residual is not a certificate.

## 5. Exact finite radius menus for all actual witnesses

A selected radius menu is not the whole deletion problem. For finite
S union Q and each x_j, list all distinct distances t of these fixed
positions from x_j with 0<=t<1/2, and use the CLOSED witness square

    W_(j,t)=Q_closed(x_j,t),
    d_(j,t)=4t^2-sigma(W_(j,t)).                          (7)

For each j choose one such t. Replace the witness rows in (3)-(6) by
these closed-square rows. Strict positivity of the same common slack
is still exact: mass greater than 4t^2 allows a slightly larger positive
open radius r>t with 4r^2 below that mass. This also handles t=0.

Conversely, every actual strict witness has captured mass at most 1,
because its captured positive sources form an eligible block under (1).
Hence its radius is less than 1/2. Its largest positive captured source
distance t is in the finite list and satisfies t<r<1/2; the complete
closed prefix at t has mass greater than 4t^2. Zero-mass extra positions
in that prefix cause no problem.

Therefore all k receivers can belong simultaneously to actual deletion
for some y on Q if and only if AT LEAST ONE finite choice of one menu
entry per receiver has epsilon_*>0. Infeasibility of one fixed menu is
not infeasibility over all radii. A certificate for every menu, or a
separate argument valid for every radius, is needed for that conclusion.
The residual constraints remain all complete residual distance prefixes
in every menu, and every menu still uses one shared y across receivers.

## 6. A finite certificate that applies to arbitrary outside positions

Let D=R^2\S be the allowed outside-position domain, and let nu be any
finite positive Borel measure supported on D. Require the measure version
of the same maximum-block constraints:

    sigma(B)+nu(B)<=max(1,d_B^2)                         (8)

for every closed bounded axis-parallel rectangle B, where d_B is its
maximum side length. Equivalently one may require (8) only for closed
axis-parallel squares, since each rectangle fits in a square of that
side length. For finite atomic inputs this is exactly (1), now without
fixing candidate positions. S has mass 1 and spans below 1 as before.
This section does not assign an unproved atomic-block definition of L
to general diffuse measures; (8) is the explicit admissibility condition.

There is a useful finite zero-cost localization. Write bbox(S) as
[a,b] times [c,d]. Its closed collar is

    J=[b-1,a+1] times [d-1,c+1].

The four closed side-1 squares obtained by choosing the x interval from
[b-1,b] or [a,a+1], and the y interval from [d-1,d] or [c,c+1], each
contain S. Their sigma-mass is 1, so (8) gives nu-mass zero for each.
In each axis the two intervals overlap because the span is below 1,
and their union is the displayed collar interval. The four squares
therefore cover the entire CLOSED J, including all boundary points,
and nu(J)=0. One may restrict the full admissible measure problem to
outside positions in R^2\J. In a finite selected-row relaxation this
restriction requires retaining these four zero-cost rows or explicitly
restricting its allowed domain; it is not inferred from omitted caps.

For arbitrary finite measures the exact residual condition is

    nu(Q_closed(x,t))<=4t^2 for every t>=0.               (9)

The necessity proof from Section 3 is unchanged. For sufficiency,
Q_open(x,r) is the increasing union of closed squares of radii t<r,
so continuity from below gives nu(Q_open(x,r))<=4r^2. Thus (9) is exactly
x not in H(nu), including atoms on boundaries. No finite radius menu
is asserted in the diffuse case. There is nevertheless an exact
countable witness reduction: (8) excludes every strict full witness at
radius r>=1/2, and any smaller strict open witness can have its radius
slightly increased to a rational r'<1/2 while preserving its mass excess.
For finitely many receivers, arbitrary witness choices are therefore
the countable union of rational-radius menus. Checking finitely many of
those menus does not exhaust this countable union.

Fix any finite witness menu, using open positive-radius squares or the
closed squares (7). Let d_j=4r_j^2-sigma(W_j), with r_j=t_j in the closed
case. Take finitely many rectangles B_a from (8), residual squares V_b
from (9) at tested receivers, and nonnegative weights alpha_a,beta_b,
lambda_j with sum lambda_j>0. Suppose the POINTWISE inequality

    sum_j lambda_j 1_(W_j)(q)
      <=sum_a alpha_a 1_(B_a)(q)+sum_b beta_b 1_(V_b)(q) (10)

holds at every possible outside position q. If

    sum_a alpha_a[max(1,d_(B_a)^2)-sigma(B_a)]
      +sum_b beta_b*4(radius(V_b))^2
      <=sum_j lambda_j*d_j,                            (11)

then no such nu can make all witnesses strict while retaining all
residual exclusions. Integrating (10) against ONE shared measure nu,
using (8)-(9), and summing the strict witness inequalities gives a
strict lower bound greater than the right side of (11), a contradiction.
This is only weak duality and is valid without any existence or
attainment theorem in an infinite-dimensional measure space.
If (10) has been checked only outside J, it can be completed globally
at no cost: add each of the four collar squares as an upper row with
weight Lambda=sum_j lambda_j. At every point of J at least one such
row contributes Lambda, while the witness multiplicity is at most
Lambda. All four row costs are zero. This supplies an exact certificate
rather than silently discarding boundary columns inside the collar.

The support condition for (10) must be explicit. Checking it only at a
finite candidate set Q proves only the fixed-Q statement. A dual found
on Q that fails at a new location cannot be promoted to arbitrary nu.
By deleting the part of nu outside the union K of the selected witness
squares, one preserves every witness, decreases all upper-constraint
masses, and preserves residual exclusion. Consequently one may restrict
the search to nu supported on K. For the nonnegative weighted upper
rows in (10), outside K the left side is zero anyway.

## 7. Complete geometric verification and exact finite row relaxations

For a proposed finite certificate (10), collect every finite x and y
boundary coordinate of its witness and upper rectangles. In each axis
partition the line into singleton boundary coordinates, the open gaps
between them, and the two open unbounded intervals. Their Cartesian
products are finitely many disjoint atoms. Every open or closed
rectangle indicator in (10) is constant on every such atom.

Thus checking one representative per atom proves (10) everywhere,
including all edges and vertices. For rational boundaries, rational
representatives suffice. Ignoring singleton-coordinate atoms can lose a
constraint: open witness and closed cap conventions must be retained.
An atom with negative upper-minus-witness slack is an exact missing
column for a fixed-Q certificate. This detects a genuine geometric gap
in the certificate, not a floating solver error. No claim is made that
repeated column and constraint generation must terminate for the full
infinite rectangle/radius families.

There is also an exact finite reduction for a FINITE SELECTED FAMILY
of upper rows and witness rows. Intersect the arrangement atoms with K
and the allowed domain D=R^2\S, discard empty intersections, and choose
one representative from each. Removing the finitely many internal source
positions is essential when exact outside support is required. With
rational data representatives can still be chosen rationally.
Any finite nu supported on K can be replaced by the atomic measure whose
mass at a representative is the nu-mass of that atom. Every selected
rectangle mass is preserved exactly. Conversely any nonnegative atom
mass vector is a finite positive measure. Therefore the finite-row
measure relaxation is exactly the LP on all these representatives.
This does not preserve unselected bounding-box caps or unselected
residual radii. Its infeasibility is a valid obstruction for the full
problem; its feasibility is only feasibility of that selected-row
relaxation. For a newly feasible atomic candidate, all exact closure caps
and all residual distance prefixes must be rechecked.

## 8. A direct union-and-overlap consequence

The same common-measure mechanism can be written without LP notation.
For two selected witness squares W_1,W_2 and positive demands d_1,d_2,
choose a closed rectangle B containing W_1 union W_2, and any selected
residual square V containing W_1 intersect W_2. Pointwise,

    1_(W_1)+1_(W_2)<=1_B+1_V.

Consequently simultaneous strict witnesses and residual exclusion are
impossible whenever

    d_1+d_2>=max(1,d_B^2)-sigma(B)+4(radius(V))^2.        (12)

If the witness squares are disjoint, omit V and its cost. More generally
any nonnegative rectangle cover of the multiplicity function
sum lambda_j 1_(W_j) supplies a valid joint-demand budget by (10)-(11).
It is the union and actual overlap of the exterior captured source
regions that are paid; the same mass cannot be chosen afresh at each
receiver. This geometric inequality is general, but no claim is made
that it pays enough receivers to imply a 3.5 deletion-area theorem.

## 9. From shared point witnesses to certified receiver area

There is an exact residual exclusion interface for a whole receiver
rectangle, which makes the next area obligation concrete. Let
C=[u,v] times [s,t] be a nonempty closed rectangle. For a nonempty
residual source subset T with bbox [a,b] times [c,d], define

    delta_C(T)=max{(b-a)/2,(d-c)/2,b-v,u-a,d-t,s-c,0}.     (13)

Then delta_C(T)=min_(x in C) max_(p in T)||x-p||_infinity.
In one coordinate, the smallest possible maximum distance from [a,b]
while the center lies in [u,v] is max{(b-a)/2,b-v,u-a}; this follows
by intersecting [b-r,a+r] with [u,v]. The two center coordinates are
independent, so their joint optimum is the maximum of the two values.
The minimum is attained because C is compact. Thus, for fixed residual
positions, exactly

    C intersect H(nu_y) is empty
    iff nu_y(T)<=4 delta_C(T)^2 for every residual T.    (14)

It suffices to use complete residual bounding-box closures: completing
a subset leaves its bounding box and delta unchanged and increases its
mass. Zero-mass positions are handled exactly as in Section 2. These
are linear mass inequalities with constants determined solely by C and
the fixed positions, retaining every residual radius and center in C.

For any chosen full source label U put

    R_C(U)=max_(x in C,p in U)||x-p||_infinity.

If (14) holds and (sigma+nu_y)(U)>4 R_C(U)^2, then every x in C is a strict full
receiver and is outside the residual halo. Consequently C is contained
in actual deletion. The last inequality is another strict linear mass
row. Distinct receiver rectangles can use different labels while all
rows share one outside mass vector. A finite disjoint-interior family
therefore certifies a deletion area at least the sum of its areas.

This full-label condition is sufficient for an entire rectangle, not
necessary for coverage by a union of several full-label receivers.
Likewise feasible cells provide an area LOWER certificate for an actual
input; they do not yield a general area upper bound. At a point where a
residual prefix is exactly tight, a two-sided residual-free neighborhood
need not exist even when the full witness has positive slack. One-sided
cells or exact residual-rectangle subtraction must handle that boundary.

An area-oriented continuation can therefore evaluate complete true
full-minus-residual halo geometry for a certified positive input and
optimize jointly realizable receiver cells, with (14) providing exact
residual exclusion. To approach a universal upper budget instead,
infeasibility cuts must control every possible receiver cell and radius
choice and must remain valid under arbitrary outside positions; a finite
selection of positive cells cannot replace that obligation.

## 10. What the computational pilot must establish

A fixed-Q, fixed-menu obstruction is meaningful only with exact source
positions and masses, complete closure and residual-prefix row lists,
a rational dual verified by (5)-(6), and a statement of which individual
receivers are feasible. It becomes an arbitrary-position obstruction
only after (10) is checked on the full arrangement. Alternatively a
failed arrangement atom supplies a new position to test. All-radius
claims require the finite menu union of Section 5 in the fixed-Q case,
or another uniform proof for general measures.

This interface preserves actual shared source geometry and exact
residual-halo subtraction. It does not yet provide a general area
budget, a maximizing-block selection theorem, or matched C2 bounds.
