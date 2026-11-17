# Complete set implies constant width: self-contained audit

Status: **verified after one necessary correction** to the proposed normal-cone
step.  The proof below removes the need to import the planar
`complete => constant width` theorem.  It uses only compactness and the
finite-dimensional strong separation theorem.  The older audit is deliberately
unchanged; this is an independently reviewable replacement lemma.

The conclusion agrees with Vrecica, *A Note on Sets of Constant Width*, Publ.
Inst. Math. 29 (1981), 289--291: its p. 289 says that in Euclidean spaces
complete sets and constant-width sets coincide, and its Theorem 1 gives a
constant-width completion.  The proof here does not invoke that assertion.

## Statement

Let `K subset R^2` be inclusion-maximal among all sets of diameter at most one.
Then `K` is a compact convex body of constant width one.  Consequently every
nonempty bounded `S` with `diam(S)<=1` is contained in such a body: use Zorn's
lemma on the family of diameter-at-most-one supersets of `S`, since a chain
union still has diameter at most one.

## Preliminary facts supplied by maximality

For a maximal `K`,

`K = intersection_{b in K} closed_ball(b,1)`.  (The right side can be adjoined
to `K` without increasing its diameter.)  It follows that `K` is closed and
convex; it is bounded because it lies in `closed_ball(b_0,1)` for any
`b_0 in K`, hence compact.

It has nonempty interior.  Otherwise the compact convex set `K` lies in a line,
so is a point or a segment of length at most one.  A sufficiently small nonzero
perpendicular displacement of its midpoint has distance strictly below one from
every point of that segment, contradicting maximality.

Choose `z in int(K)`.  Then

`rho := max_{b in K}|z-b| < 1`.

Indeed equality for some `b` would allow a sufficiently short displacement of
the interior point `z` directly away from `b`, still within `K`, at distance
greater than one from `b`.  Compactness gives the displayed uniform strict
inequality.

## The active-ball normal lemma

Let `a in bd(K)` and let `u` be a unit outer support normal at `a`, i.e.

`u dot (x-a) <= 0` for all `x in K`.

The active centres

`A(a) = {b in K : |a-b|=1}`

form a nonempty compact set.  Nonemptiness follows because otherwise all ball
constraints are uniformly strict at `a`, putting `a` in the interior of their
intersection, which is `K`.  Put

`V(a)={a-b:b in A(a)}` and `C=closure(pos V(a))`, where `pos` means finite
nonnegative linear combinations.

**Claim.** `u in C`.

If not, finite-dimensional strong separation of the point `u` from the closed
convex cone `C` gives `w` such that

`u dot w>0` and `v dot w<=0` for every `v in C`.

Put `q=z-a`.  For every active `v=a-b`, the strict interior point estimate
gives

`1+2 v dot q+|q|^2 = |z-b|^2 <= rho^2`,

so `v dot q <= -(1+|q|^2-rho^2)/2 =: -delta <0`.  Choose a sufficiently small
`eta>0` so that `u dot (w+eta q)>0`, and write `y=w+eta q`.  On active centres,
`(a-b) dot y <= -eta delta`.

More explicitly, by compactness and continuity there is a relatively open
neighbourhood `U` of `A(a)` in `K` on which `(a-b)dot y<=-gamma<0`.  For
`b in U`,

`|a+t y-b|^2-1 <= -2 gamma t+t^2|y|^2 <=0`

when `t>0` is sufficiently small.  On the complementary compact set `K minus
U`, the constraint `|a-b|^2-1` has maximum `-kappa<0`; its perturbation is
bounded by `2t|y|+t^2|y|^2`.  Reducing the same positive `t` if necessary makes
this at most `kappa`.  Thus

`|a+t y-b|<=1` for every `b in K`.

The ball-intersection identity puts `a+t y` in `K`, whereas
`u dot(a+t y)>u dot a`, contradicting that `u` supports `K` at `a`.  This
proves the claim.

This is the point where the proposed argument needs care: an infinite-ball
normal-cone rule should first give membership in the **closed** positive hull.
The preceding separation proof establishes exactly that, with a justified
uniform Slater point `z`.  It does not prematurely assume an exact finite
multiplier formula.

## Turning the closed-cone statement into the desired antipode

Choose nonzero finite combinations

`c_n = sum_i lambda_{n,i} v_{n,i} in pos V(a)`, `lambda_{n,i}>=0`,

with `c_n -> u`, and normalize `w_n=c_n/|c_n|`.  Thus `w_n -> u` and

`w_n=sum_i mu_{n,i}v_{n,i}`, `mu_{n,i}>=0`, `L_n:=sum_i mu_{n,i}>=1`,

because every `v_{n,i}` has length one.  For every `x in K` and active centre
`b`, expansion of `|x-b|^2<=|a-b|^2=1` yields

`|x-a|^2+2(x-a) dot (a-b)<=0`.

Multiply by the `mu_{n,i}` and sum.  Since `L_n>=1`,

`|x-a|^2+2(x-a) dot w_n<=0`.

Passing to the limit gives

`|x-a|^2+2(x-a) dot u<=0` for every `x in K`.

Set `d=a-u`.  The last inequality says

`|x-d|^2=|x-a|^2+2(x-a) dot u+1<=1` for all `x in K`.

The ball-intersection identity therefore gives `d in K`.  Since `u dot d=u dot
a-1` while the diameter bound gives `u dot a-u dot x<=1` for every `x in K`,
the width of `K` in direction `u` is exactly one.  Every unit direction is an
outer normal at some support point of compact convex `K`; therefore all widths
equal one.

## Audit of the proposed proof

The proposed inequality calculation is correct, but it suppressed two points:

1. From `u=sum lambda_i(a-b_i)`, one only gets
   `L|x-a|^2+2(x-a)dot u<=0`, where `L=sum lambda_i`.  The desired version
   follows because `L>=|u|=1`; it is not literal termwise cancellation.
2. For infinitely many balls, directly asserting an exact finite multiplier
   formula needs a normal-cone theorem plus a constraint qualification.  The
   proof above supplies the qualification using `z in int(K)`, proves the
   weaker but sufficient closed-cone statement by separation, and then uses
   normalized finite approximants.  Thus no unproved infinite-intersection rule
   remains.

The argument is Euclidean: the squared-norm expansion and the unit-vector
normalization are essential.  It should not be exported unchanged to a general
Minkowski normed plane.
