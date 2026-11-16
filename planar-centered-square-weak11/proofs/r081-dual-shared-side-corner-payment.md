# R81: shared-side capacity pays a fixed four-corner family

2026-11-16. Centered axis-parallel common-radius squares, strict threshold
1, finite positive atomic sources. This note gives a new joint payment for
one fixed twelve-source one-parameter family. It does not claim an arbitrary
maximum-block theorem. Frozen R79 and R80 files are unchanged.

## 1. General bookkeeping and the inherited boundary

Normalize the maximum eligible mass to 1, with an old maximizing block
sigma of mass 1 and bbox [0,w]x[0,h], w,h<1. The maximum-block property gives
all closed square caps mu(Q)<=max(1,side(Q)²), hence every strict full and
residual witness has radius r<1/2. Write D=H(sigma+nu)\H(nu), with one actual
outside nu throughout. As in R63, D lies inside

    C=(-1/2,w+1/2)x(-1/2,h+1/2).

For w,h>1/2, let the four corner boxes be the products of the two horizontal
intervals (-1/2,w-1),(1,w+1/2) and corresponding vertical intervals.
Their total area is 4(w-1/2)(h-1/2); their complement N in C has area

    |N|=3(w+h-wh).

Thus any joint corner certificate with

    sum_corners |D intersect Gamma| - |N\D|
       <=1/2+3(1-w)(1-h)                              (1)

pays |D|<=7/2. Here N\D is the disjoint union of the uncovered-full region
N\H(sigma+nu) and the complete residual region N intersect H(nu). Equation
(1) is bookkeeping from R63, not a new estimate.

The results read before choosing the new family were R63's exact deficit
identity; R68's support-tree and residual-free central rectangle; R69's
single-residual halfplane credit; R70's balanced-square carrier theorem and
near-four source-only obstruction; and R74's one-effective-outside theorem
and no-independent-pasting example. The present argument goes beyond that
one-position class: each effective corner has two positive outside atoms,
the source square has side .999, and the two strong corners at the upper
endpoint have the R80 actual half-local-mass obstruction. The new ingredient
is a bound on the *sum* of four joint-label costs using their shared side
capacities. We do not impose an equal receiver quota on the corners.

## 2. Fixed actual family and side capacity

From now on fix

    w=999/1000, e=1/1000, m=1/4, d=71/200, D0=171/200,
    sigma=(delta_(0,0)+delta_(w,0)+delta_(w,w)+delta_(0,w))/4.

Let T(x,y)=(w-y,x), the quarter-turn about the source-square center. For
i=0,1,2,3 put outside sources at

    q_i=T^i(-d,0), z_i=T^i(0,-D0),
    nu_t=sum_i [a_i delta_(q_i)+b_i delta_(z_i)],

where 1/4<=t<=499/1000 and

    i even: a_i=1/4, b_i=t;
    i odd:  a_i=1/2-t, b_i=1/4.                       (2)

All twelve source positions are distinct and all masses are positive.
The theorem is conditional on the full input preserving all square caps.
It neither assumes independent corner feasibility nor claims a unique
maximizer. The producer's separate exact cap-polytope certificate establishes
that the entire requested interval is feasible; both computed endpoints have
five maximum-block ties, one of them the old four-point block S. That finite qualification is
separate from the all-t analytic area proof below.

Each physical side has one q_i and one z_(i-1). A closed unit square
containing that side's two old endpoints and both those outside points exists:
on the left side use [-D0,1-D0]x[0,1], with rotations for the others. Thus
global feasibility requires the four shared side capacities

    a_i+b_(i-1)<=1/2, indices modulo 4.                (3)

The masses (2) saturate all four, giving total outside mass 2. They cannot
be chosen independently at four corners. Their local total masses are

    B_even=1/4+t, B_odd=3/4-t, sum_i B_i=2.            (4)

## 3. Actual local menus and the near contribution

The effective source window for the southwest corner Gamma=[-1/2,-e]² is
P=(-1,w-1/2)². Its rotations are pairwise disjoint. Every valid square at
a corner is strictly contained in its corresponding P because r<1/2.
Each P contains exactly one old corner atom of mass m and its actual two
outside atoms q_i,z_i. The local full and residual halos at all radii are
therefore the true local three-source and two-source halos. No old mass has
been collapsed. All other source positions remain in the global cap tests.

Use positive receiver coordinates X,Y and crop [e,1/2]² at one corner. Write
a for its near mass and b for its far mass, and put

    v=1/4, u=sqrt(m+a)/2, z=sqrt(a)/2,
    U=sqrt(m+a+b)/2, V=sqrt(a+b)/2.

The full seven labels are the four old-containing labels plus the three
pure-outside labels. The residual is exactly the latter three. The pure-old
receiver is (0,v)², the near mixed receiver is (d-u,u)x(0,u), and the joint
receiver after crop, if nonempty, is

    [e,U]x[D0-U,U].                                   (5)

Indeed joint activity requires U>D0/2>d, giving the stated X crop. All
local full radii are below 1/2 under (2). The far mixed label contributes
no deletion: at even corners a=m, so its radius equals the residual pair's
radius V and their cropped receivers are identical; at odd corners b=m,
its radius is sqrt(1/8)<D0/2, so it is empty. The far residual singleton
misses the crop because b<=.499 and sqrt(b)/2<D0-1/2=d. The pair residual
cannot meet the pure-old or near mixed receivers: its lower Y endpoint,
when it is nonempty, is D0-V>.855-.433>.353554>=u. Also D0-U>=.355>u,
so the joint full receiver lies above the entire near part.

For every a in [.001,.25], u<=sqrt(1/8)<d-e, z>e, and u>=v. The actual
pure-old-plus-near deletion therefore has exact area

    G(a)=.213501+a/2-.606u-(u+z-d)_+(z-e).             (6)

Here and below terminating decimals denote exact rationals. To derive it,
the old-plus-near full union has area

    (2u-d)(u-e)+(d-u-e)(v-e)
       =(m+a)/2+dv-ve+e²-(d+v+e)u.

Its intersection with the complete local residual is just the near singleton
overlap, of area (u+z-d)_+(z-e). The pure-old thin strip is included. Thus
(6) subtracts the actual residual and is not a full-halo proxy.

At even corners a=m, the R80 first-component expression gives

    G(m)=182323/500000-(171/200)sqrt(1/8)<63/1000.     (7)

For odd corners, the uniform estimate G(a)<1/10 holds for every permitted a.
Here is an analytic proof, rather than a parameter grid check. If u+z<=d,
then z<=d-v=.105, so a=4z²<=.0441 and (6) gives

    G(a)<=.213501+.0441/2-.606/4=.084051<.1.

If u+z>d, expand the last term in (6), using z²=a/4:

    G(a)=.213146+z²+.356z-(.605+z)u
        <=.213146-.605u+.356z.

The identity u²-z²=1/16 and reverse Cauchy give

    .605u-.356z >=(1/4)sqrt(.605²-.356²)>.12,

because .605²-.356²=.239289>.48². The reverse Cauchy inequality follows
by squaring positive quantities; its squared surplus is
(.605z-.356u)²>=0. Consequently G(a)<.093146<.1. The four near components
are thereby paid jointly at

    sum_i G(a_i)<2(.063)+2(.1)=.326.                  (8)

## 4. Shared budget for the four joint labels

Let B=a+b and define B0=D0²-m=19241/40000=.481025. The joint full receiver
is empty if B<=B0, including equality. Otherwise its cropped area is

    f(B)=(U-e)(2U-D0), U=sqrt(m+B)/2.

For B in (B0,.75], differentiation is legitimate and gives

    f'(B)=1/2-(D0+2e)/(4sqrt(m+B))
          <=1/2-(D0+2e)/4=1143/4000=:kappa.

Since f(B0)=0, its area is bounded by kappa(B-B0)_+. The actual joint
deletion is at most this *full* joint area; its same residual pair can only
reduce it. This safe upper estimate does not replace the residual by a
chosen mass or omit that pair from the actual halo definition.

Write eta=t-1/4, so 0<=eta<=249/1000, and delta=1/2-B0=759/40000. By (4),

    sum_i (B_i-B0)_+
      =2[(delta+eta)+(delta-eta)_+]
      <=2(delta+249/1000)=10719/20000.

For eta<=delta the expression is the smaller constant 4delta; above delta
it is increasing, so the endpoint is the stated global maximum. Therefore

    sum_i actual_joint_deletion
       <=(1143/4000)(10719/20000)<31/200=.155.         (9)

This is the cross-corner budget: two strong joint masses force the other
two joint masses down, by the same side-capacity assignments. Maximizing
all four local joint costs separately would lose (9).

Combining (8)-(9), all t in (2) satisfy

    sum_corners |D intersect Gamma_i|<481/1000.

The noncorner container has area 3-3e²=2.999997. With no additional
noncorner credit needed, the complete actual deletion is therefore paid:

    |D|<3-3e²+.481=3480997/1000000<7/2.               (10)

All labels used the same nu_t. Full and residual outside-corner radii were
not independently chosen, and no per-corner m_P/2 claim was used.

## 5. A tempting affine local cost is not a free lemma

The initial uncapped candidate G_m(B)<=m/4+B/8 is false even with one outside
atom. Take m=1 at origin, outside mass B=1/100 at (-11/20,0), and the whole
southwest orthant. With u=sqrt(101)/20 and residual radius 1/20, direct
old/near rectangle subtraction gives the true deletion

    (161-11sqrt(101))/200 >201/800=1/4+B/8.

The strict comparison follows from sqrt(101)<443/44, certified by squaring.
This is an *uncapped* local obstruction: the closed unit square containing
both atoms has mass 1.01, so the old singleton is not the maximum block of
mass 1. It does not disprove a version requiring full global caps and actual
internal source positions. The successful proof above instead uses its
specific near geometry and the shared side masses; it does not invoke this
false affine shortcut.

## 6. Verification and remaining general problem

The author runs only `computation/round081-dual/check_shared_side_constants.py`
on `computation/inputs/r081-dual/shared-side-constants.json`. It uses exact
Python Fraction arithmetic to verify the displayed constants, the radical
comparison in (7), the reverse-Cauchy rational gap, the two hinge branches,
and the strict aggregate and affine-counterexample comparisons. Its result
and execution receipt are in `computation/results/round081-dual/checks.json`.
It does not enumerate outside positions, solve a cap polytope, integrate a
halo union, or perform a uniform subdivision or B&B in t. The analytic
inequalities (6)-(9) prove the continuum statement.

The separate producer owns the complete all-t cap-feasibility certificate
and two actual endpoint halo computations. Those sources are pinned in the
component manifest as imported qualification/evidence, not rerun author
work. Maximum-block ties are retained. The conditional area proof requires
only that S be a maximum block, and applies to that S despite the ties.

For arbitrary source bboxes, internal mass profiles, outside positions or
unrestricted outside rank, (1) remains an unpaid general requirement.
Neither (6)'s near geometry nor (4)'s two-pair mass relation is universal.
The general maximum-block selection theorem and the resulting 7/2 weak
bound remain open. No MATLAB, agent creation, commit or push is performed.
