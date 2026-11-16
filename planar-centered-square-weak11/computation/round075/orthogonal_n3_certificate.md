# R75: complete n=3 orthogonal equal-outside parameter certificate

Let S3={0,4/9,8/9}^2, each mass1/9, and keep the full grid sigma3.
The outside measure is

    nu_(a,b)=b delta_(-a,0)+b delta_(0,-a),
    0<a<=1, 0<=b<=B3(a),
    B3(a)=1/2 min_(k=1,2,3){max(1,[a+(k-1)4/9]^2)-k^2/9}.

The separately frozen core theorem (SHA7a0e6c14...) proves this is exactly
the all-square-cap condition keeping S3 maximum of mass1. This component
retains that condition; it does not use the false uncapped two-outside
collapse inequality. On Gamma=[-1/2,-1/9]^2, the only effective internal
source is the origin, mass m=1/9. Strict full/residual witnesses all have
radius<1/2 by the full caps.

The completed exact parameter tree certifies

    |[H(sigma3+nu_(a,b)) minus H(nu_(a,b))] intersect Gamma|
       <= Phi < 1/18,

uniformly for every admissible parameter, where

    Phi=1563026333441017474428935886295754233898431
          /28224380641890319913506143438820582531006464
        approximately .055378587515263.

Phi is a certified uniform upper bound, not an optimized or attained
maximum. The maximum sampled midpoint area was only about.0274303432.

## Midpoint screening and exact caps

One31 by31 midpoint batch sampled a in(1/9,1) and relative mass
t=b/B3(a) in(0,1). It therefore tested many b values at each a and did
not assume deletion is monotone in b. Every point's65 complete global
S3-plus-outside bbox caps was checked with Fraction arithmetic,62465
checks total. All7 local full-source and3 residual source subsets were
used in the floating complete-halo corner evaluation. The largest
displayed area was.0274303431731025, at a187/279 and b3789625/9652284.
The display's margin to1/18 was.0281252123824531. These961 samples
were a screen, not a continuum proof, and were not independently replayed.

## Uniform boxes and exact endpoint directions

The root parameter rectangle is[1/9,1]x[0,4/9]. Because B3 is nondecreasing,
a box[aL,aU]x[bL,bU] is infeasible if bL>B3(aU). Otherwise bU is clipped
to min(bU,B3(aU)). Every removed upper-b strip is recorded with its
lower edge open, keeping the equality b=B3(aU) in the retained box.
This uses the complete min-max cap, including nonunit square rows.

For each retained parameter box, dyadic directed72bit bounds are used:
uU bounds sqrt(m+bU)/2 above, vU bounds sqrt(m+2bU)/2 above,
sL bounds sqrt(bL)/2 below, and tL bounds sqrt(2bL)/2 below.
The full-label outer rectangles are

    [-1/6,1/6]^2,
    [-uU,uU-aL]x[-uU,uU],
    [-uU,uU]x[-uU,uU-aL],
    [-vU,vU-aL]^2.

The common residual inner rectangles are

    [-aL-sL,-aU+sL]x[-sL,sL],
    [-sL,sL]x[-aL-sL,-aU+sL],
    [-tL,tL-aU]^2.

All are clipped to Gamma; empty rectangles are omitted. Their exact
full-union minus residual-union area is a uniform upper bound for every
actual parameter in the box. The core proof establishes the pointwise
inclusions: this is an enclosing calculation, not an independent choice
of full and residual outside measures. All decisions about cap clipping,
roots, areas and comparison with1/18 use integer/Fraction arithmetic.
Floating heap priorities only determine which unresolved box is visited next.

Closed child boxes bisect the longer coordinate and cover their retained
parent, sharing only the bisection boundary. The certificate gives every
node's original/retained box, cap, removed strip, split midpoint, children,
four root records, clipped rectangles, exact upper area and status. Thus
the entire parameter domain is traceable rather than represented by a
list of successful leaves. For0<a<1/9, the cap forces b=0; the pure
internal Gamma area is exactly1/324 and is handled separately.

## Actual complete cover

The single authorized run used71 nodes:35 splits and36 certified leaves.
Seven upper-b strips were removed by exact caps. Maximum depth was9;
there were no unresolved branches and no node/depth budget exhaustion.
The predeclared limits5000 nodes/depth24 were not extended. Computation
stopped when every terminal bound was<=1/18, without refining toward a
numerical maximum. The smallest target margin is

    4994813330666965210294304749833684490817
      /28224380641890319913506143438820582531006464 >0.

The earlier uncapped d2 example remains frozen and outside this class.
For example a=2/3 requires b<=32/81 by the full cap; the naively scaled
uncapped choice b=4/9 violates a nonunit-square row. This positive family
certificate neither reverses the uncapped counterexample nor proves an
arbitrary outside-position or unequal-mass corner bound.

## Independent replay and global consequence

The audit agent independently rebuilt the full71-node tree,35 splits,
all7 cap strips,284 root bounds, every node rectangle, and all71 area
bounds using horizontal y-slabs and merged x intervals rather than this
producer's vertical x-event sweep. It did not import producer code,
regenerate the adaptive tree, or rerun the961-point probe. All572 audit
slabs and final Fraction bounds agreed; the separate audit output is
source-pinned in this component manifest.

For the specified two-outside family, the other three corners have no
effective outside source and each has pure internal area1/324. The
remaining part of the deletion container may be paid by its whole area:

    |C minus the4corners|=80/27,
    complete actual Delta <=80/27+1/18+3/324=109/36<7/2.

This is a payment for this full-grid n3 orthogonal equal-mass family,
without an added connectivity/nonterminal hypothesis. The consequence
uses the frozen capacity/localization theorem and the verified parameter
cover. It is not a new general C2 bound or a result for general two
outside positions, unequal masses, or n>=5.

## Files and execution

The producer uses Python3.13.18 and standard-library integer/Fraction
arithmetic. The midpoint run took.173633 seconds internally and exited0;
the one parameter-tree run took.023827 seconds and exited0. PIDs16656
and23767 were verified absent. No agent, MATLAB, LP, parameter search
extension, old-bundle edit, commit, or push was added by this component.

Reproduce after creating output directories:

    python computation/round075/orthogonal_n3_midpoint_probe.py --input computation/inputs/r075/orthogonal_n3_probe.json --out computation/results/round075-orthogonal-probe/midpoint961.json
    python computation/round075/orthogonal_n3_parameter_boxes.py --input computation/inputs/r075/orthogonal_n3_boxes.json --theory /absolute/path/to/r075-orthogonal-two-outside-full-grid-caps.md --out computation/results/round075-orthogonal-boxes/summary.json --certificate computation/results/round075-orthogonal-boxes/parameter-tree.json.gz

The source-pin hash ensures the exact core interface used by the run.
The compressed tree is the full certificate, not a cache or a sample.
