# R74 corner localization and one finite mass-variation pilot

Actual date:2026-11-16, Asia/Shanghai. The earlier Borel clipping
preliminary package with manifest SHA3a52e864... remains byte-for-byte
frozen. This new component contains the requested finite source-window
validation and one small true-corner mass experiment. It does not
continue the R72 whole-cell MIP or expand positions, starts, or n.

## Strict source windows

For S_n, e=1/n^2 and w=1-e, the southwest receiver corner is
`[-1/2,-e]^2` up to area-zero boundary assignment, and its source
window is the **open** square `P_SW=(-1,1/2-e)^2`.
Every strict full or residual witness has radius r<1/2: a witness
at r>=1/2 would violate the common square cap at side2r>=1.
For x in the closed receiver corner and a captured source p,

    -1 < x_i-r < p_i < x_i+r < 1/2-e

with the outer inequalities interpreted as weak x bounds followed by
the strict r<1/2 bound. Thus every actually captured source is strictly
in P_SW. Sources on its boundary cannot participate in a strict witness.
The other three P windows are reflections about w/2; their one-axis
intervals are `(-1,1/2-e)` and `(1/2,2-e)`, separated by a positive gap e.

Only outside mass is clipped; the full S is retained. Removing outside
mass preserves every square cap, hence cannot create a large-radius
witness. It neither removes nor adds any full or residual witness in
the corresponding corner. Pointwise localization is the geometric
argument; overlapping numerical area enclosures alone would not prove it.

Four inherited actual controls (S25, R71 positive33, R72 aligned38,
R72 contracted unique38) were tested in all four corners. There are
12 distinct restricted exact inputs, with11 new complete geometry
evaluations and reuse of identical S-only geometry. In all16 cases,
the restricted full-halo, residual-halo and deletion intervals overlap
their original certified corner intervals. The four SW deletion values
are approximately.0036,.0081683372,.0124298247,.0124298245. Every old
control has zero positive sources on the P boundaries, explicitly
recorded; the strict boundary assertion is supported by the radius
argument, not falsely attributed to these boundary-free controls.

The frozen R72 complete geometry engine and frozen R73 v2 clipping
functions are read-only imports with bytecode writing disabled.
This is new geometry only for genuinely restricted actual inputs;
it is not described as an independent complete-halo audit.

## Two explicit finite seeds

Global caps always retain S25 and its mass1. Two seeds selected with
Astra were completely evaluated before any mass directions:

| outside atoms | masses | true SW deletion |
|---|---|---:|
| (-12/25,-12/25) |16/25|.0223101428319403|
| (-12/25,6/25),(6/25,-12/25) |8/25 each|.0287820999749009|

All actual bbox caps were checked exactly. Both seeds have L=1 and
are connected/nonterminal. Neither violates the stronger auxiliary
threshold `m(sigma intersect P_SW)/2=2/25=.08` or the actual single-corner
target `1/8+3/(4*5^4)=631/5000=.1262`.

## One selected seed, exact true-area acceptance

Only the better cross seed was varied. Fixed outside positions are the
32 R71 grid positions strictly inside P_SW plus q_new=(-3/10,47/200),
33 total. All global maxblock constraints against the full S25 give
487 deduplicated linear cap rows. The full fixed-source closure catalog
has1648 members. Complete full and residual halo unions are used for
the **actual** SW deletion objective, not a fixed-cell coverage score.

Finite differences at step1e-5 propose a gradient. One HiGHS LP per
direction, trust radius1/16, proposes feasible mass changes. Four
rational step fractions1,1/2,1/4,1/8 are tested. Masses are rounded to a
common denominator1e8 and uniformly shrunk if an exact cap requires it;
every cap is then checked with Fraction arithmetic. Each rational line
point is evaluated by the complete72bit full/residual geometry. A move
is accepted only if its new lower area bound exceeds the current upper.
The floating gradient and LP success status carry no proof of increase
or optimality.

The three-direction delivered batch accepts two1/8 steps:

    .0287820999749009 -> .0288947144617460 -> .0291471246353845.

All four steps of the third proposed direction fail strict improvement.
The full steps of the first two directions also reduce the true area;
this is a concrete nonlinear overstep in the linear proposals, rather
than LP infeasibility. The selected small batch was stopped. It does
not certify local optimality or exclude smaller steps or other inputs.
The gain supplies no near-threshold signal for.08 or.1262, so no
positions, starts, or n9 experiments were added.

The final exact input has S25 plus five positive outside atoms:

| position | mass |
|---|---|
|(-12/25,6/25)|487/1600|
|(-6/25,-6/25)|1/128|
|(-6/25,0)|1/128|
|(0,-6/25)|12207/6250000|
|(6/25,-12/25)|487/1600|

It has30 positive sources,447 bbox closures,156 eligible closures,
one eligible component, and **unique** maximum S of mass1. Actual span
is36/25, total mass40658203/25000000, and nonterminal gap
57536579/165000000>0. Its corner full halo is about.0434190450,
residual about.0142719204, and true deletion about.0291471246.
The full global S-deletion is about1.4880770442. All decimals are
displays of saved directed rational bounds. This is a new explicit
finite control, not a violation of either proposed corner budget or
a new C2 lower bound. There is no strict threshold violation requiring
an independent second area algorithm to validate a counterexample.

## Actual execution and provenance

Python3.13.18 and, for mass proposals, NumPy2.5.4/SciPy1.18.1 were used.
LP feasibility tolerances were1e-9, requested HiGHS threads1, and
BLAS/OMP one thread. All source files and frozen imports are SHA-pinned.

The first mass-variation run completed its calculations but failed to
write the final output because its result directory had been created
in the wrong workspace parent. After correcting the directory, the
unchanged script was replayed successfully. Both actual executions
are retained, including the failure receipt; accepted waypoint JSON
records match exactly. Thus the delivered experiment has3 direction
calls/12 line evaluations, while actual session totals include the
replay and are6 LP calls/24 line evaluations. The failed run is not
presented as another seed or independent confirmation.

All recorded processes have exited. The accidentally created empty
parent directories are removed only if still empty; other files are
untouched. R73 frozen v2 and the frozen Borel preliminary package are
hash-checked unchanged. No MATLAB, commit, push, reset, rebase, or old
round edit was performed.

Reproduce from topic root with data temporary/cache settings:

    python computation/round074/verify_corner_source_windows.py --archive /absolute/integrated/topic --out computation/results/round074-localization/four-control-window-check.json
    python computation/round074/single_corner_two_seed_pilot.py --archive /absolute/integrated/topic --out computation/results/round074-pilot/two-seed-results.json
    python computation/round074/corner_mass_variation.py --archive /absolute/integrated/topic --out computation/results/round074-mass-variation/one-seed-three-direction.json

Create the corresponding output directories first. Input/output JSON
and the component manifest give exact atoms, all accepted and rejected
trial weights, true-area endpoints, caps, parameters and receipts.
