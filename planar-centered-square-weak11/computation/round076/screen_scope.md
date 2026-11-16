# R76 one finite two-outside screen with full S3 capacities

Actual date:2026-11-16, Asia/Shanghai. One purpose-designed batch tested
the auxiliary southwest-corner budget Delta_Gamma<=1/18 for S3 plus
two arbitrary candidate outside positions. No threshold signal occurred.
The batch was stopped; no repeat, parameter refinement, LP or proof
of the general auxiliary budget was inferred.

The internal measure is the full S3={0,4/9,8/9}^2, each mass1/9.
Candidate points lie strictly in P=(-1,7/18)^2 and outside closed
J=[-1/9,1]^2. Points are distinct and never coincide with S3.
The full caps keep S3 maximum of mass1. Consequently every strict
full/residual witness has radius<1/2, and only the internal origin
can contribute on Gamma=[-1/2,-1/9]^2. The local origin mass1/9 is
not treated as a global maximum block.

## Position and mass coverage

Coordinates use denominator900. The negative-coordinate levels are

    -899,-875,-750,-600,-500,-450,-400,-300,-225,-150,-125,-101,

and the other arm-coordinate levels are

    -75,0,75,150,225,300,325,349.

Thus the336-point pool consists of96 left-arm points,96 bottom-arm
points and144 points with both coordinates below the collar. It includes
critical scales and coordinates1/900 inside the open P boundary and
outside closed J. Pair seed760076 selects unique unordered pairs:
96 mirror cross pairs,480 further cross pairs,128 same-left,128
same-bottom,128 corner-plus-arm,64 corner-corner. There are1024 pairs.
This is a finite rational pool, not a continuous-position sample or
exhaustive pool-pair enumeration.

For each fixed pair, every nonempty full S3-plus-outside bbox closure
is enumerated. A row with outside tag0,1,2,3 has variable coefficient
0,b1,b2,b1+b2, respectively. The tightest exact right endpoints produce

    0<=b1<=V1, 0<=b2<=V2, b1+b2<=U.

No unit-square-only replacement is made. Every row uses the actual
maximum bbox span and right side max(1,span^2)-internal_mass. All exact
vertices of this polygon are retained. Its vertex centroid, half
centroid, and the two midpoint directions toward coordinate-biased
vertices add at most4 interior mass points; duplicates are removed.
Each pair has at most9 mass points. Every sampled point is rechecked
against **all** original bbox rows, not only the three compressed rows.
The actual batch has8143 mass points and581126 exact cap checks.

## Actual halo calculation and result

For each mass point, all7 nonempty local full-source subsets and all3
residual subsets are used. Eligibility mass>span_squared is decided
exactly; square roots, receiver endpoints and clipped union areas are
floating. The same outside vector defines the complete residual halo.
Every rectangle is clipped to the whole Gamma. No R72 whole-cell
surrogate, chosen-radius menu, or omitted triple is used.

Duplicate *actual measures* are cached after zero masses are removed
from the cache key. The8143 candidates correspond to5392 distinct
actual measures;2751 repeated area evaluations were avoided. The
remaining exact-cap checks still occur for every candidate. Positive
mass support, not unused candidate positions, defines the cache key.

The floating trigger was1/18+1/100000. None occurred. The largest
displayed deletion area is

    .03812305031872272,

with displayed margin.01743250523683283 to1/18. Its exact sampled
positions and masses are

    q1=(-2/3,-1/12), b1=4/9,
    q2=(1/6,-5/9),   b2=1/3.

The exact polygon bounds there are V1=46/81,V2=1/3,U=7/9. The
floating full/residual corner areas are.11031369414423534 and
.07219064382551264. These are screening values, not a certified
new area lower bound or an extremum. Since there was no trigger,
no further search or counterexample-certificate computation was run.

## Triple-label grouping

For the bbox[a,b]x[c,d] of origin and the two candidate points, define

    delta=max((b-a)/2,(d-c)/2,b+1/9,d+1/9,
              -1/2-a,-1/2-c,0).

The formal triple label has a strict Gamma witness exactly when
1/9+b1+b2>4delta^2. This predicate is recorded exactly, but **does
not filter** any numeric sample. For746 sampled position pairs,
1/9+min(U,V1+V2)<=4delta^2, so their whole feasible mass polygon
has no such formal triple witness. This is an exact classification
of those fixed position pairs, not all possible positions.

The saved grouping separates actual positive outside arity, since
zero-mass candidate locations can make an irrelevant formal label
active. Counts are:

| positive outside arity | triple null | triple active |
|---|---:|---:|
|0|1009|15|
|1|1731|317|
|2|4136|935|

The two-positive active group gives the overall maximum above.
The two-positive null group's largest display is.018459142325719626.
The independently proved triple-null payment is owned by continuous/core;
these finite floating comparisons are not presented as its proof.
Grouping rereads saved data only, with no new pair, mass or halo evaluation.

## Files and execution

Python3.13.18 and standard-library integer/Fraction arithmetic were used
for caps and generation. The batch took1.068219 seconds internally,
PID42374 exited0; the separate grouping PID46453 also exited0. Both
were checked absent. A single process was sufficient for this small
workload; no additional agents, MATLAB, workers or optimization were used.

The parameter JSON specifies all coordinate levels, seed, category
counts, mass rule and trigger. The compressed data file contains every
actual pair's V1,V2,U, vertices, mass candidates and area displays.
No old R75/R74/R73 file was modified; previous manifests and files are
hash-checked unchanged. No commit or push was performed.

After creating output directories, reproduce from the topic root:

    python computation/round076/two_outside_capped_screen.py --input computation/inputs/r076/two_outside_screen_parameters.json --out computation/results/round076-two-outside/summary.json --data computation/results/round076-two-outside/screen-data.json.gz
    python computation/round076/group_screen_statistics.py --input computation/results/round076-two-outside/screen-data.json.gz --out computation/results/round076-two-outside/group-statistics.json

This single negative screen does not prove the auxiliary bound,
the four-corner budget, the general7/2 deletion budget, or exact C2.
