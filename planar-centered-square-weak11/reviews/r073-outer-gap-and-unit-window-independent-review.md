# R73 independent review: outer credit, exact caps and bounded window weights

2026-11-16. Sol, pcw-sol-cover-structure. Only R73 independent files were
written. R72 publication 7c3d261555784e16ee5cf6f2d757612589dd0a92 is the
inherited archive; old inputs/certificates were not modified or rerun merely
for provenance. No agent was spawned and no commit/push performed here.

## 1. Target and the missing interior credit

For the uniform S_n of mass1, write e=1/n², w=1-e,
C=[-1/2,w+1/2]², J=[-e,1]², G=C minus J, and
D=H(sigma+nu) minus H(nu). Maximal eligible mass1 gives D subset C.
The stronger proposed statement is

    |G minus D| >= |C|-7/2.                            (1)

It remains unproved and unrefuted. For every n, exact accounting is

    |C|-|D|=|G minus D|+|J minus Hfull|+|J intersect Hres|.

Both missing-deletion components in each region are disjoint: vacancy
outside Hfull and residual halo inside Hfull. For n>=5, J subset H(sigma),
so its vacancy is zero. The original budget permits the positive interior
residual credit; (1) discards it. Nu(J)=0 does not imply H(nu) misses J.

I independently derived strict J coverage for n>=5. Use open singleton
bands of half-width1/(2n). The intervening closed gaps have width e.
One gap is covered by a two-source label because n+2<sqrt(2)n; two gaps
are covered by the four-source square because n+2<2n. Outer J endpoints
are strictly in the first/last singleton bands. Thus there is no lost
closed boundary. For n3 and n4 the new exact controls (9/50,0) and
(13/100,0) have strictly positive nonhalo slack at EVERY complete old
distance prefix, proving genuine interior holes. No J-coverage assertion
for those n is imported from the large-n case.

A new admissible27-atom example shows the interior residual term is real.
Use S5 and add q1=(-1/20,12/25), q2=(6/5,12/25), each mass1/50. Any proper
old subset plus both outside atoms has mass<=1, with equality impossible
for an eligible set since the outside pair span is5/4. Full S plus either
or both outside atoms is ineligible by exact spans. Hence S5 is the unique
maximum. Eligible attachment pairs connect the input; its nonterminal gap
is1213/2640>0. The residual halo consists of two disjoint singleton squares
of radius sqrt(2)/20. Only q1 enters J, in exact area

    |J intersect H(nu)|=(10-sqrt(2))/1000>0.

The new Fraction run checks all326 bbox closures, uniqueness, connectivity,
nonterminality, strict coverage inequalities and the two small-n holes.
PID26195 exited0 and was absent. This example refutes the zero-interior-
credit inference, not (1) or the actual7/2 target.

## 2. All-parameter cap proof and the local-full-halo obstruction

I read the complete author proof
`r073-cap-preserving-atomization-and-local-full-halo-obstruction.md`, final
SHA53f99433f9fcdb2d8cc5ebb46ec505c92bbe088c8e8b3b223d303495e1100349.
It is correct, including all edge cases. Selected complete delta cells
avoid closed J_delta; center capture is bounded by uniform mass in the
side-(R+2delta) expansion. If old axis counts are positive, expanding both
intersecting intervals increases overlap by2delta. With c=t/(1+2delta),
G0=(n-k)(n-l)/n³, the exact excess factorization is

    -(1-2n delta)G0/(1+2delta)
      -(R-1)[(1-c)(R+1)-4c delta].

The second bracket equals(1-c)(R-1)+2(1-t)>0. Nonpositive excess follows
for R>=1 and delta<=1/(2n). Zero old capture uses t(1+2delta)<=1; smaller
squares are included in a unit square. Closed ties are retained. The
argument works for any bounded truncation, not only the original square.
It proves exact maximum mass1 at each mesh; outside contraction makes S
unique, not merely approximately maximizing.

Independent Python3.13.18/SymPy1.14.0 execution checked six formal rational
identities by exact zero polynomial numerators, and positive-polynomial
forms for n>=5 and n>=7. It also exactly recomputed the explicit n9
counts, mass, nonterminal gap and area lower bound; PID29335 exited0 and
was absent. This is distinct from the author's3092/9276/42 finite Fraction
controls, which were source/code/result reviewed, not repeated here.
General signs and geometry were mathematically reviewed, not inferred
from finitely many parameter samples.

The n9 formula specifies79056 actual positive atoms. Its local full-halo
lower bound7832077922743/2230740000000>7/2 is valid, as are unique S9,
connectedness and nonterminality. Neither I nor the author enumerated
these79056 atoms, their subsets or a complete halo. The same explicit
construction already has exterior residual singleton area
28179*29997/340000000, which exceeds |C|-7/2=5915/13122. Thus it SATISFIES
(1). It rules out a local-full-halo replacement, not actual deletion or
the stronger outer-gap statement. A useful sufficient certificate for
(1) is simply a certified residual inner area in G at least |C|-7/2;
adding a disjoint certified full-halo vacancy only improves it.

## 3. Shared kernels and the full finite-domain qualification

I read all six sections of
`r073-shared-unit-window-integral-kernels.md`, final proof
SHA8b03a1777f86d50fb33eb73e6f0e980b4d0ca4f86d3eedd51b7d1e38b3b18b10.
The exact nonnested L decomposition, Tonelli kernel budget, tent kernel,
reflected-family overlap correction and admissible non-CDF example are
valid. The budget keeps ONE nu throughout; it does not give independently
chosen corner masses. The26-atom example's counts/margins are correct.

The new finite construction in Section5 is fully qualified. It uses
Brect=[-1,3]x[-1,2], delta=1/(100n²), lambda=1-delta², n>=5. The expanded
hole lies strictly inside the lattice rectangle. Its complement is a
four-neighbor connected frame; 2lambda c>1 makes neighboring pairs
eligible. The selected left bridge has distance41/(40n²)<1/n from an
internal source. Total mass<=13 and span>=399/100 give a uniform positive
nonterminal margin. The original square B, with n growing, could not
justify this domain assertion; the rectangular repair is essential.

The closed-unit-window estimate was independently checked in detail.
Center replacement changes mass only in the centered side-(1-delta) to
side-(1+delta) boundary strip, area4delta, including centers on a closed
edge. Missing cells are in the expanded J boundary layer, area<=14delta;
the density change costs<=3delta. Thus the uniform21delta error is valid.
A unit interval cannot have both endpoints strictly in(0,1); the grid
CDF and its left limit therefore give count/length error<=1/n even at
closed ties. Changing[0,1] to[-e,1] adds at most e. These estimates give

    0<=g_n(z)<=4/n+21/(100n²), every z in C_infinity=[-1/2,3/2]².

Consequently ANY finite positive window-parameter measure rho_n supported
there, of total mass<=M, has slack integral tending0, even when rho_n
depends on n and the entire input and is atomic. Pointwise control comes
before integration; no unjustified transfer of discontinuous kernels is
used. The center-domain and bounded-total-weight restrictions matter.
This excludes only extracting a uniform positive slack using that class
of weights; it does not exclude growing weights, radius information or
actual residual-halo credit.

The author's441 product cells and961 corner limits were code/result/source
reviewed and NOT independently rerun. Its n5 enclosure[0,9/25] is a
corner-limit upper enclosure, not an assertion of actual attainment.
The68,219,184-source construction is an index formula; no huge source
array or halo catalogue was generated. No new Borel or R72 leaf run was
performed here. The inherited R72 strict Borel lower bound remains an
explicit source/certificate dependency.

## 4. New regional and complete-geometry executions

For one selected control, the R72 contracted unique38-atom input, I newly
reconstructed exact closures and directed roots from the RAW ATOMS and
clipped them to all nine C/J regions. Horizontal y-slabs with x-interval
merging matched all90 regional area endpoints exactly. Independent
source-origin projection masks reconstructed all163 nonempty closed unit
capture patterns, with one saturating pattern. PID28902 exited0 and was
absent. This was a new regional audit, not an unclipped R72 full-area rerun.
The other three regional controls were source/code reviewed and bound;
their region areas were not independently recomputed here.

I identified the producer's original orientation label as wrong: its
events are x endpoints with merged y intervals. Final v2 correctly calls
this vertical x-event interval merging. Its actual v1/v2 semantic-equality
record is source-bound; only v2 is archived. My area calculation uses the
other orientation. This was a method-description correction, not a change
in geometric formulas or results.

Root's new29-atom four-axis input was independently evaluated in full
using the READ-ONLY frozen R72 horizontal-slab checker, as requested.
All476 closures/171 eligible closures, the unique S5, connectivity,
nonterminal gap733950007/515625000, directed rectangles and every exact
full/residual/internal/deletion area endpoint matched. PID26699 exited0
and was absent. The actual deletion is approximately1.5826878638260273;
its exact lower bound is strictly above79/50. Input/source/result pins and
the post-execution exact bound comparison are saved in the R73 receipt.
This is a new candidate execution, not an old certificate relabelled.
It does not improve the general C2 bound or threaten7/2.

## 5. Final evidence scope and reproduction

The final report version reviewed is
eede9b61ed3fbd4eeb9d937fe3f1e5e11efd0784148f98eb9994a6be0cd8bdfb.
It correctly limits the weighted-slack theorem to the relevant center
domain and bounded weight mass, distinguishes actual deletion from local
full halo, and records formulas rather than enumeration for large inputs.
No mathematical changes remain requested. Hash binding is only version
verification, separate from these readings and actual executions.

All four new independent runs retain source hashes, actual software
versions, exact outputs, commands and exit/PID receipts. Root's D4 probes,
Wolfram checks and the two authors' small kernels are source-bound evidence;
I did not rerun or claim them as my executions. No MATLAB/Wolfram process
was launched by this workspace. The stronger outer-gap and corner-only
actual-deletion inequalities, arbitrary finite-outside7/2 budget and exact
C2 remain open.

From the integrated topic directory, use explicit new data outputs:

    python3 -B computation/review_r073_core_credit.py --input computation/inputs/r073-core-credit.json --output /Volumes/data/CodexTemp/T2/r073-replay/core-credit.json

    python3 -B computation/review_r073_cap_identities.py --output /Volumes/data/CodexTemp/T2/r073-replay/cap-identities.json

The second command requires SymPy; the actual run used the preinstalled
main environment recorded in its receipt. Regional audit:

    python3 -B computation/review_r073_regional.py --input computation/inputs/r072/contracted-unique-actual-input.json --saved computation/results/round073-regional/four-control-regions-v2.json --producer computation/round073/regional_credit_diagnostic_v2.py --helper computation/review_r072_actual_halo_y_slabs.py --case contracted-unique --output /Volumes/data/CodexTemp/T2/r073-replay/regions.json

The29-atom audit uses the unchanged archived R72 checker on the NEW input:

    python3 -B computation/review_r072_actual_halo_y_slabs.py --input computation/inputs/r073/s25-four-axis-unique-input.json --saved computation/results/round073-root/s25-four-axis-unique-area.json --producer computation/round072/actual_halo_geometry.py --bits 72 --output /Volumes/data/CodexTemp/T2/r073-replay/four-axis.json

Root owns final publication and remote verification.
