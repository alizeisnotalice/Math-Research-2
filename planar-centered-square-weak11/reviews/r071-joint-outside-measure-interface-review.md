# R71 independent interface audit: one common outside measure

2026-11-16. Sol independent reviewer, pcw-sol-cover-structure. R70 and all
earlier rounds are frozen. Root's published R70 archive is authoritative;
this review owns no author/producer file, commit or push.

## 1. Quantifiers and exact finite formulation

Fix a finite positive old measure mu supported on S, with mu(S)=1 and
strictly eligible S (its two bbox spans are less than 1). Fix distinct
candidate outside positions q_j not in S; duplicate candidates must first
be merged. Unknown z_j>=0 define one common nu=sum_j z_j delta_(q_j).
Zero candidate masses are allowed. Every geometric domain conclusion must
use the resulting positive support, not all candidate locations.

For every complete source-coordinate bbox closure B in S union {q_j},
let d_B be its tight L-infinity diameter. The full constraints

    sum_(q_j in B) z_j <= max(1,d_B^2)-mu(B)              (1)

are exactly equivalent to maximum eligible mass at most 1, by R64. Since
the fixed positive S is eligible with mass 1, (1) is equivalent to its
remaining a maximum block of mass exactly 1. Closures of *all* subsets,
including those ineligible at a seed, are needed. Zero candidate locations
can remain in the closure list: completing a violating positive subset
to its tight bbox closure preserves its bbox and only increases mass.
They cannot justify connectivity, fatness or nonterminality from the
candidate bbox. Unique maximality is a separate optional requirement.

For a prescribed receiver x let R range over all distinct distances
||q_j-x||_infinity, including R=0 when present. Then

    x outside H(nu)
    iff sum_(||q_j-x||<=R) z_j <=4R^2 for every R.       (2)

All ties are included in the same closed prefix. Necessity follows because
a strict excess at R remains a strict excess for a slightly larger open
radius. Conversely, a strict open-radius halo witness contains a finite
positive subset whose largest distance R is strictly smaller than its
witness radius; its complete prefix violates the corresponding bound.
Equivalently use the largest source distance in that captured set. At
R=0, (2) forces all outside mass exactly at x to vanish. Omitting zero
distance is invalid; omitting empty prefixes is harmless. Repeated
distance events may be deduplicated only after retaining every tied mass.

Fix one rational radius 0<r_k<1/2 for each of finitely many receivers x_k,
and define the open witness square W_k={p:||p-x_k||<r_k}. With
a_k=mu(W_k), d_k=4r_k^2-a_k, require

    sum_(q_j in W_k) z_j >=d_k+eta for every k.         (3)

If eta>0, (3) supplies an actual strict full-input halo witness at r_k;
(2) makes each point a true deletion receiver. If the implementation uses
closed capture ||p-x||<=r_k instead, a positive margin still gives a true
halo witness at some slightly larger radius r'_k<1/2. It does not prove
that the original open radius r_k captures its boundary sources. State
which semantics is used and retain the strict margin in either case.

For the fixed support and radii, (1)-(3), z>=0, eta>0 are necessary and
sufficient for simultaneous strict deletion membership with those open
witness radii. Sufficiency is constructive from the certified masses.
Necessity at those radii follows by taking the minimum of finitely many
positive witness margins. This is not necessity at arbitrary radii. For
arbitrary actual deletion points and this support, each has some witness
radius below 1/2, since any captured mass is at most 1; the quantifier is
"there exist radii", not a prescribed-radius menu. Infeasibility of one
fixed radius tuple or finite tuple menu does not eliminate other radii,
source positions, supports, or general measures unless the finite radius
menu has separately been proved complete for those fixed positions.

There is an exact complete finite radius menu for fixed S union Q:
for each receiver list every distinct source distance 0<=t<1/2, using
the CLOSED capture square at t and requiring strict excess over 4t^2.
Any true full witness captures at most mass 1 under (1), hence has radius
less than 1/2; its largest positive captured distance t is one of these
entries and has strict excess. Conversely a strict closed excess yields
a slightly larger positive open witness, still below 1/2. Thus the
finite union over one complete-menu choice per receiver is exact for
all actual radii on this fixed support. A certificate excluding EVERY
tuple of this complete menu can make an all-radius fixed-support claim.
This differs from testing an arbitrarily selected finite radius menu.

## 2. Exact Farkas/dual interface

Stack (1)-(2) as Bz<=b, (3) as Cz>=d+eta*1. All coefficients are rational
for rational positions, receivers and radii. With eta free in R, consider

    maximize eta, subject to z>=0, Bz<=b, Cz>=d+eta*1.  (P)

For a legal fixed old block, b>=0 and z=0 with sufficiently negative eta
is feasible. Every candidate singleton gives z_j<=1, so the mass domain
is compact and eta is bounded above. Hence the optimum is finite and
attained. The strict system is feasible iff eta_star>0. A floating solver
failure or tolerance-level eta is not such a conclusion.

Any rational y>=0, lambda>=0 satisfying

    sum_k lambda_k=1, B^T y>=C^T lambda               (4)

is an independently checkable dual certificate, with

    eta <=b^T y-d^T lambda.                           (5)

Indeed multiply the lower rows by lambda and the upper rows by y;
nonnegative z and the column dominance in (4) give (5) directly. An
exact upper bound <=0 excludes strict feasibility. Positive feasible
masses and eta must instead be substituted into every row, with exact
arithmetic. A primal/dual equality certifies an optimum but is unnecessary
for an exclusion. For a specified required margin eta0, the unnormalized
Farkas alternative uses y,lambda>=0 with B^Ty>=C^Tlambda and

    b^Ty < lambda^T(d+eta0*1).                        (6)

The difference in (6) is the certified contradiction. Column equality is
not required because z is nonnegative. If the common margin variable is
constrained nonnegative, its dual normalization changes; do not copy (4)
without matching the primal convention. (P) uses a free margin explicitly.

## 3. A support-independent shared-measure cut

The same weak-duality argument can give an actual geometric compatibility
statement. Choose finitely many global axis-parallel closed rectangles
B_l, with diameter D_l, and residual centered closed squares Q(x_k,R),
with nonnegative weights y_l and v_(k,R). Let lambda_k>=0, with at least
one positive lambda. Suppose the *pointwise* inequality

    sum_k lambda_k 1_(W_k)(q)
      <=sum_l y_l 1_(B_l)(q)
         +sum_(k,R) v_(k,R) 1_(Q(x_k,R))(q)            (7)

holds for every allowed outside location q, not just a finite candidate
list. Then every common nonnegative nu retaining maximum mass 1 and
excluding all x_k from H(nu), while making all W_k strict full witnesses,
must satisfy

    sum_k lambda_k(4r_k^2-mu(W_k))
      <sum_l y_l(max(1,D_l^2)-mu(B_l))
         +sum_(k,R) v_(k,R)4R^2.                     (8)

Integrate (7) against the same nu. Global rectangle caps follow from
maximality as in R64, and residual caps from halo exclusion. At least one
positive lambda makes the witness inequality strict. Thus a certificate
with left side >= right side rules out *every* common outside measure on
the allowed location domain at these witness radii. It is not restricted
to a candidate grid. The radius tuple remains fixed unless separately
covered. No finite or infinite strong-duality claim is used for (8).

This form applies to arbitrary finite outside measures, and also to any
positive Borel measure for which the stated rectangle caps and witnesses
hold. A finite candidate dual certifies (7) only at those candidates. To
promote it to (8), verify all nonempty arrangement cells of the selected
rectangles and squares, including edges and vertices: open witness
boundaries and closed upper-set boundaries have different membership.
Open two-dimensional cells alone do not control atomic boundary mass.

Single-receiver feasibility checks are essential diagnostics. A joint
negative certificate whose lambda has only one nonzero entry is already
a one-receiver obstruction. Even a certificate involving several entries
can be redundant; explicitly certify each requested single-receiver
system feasible before describing a new shared-measure incompatibility.
The algebraic toy z1+z2<=1, z1>=3/5+eta, z2>=3/5+eta has separately
positive-margin feasibility but joint eta<=-1/10. It illustrates dual
signs only and is not claimed as a geometric source configuration.

## 4. Independent exact controls and pending producer audit

`computation/review_r071_joint_measure_boundaries.py` independently tests
the bbox-cap and residual-prefix equivalences by literal subsets and
Fraction arithmetic, including zero masses, all ties and zero-distance
events. It also checks a boundary-capture example and the algebraic dual
toy. It imports no producer LP or geometric module. Its saved output and
execution receipt record the actual test counts and source hash.
The actual run checked 256 mass assignments, 3471 literal positive
subset eligibility tests, 2304 closure-cap rows, 1536 receiver cases,
20826 literal positive halo subset tests, 3072 complete tied prefixes,
and 384 positive zero-distance cases. It used Python 3.15.0 and exited
zero, with PID 8725 absent after wait. These are executed loop counts,
not the number of all potential subsets before short-circuiting.

The core author proof `r071-shared-outside-mass-and-residual-duality.md`
was independently read in full (sections 1-10). Its cap equivalence,
residual zero/tie conventions, free-margin dual, exact complete finite
radius menu, diffuse-measure weak cut and complete arrangement atoms
are valid. The finite-row atom reduction must intersect each atom with
the permitted domain K minus S before choosing representatives if the
outside support is required disjoint from S. Core incorporated this
domain clarification; checking the stronger whole-plane covering is
also sufficient. Its later closed-collar reduction is valid: the collar
is the union of four side-one closed squares containing all S, each with
zero available outside mass. Missing columns wholly within that collar
can be covered at zero dual cost. The added receiver-cell minimum-radius
formula and exact all-residual-closure exclusion are also correct; their
full-label maximum-radius condition supplies an area lower certificate,
not a general upper bound. A positive full witness at a point does not
by itself imply a residual-free neighborhood.
The frozen final author proof is SHA-256
`d61554809a35ef393b735985fea4c8732a463df47ba51ff69ef837ca86fc9770`.
Its final core packet is
`c7e60eda2c36787889ab43d0fc4269b69e22369bd5089ff6ccb432a9042c3c06`.
The author's corrected two finite interface sweeps and receiver-cell
minimax kernel were code/result/receipt reviewed and source-bound, not
independently rerun here. The initial author's short-circuit count was
identified during this review; final source explicitly evaluates every
subset before applying all, and the two new full runs are separate from
the retained initial materials. Core records that correction rather
than relabeling the original runs as complete sweeps.

Inherited source actually read: R64 Section 2, lines 34-54 of
`proofs/r064-source-state-mass-polytope-variation.md`, SHA-256
`ef14cc7ee5b1384c241eb69bf78b6ae7639f1916afe58723e9187bef3a5cbdfe`,
from root integration at published R70 commit
`0d1272a4af4eb81a77842012dd56db794656f6a0`. This is a reading of the exact
cap proof, not a rerun of R64's optimization or geometry.

## 5. Actual production certificate audit

`review_r071_production_certificates.py` imports no producer module and
does not call an LP solver. Its final actual run used Python 3.15.0,
took 4.988323 seconds, exited zero and had owned PID 11762 absent after
wait. The initial version of this check also passed; the final run adds
the saved single-midpoint dual check. Initial output/receipt are retained
in data temporary storage and are not reported as the final run.

The checker starts from the exact source coordinates. S is the 5-by-5
grid of spacing 6/25 on [0,24/25]^2, with old mass 1/25 at each point.
The original outside list has 144 positions. The six receivers are
(-1/10,6/25), (-1/10,12/25), (-1/10,24/25), (0,53/50),
(12/25,53/50), (18/25,53/50), with fixed CLOSED witness radius 1/4.

For the original six-receiver model I reconstructed all 8056 source
coordinate rectangles having an outside capture, and their 6856
distinct minimal outside-mask caps. Every exported row is checked back
against exact membership in its actual rectangle, its old captured mass,
its rhs and coefficient signs. The 108 residual rows are exactly every
distinct complete distance event at the six receivers. With six witness
rows, all 8170 exported rows are verified. The saved rational primal and
dual are both feasible with the same eta=-1/100; every dual column is
nonnegative and its free-eta column is exactly 1. Thus this is an exact
fixed-grid, fixed-radius strict obstruction, not merely solver status.

Each of the six receivers is individually positive feasible. I checked
the exact original source/candidate D4 symmetries and transformed three
saved single-point primals, then substituted each transformed input into
all actual positive-support caps and every residual distance event.
The six margins are 7/100, 7/100, 151/2500, 151/2500, 7/100, 7/100.
The saved midpoint primal and sparse dual were additionally checked in
its independently reconstructed single-point model and match eta=7/100.
I did not independently solve or certificate-check the six leave-one-out
models, the other pilot menus, or their numerical optimizations.

The original sparse dual does not hold at arbitrary outside positions.
An independent arrangement partition, including unbounded atoms, all
boundary coordinate singletons and all forbidden internal positions,
checked 816 permitted atoms. Its minimum slack is -1/3, matching the
production oracle. The exact reported point (-7/20,23/100) and the
two-dimensional interior point (-3/10,47/200) each have two witness
weights 1/6 and no covering upper weight. The interior point lies outside
the closed collar [-1/25,1]^2, so a zero-cost collar repair does not fix
this missing column. It remains a negative column for open versions of
these witnesses as well. The 236 negative atoms in this independent
partition are not comparable as a raw count with the oracle's 172 atoms
in a coarser partition restricted to the witness union.

Adding exactly q=(-3/10,47/200), and rebuilding ALL caps and residual
events, reverses the conclusion. The independent reconstruction checks
10613 coordinate rectangles, 7816 minimal cap masks, four genuine
zero-cost collar rows, 114 residual prefixes and six witness rows:
7940 exported rows. The saved rational primal and dual match at
eta=101/10000>0. No missing old cap is carried over from the first grid.

The resulting input has EIGHT positive outside sources (145 candidate
positions) and 33 positive atoms in total. From these actual positive
positions alone I separately reconstructed 380 coordinate rectangles,
42 distinct minimal outside-cap masks, and respectively 5,6,7,7,6,6
complete residual distance events at the six receivers. All pass. S
retains maximum eligible mass exactly 1; each receiver is outside both
the old halo and the actual residual halo. Each closed radius-1/4 label
has full mass 2601/10000. Consequently the rational radius
40101/160000, strictly above 1/4 and below 1/2, is a genuine OPEN full
witness for every receiver, with strict threshold verified exactly.
This is a real common outside input, rather than separate mass choices.

All exact input hashes, row counts, source hashes, arithmetic outputs
and execution receipts are in `computation/results/round071-independent`.
This finite shared positive example invalidates any arbitrary-position
interpretation of the original fixed-grid obstruction. Neither six
positive receivers nor the LP objective is an area certificate. There is
no new general deletion budget or bound on C2.

The separately owned `r071-independent-robust-deletion-cells.md` was
read in full. Its displaced-prefix lemma retains every residual radius,
requires epsilon below the first positive distance, and correctly uses
the same full-label radius at every point of each cell. The five closed
cells of half-side 101/40000 and common open full radius 253/1000 have
area 10201/80000000 contained in actual deletion. Its separate kernel
checked all 6*255 residual subsets by the exact cell minimum formula,
all ten pairwise separations and an independent union-area sweep. I
source/code/result reviewed and bound that kernel; I did not repeat its
1530 subset/cell tests. The sixth point has exact residual equality at
distance 7/50 and cannot contribute a positive symmetric residual-free
cell. The proof correctly avoids a purely-mixed-region containment claim.

## 6. Final versions and precise remaining scope

The final author core packet, production packet and continuous/robust
packet are respectively bound at SHA-256
`c7e60eda2c36787889ab43d0fc4269b69e22369bd5089ff6ccb432a9042c3c06`,
`a309ecde3b73e27d8eac130bf676dd10c8e779f6ff23017625b3eff752c6e21b`,
and `e41aeee6f9710a53c3876759d5bbb9e04f73f16f16b9afc8535dbe18d44195ad`.
The final round report is pinned at
`d5e795310c1691c0721a1cc22507d22bec3555f95f2a3fb6568521b633d46a45`.
Presence/byte checks are separate from the proof readings and exact
mathematical checks described above.

The production metadata correction (historical static radius49 scope
text for actual radius25 cases) does not change historical case geometry.
Its old executed source snapshot and correction record are retained;
the new radius25 smoke is separately recorded. This review's mathematical
checks use the per-case radius and exact exported rows, not that stale
scope string. Producer's 66 LP calls include one metadata-repeat call
and are not 66 distinct inputs. Root's three MATLAB numerical solves and
Wolfram exact scalar/QE checks are separately owned evidence; I neither
launched those programs nor present their results as my own runs.

The two saved single-receiver input files `menu.json` and
`single-midpoint.json` are included byte-for-byte in this independent
packet because the portable checker needs them. Their original source
hashes remain unchanged. Other omitted exploratory outputs are not
needed by this checker. No final input relies on an unarchived workspace.

Exact C2, general deletion and Wprime selection, nonterminal-domain
maximizing-block choices and a universal area payment remain unresolved.
No commit or push is made here; root owns publication and remote checks.
The final byte/provenance checker actually passed from an independently
staged archive root: three component manifests, 61 payload/dependency/
external-source rows (including repetitions), both independent result
and receipt chains, and the frozen round report. PID 13452 exited zero
and was absent after wait. This establishes byte bindings and archive
portability only; it is not an additional mathematical computation.

Portable boundary controls from the integrated topic directory:

    python3 -B computation/review_r071_joint_measure_boundaries.py --output /Volumes/data/CodexTemp/T2/r071-replay/boundaries.json

Portable certificate audit, after all component files are integrated:

    python3 -B computation/review_r071_production_certificates.py --model computation/inputs/r071/six-active-exact.json --results computation/results/round071-joint/six-receiver-check.json --menu computation/results/round071-joint/menu.json --single computation/results/round071-joint/single-midpoint.json --augment-model computation/results/round071-joint/one-column-model.json --augment-result computation/results/round071-joint/one-column-result.json --geometry computation/inputs/r071/six-active-dual-geometry.json --oracle computation/results/r071-six-active-continuous.json --output /Volumes/data/CodexTemp/T2/r071-replay/production-certificates.json

Portable byte/provenance binding check:

    python3 -B computation/review_r071_source_bindings.py --packet computation/results/round071-core-integration-manifest.json --packet computation/results/round071-joint/component-manifest.json --packet computation/results/r071-continuous-robust-component-manifest.json --review-root . --archive-root .. --report rounds/071.md --output /Volumes/data/CodexTemp/T2/r071-replay/source-bindings.json
