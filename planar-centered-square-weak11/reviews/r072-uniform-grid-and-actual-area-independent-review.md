# R72 independent review: actual halos and the uniform-grid question

2026-11-16. Owner: pcw-sol-cover-structure. This workspace changes only
R72 independent files; earlier rounds are frozen. Root's R71 publication
9050e6cbe03a94237591608d5539304882f7b737 is the inherited archive.

## 1. Fixed question and admissible evidence

The finite-atomic question fixes S_n, n-by-n equal positive masses 1/n²,
on [0,w_n]², w_n=1-1/n², spacing h_n=(n+1)/n². For a single actual finite
positive atomic outside measure nu retaining maximum eligible mass 1,
the object is Delta=H(S_n+nu) minus H(nu), with all centered common radii
and strict threshold 1. The target is a uniform Delta area bound 7/2 or a
strict actual counterexample. R70's source-only envelopes tending to 4 do
not settle this question. A Borel version permits diffuse nu and is a
separate question unless an explicit halo/deletion transfer is proved.

R71's whole-cell interface was reread from the published archive:
`proofs/r071-shared-outside-mass-and-residual-duality.md`, Section 9.
For a closed receiver rectangle C and every residual bbox closure T,
nu(T)<=4 delta_C(T)² is exactly the exclusion of the entire residual
halo from C, including all radii and tied/zero distances. A chosen full
label satisfying mass(U)>4 R_C(U)² is sufficient for all C to be in the
full halo. It is not necessary for coverage by a union of different
labels. Infeasibility using one full label per C cannot exclude that
more general coverage. A positive collection of fixed-Q cells is an
area lower certificate for one shared input, not a universal upper bound.

Zero-mass candidate locations are omitted from the actual bbox and
support. Every maximum eligible block tie must be evaluated if the
conclusion concerns maximum-block selection. A closed radius with a
strict mass margin yields a slightly larger true open witness; equal
mass/area alone does not. Finite candidate positions and fixed full
radius menus do not cover arbitrary outside positions or other radii.

## 2. Independent analytic internal-halo localization

Every nonempty bbox closure in S_n is a consecutive k-by-l grid rectangle.
It has mass kl/n² and largest coordinate span h_n max(k-1,l-1). Its halo
receiver is the open rectangle obtained by expanding inward from its
opposite bbox endpoints by sqrt(kl)/(2n). Completing any source subset
to its bbox closure preserves its bbox and only increases mass, so these
closures give the complete internal halo I_n=H(S_n).

Exactly the following rectangle shapes are eligible:

* k=l, always;
* |k-l|=1, with min(k,l)(2n+1)<n².

Indeed, if k>=l+2 then kl<=(k-1)²-1 and h_n>1/n rules out eligibility.
For squares, k/n>h_n(k-1) is equivalent to k<n+1. For the adjacent shapes
k-by-(k+1), direct squaring reduces to k(2n+1)<n². Equality is excluded.

For a receiver extending to the left of the grid, its greatest possible
extension from a k-column label is

    sqrt(kl)/(2n)-(k-1)h_n.

For an eligible shape l<=k+1. At k=1 its numerator sqrt(kl) is at most
sqrt(2). For k>=2, sqrt(k(k+1))<k+1/2<=2(k-1)+sqrt(2). Together with
h_n>1/n, this proves every extension is at most sqrt(2)/(2n). Reflections
give the other three faces. Therefore, for every n>=2,

    I_n subset [-sqrt(2)/(2n), w_n+sqrt(2)/(2n)]²,
    |I_n| <=(1-1/n²+sqrt(2)/n)².                       (1)

The n² singleton halos are disjoint squares of side 1/n because their
spacing h_n is strictly larger than 1/n. Their total area is exactly 1.
Consequently

    1<=|I_n|<=(1-1/n²+sqrt(2)/n)² ->1.                 (2)

For n>=3 a 1-by-2 label is eligible and achieves the one-face extension
sqrt(2)/(2n); at n=2 the largest extension is 1/(2n). Equation (2) concerns
the actual INTERNAL halo, not mixed receivers or the general Delta.

As a scope control, n=2 itself is already paid by the simple source
carrier: |Delta|<=|bbox(S_2)+[-1/2,1/2]²|=49/16<7/2 for every admissible
outside atomic measure. At n=3 that support carrier exceeds 7/2 by only
11/162; its excess tends to 1/2 as n grows. Filling this gap requires
actual shared outside geometry and residual credit.

The new Fraction checker `review_r072_uniform_grid_internal.py` actually
tested all 22379 ordered shapes for n=2,3,4,5,9,17,33,65,129. It checked
the exact strict classification, isolated-radical extension inequality,
and diffuse cap scalar identities below. Directed 96-bit root endpoints
are retained for the containing-square bound, not labelled exact halo
areas. PID 17737 exited zero and was absent after wait. The finite
controls support, but do not replace, the all-n analytic proof above.

## 3. Diffuse outside-family cap review and transfer boundary

Root proposed the finite Borel family

    J=[-1/n²,1]², t=n/(n+1), B=[-1,2]²,
    nu=t 1_(B minus J) dx.

I independently checked its maximum-mass rectangle-cap argument. For
a one-dimensional interval [x,x+R] of length R>=1 capturing k>=1 grid
coordinates, let a_k=(k-1)h_n+1/n². Its intersection length with
[-1/n²,1] is at least a_k: the overlap is the nonnegative minimum of
R, 1+1/n², x+R+1/n² and 1-x. The first two are at least 1>=a_k;
capture implies x+R>=(k-1)h_n and x<=(n-k)h_n, giving the last two
bounds. Closed endpoint captures cause no problem for Lebesgue overlap.

For a square of side R>=1 capturing k-by-l old points, with k,l>=1,
its total mass is at most tR²+f(k,l), where

    f(k,l)=kl/n²-t a_k a_l,
    f(k+1,l)-f(k,l)=(1-l/n)/n²>=0,
    a_n=1, f(n,n)=1-t.

Thus mass<=tR²+1-t<=R². If k or l is zero the direct tR² bound suffices.
For R<=1, containment in a closed unit square yields mass<=1; it does
not imply mass<=R². These are exactly the Borel admissibility caps
mass(Q)<=max(1,R²). Truncating the outside background to finite B only
decreases mass. Its density t<1 gives H(nu)=empty at every radius.

This family is not finite atomic outside input. Its admissibility and
empty residual halo cannot be transferred by unqualified weak convergence
or grid discretization. Equation (2) gives a direct warning: the atomic
S_n converge weakly to unit density on [0,1]², whose strict halo is empty,
while their own halo areas tend to 1. More generally a grid of N² atoms
of total mass t<1 on [0,1]² has disjoint singleton halos of total area t,
even though its weak limit has density t and empty halo. A finite-atomic
counterexample requires its actual residual halo to be subtracted.

For root's Borel area helper, inside C_n=(-1/2,w_n+1/2)² every candidate
r<1/2 square lies in B. The exact excess is

    old_closed_prefix_mass
      -t |J intersection Q(x,r)|-4(1-t)r².             (3)

Between old distance events the prefix mass is constant and the two
cost terms are nondecreasing. Complete closed old-distance events,
including zero and all ties, therefore suffice at a fixed receiver.
Strict positivity permits a slightly larger open witness. For a receiver
cell, possible old-source counts plus the minimum J overlap give a valid
upper envelope; the count can change only at source-to-cell minimum
distance events. Guaranteed old-source counts plus maximum overlap at
a FIXED radius give a sufficient whole-cell lower witness. In one axis
overlap is maximized by projecting the interval midpoint onto the receiver
interval, and minimized at one of its endpoints. These bounds must remain
within r<1/2 and C_n; no finite point sampling certifies area.

The Borel helper code was independently read and these classification
lemmas checked in full. The saved certificates were then independently
replayed without importing the adaptive classifier or repeating its search.
For n=5 all 198139 leaves form a complete dyadic partition. Every 30669
excluded and 38747 inner leaf was checked by literal TWO-dimensional
source distances, rather than the producer's separately sorted 1D lists.
That executed 1735400 atom-distance tests and 54705 actual 2D possible
distance events. All exact areas and margins matched, independently
certifying the entire bracket

    27886815467/10485760000 <= |Delta_Borel,n5|
      <=28096528411/10485760000.

PID 20942 exited zero and was absent after wait; internal time 1.347031s.
For n=17 I checked the complete 306433-leaf partition and exact area
bookkeeping, plus 16 deterministic quantile leaves of EACH classified
state. PID 21181 exited zero and was absent; time 1.427957s. I did NOT
independently replay all n17 classifications. Its reported bracket
approximately [2.67133,2.73836] remains producer evidence with proof/code
review and sampled labels; width approximately .067 did not reach .02
before the stated 200000 processed-cell budget. Unresolved cells remain
in the upper area. No adaptive search or further refinement was run here.

## 4. Actual new complete baseline area audits

The standalone `review_r072_actual_halo_y_slabs.py` imports no producer
geometry. It starts from exact positive atoms, enumerates all coordinate
rectangles by direct point-membership scans, tightens and deduplicates
closures, tests strict eligibility and every maximum tie, and recomputes
directed integer square-root bounds at the producer's 72-bit scale.
It then uses horizontal y-slabs, literal active-rectangle scans,+x-interval merging and DIRECT interval difference for area. Production
uses vertical x-events and a three-channel segment tree. These are
different area algorithms and different closure implementations.

Two new actual runs matched every full closure catalog entry, every
eligible directed rectangle, every maximum block and all exact rational
lower/upper halo and deletion area endpoints:

| Input | Positive atoms | Maximum ties | Full halo (decimal display) | Actual Delta for all maxima (decimal displays) |
| --- | ---: | ---: | ---: | --- |
| R71 actual input, newly evaluated | 33 | 2 | 1.9719531518213766 | 1.4419810899674625; 1.3941316996330793 |
| S25 only, newly evaluated | 25 | 1 | 1.373993652089341 | 1.373993652089341 |
| Aligned225 exact positive input | 38 | 4 | 3.3050127813748214 | 1.5551882180744014; 1.3390731136614293; 1.3109816138066264; 1.2963772190195557 |
| Same input, outside masses contracted | 38 | 1 | 3.305012762967794 | 1.5551882171496079 |

The residual and internal halos for both maximum blocks were independently
recomputed as well. The R71 actual input remains terminal, with exact
span²-(35/33)totalmass=-371/15000. It is not evidence in the remaining
nonterminal domain. Its input provenance is from R71; the complete halo
and area computations above were newly executed in R72, not inherited
old execution. PIDs 16656 and 16657 exited zero and were absent after
wait. Exact areas, software versions, input/source hashes and actual
receipts are retained in `computation/results/round072-independent`.
The two new 38-atom candidates were also fully rechecked from atoms with
the same independent horizontal algorithm, including all four original
maximum ties and the contracted unique maximum. Both are connected and
strictly nonterminal. Every catalog, directed rectangle and exact area
endpoint matched. PIDs 18686 and 20019 exited zero and were absent after
wait. The large full halo approximately 3.305 is not the deletion area:
the complete residual halo must be subtracted.
The contracted unique maximum's exact lower endpoint exceeds 3/2 and
the R71 chosen block's exact upper endpoint. Thus this actual connected
nonterminal input disproves universal Delta<=3L/2 payment in that domain,
including when the maximum S_5 is unique. This does not improve the known
general C2 lower bound and does not threaten the 3.5 target.

## 5. Independent whole-cell positive certificate

For the aligned225 input, I separately read the original candidate
positions, rational weights and full-label certificate; no MIP optimizer
was imported or rerun. The actual residual has THIRTEEN positive atoms.
Direct geometric closure reconstruction gives 120 actual residual
closures. For each of the 89 selected closed cells, I checked all 120
exact delta_C inequalities, 10680 actual source/cell comparisons. I also
checked every chosen full label against its farthest distance on the
entire closed cell, retaining zero-weight candidate positions when they
were included in that label. Every original full margin is at least
1/1000000. Thus this is a whole-cell certificate with all residual radii,
not a sample of receiver points.

The complete 225-cell partition passed 25200 pairwise disjoint-interior
tests and an exact area-union check. Horizontal merging independently
gave the selected union area

    |union_(89 selected C) C|=345909/250000=1.383636,

contained in actual deletion. Contracting outside masses by
99999999/100000000 preserves all selected residual constraints and leaves
minimum full margin 19907871/20000000000000>0. The SAME 89 cells therefore
certify the same area for the unique-max contracted input. This exact run
exited zero, PID 20516 absent after wait. It certifies feasible lower area,
not a rigorous MIP optimum or upper bound, even on the fixed candidate
family. A complete critical R_C menu is exact for one full label covering
C; union coverage by several full labels can still cover cells that this
menu formulation does not certify.

## 6. General reductions: full independent mathematical reading

The two new author proofs were independently read in full, including
their boundary and quantifier qualifications:

* `r072-uniform-grid-borel-family-and-atomic-discontinuity.md`;
* `r072-maximizer-tie-breaking-reduction.md`.

The constant density t is sharp within its specified B minus J family:
the unit square [-a,1-a] times [0,1], a=(n+2)/n², has right endpoint
(n-2)h, captures exactly n-1 closed old columns and n rows, and has
outside area h. Its cap forces c<=t, including n=2. This is not optimality
among arbitrary outside measures.

For canonical atomization, take ALL complete delta grid cells in B minus
J and put mass t delta² at each center. A square of side R captures
atomic mass no larger than the original Borel measure in its concentric
side-(R+2delta) expansion. All positions lie in B of span 3. For an
eligible subset of mass m and span R>1, the affine Borel cap yields

    m<=t(R+2delta)²+1-t
      <t m+1-t+t(12delta+4delta²),

so m<1+n(12delta+4delta²). For R<=1, the ordinary cap gives
m<=(1+2delta)², bounded by the same expression. Since S_n stays eligible
with mass 1, the actual global maximum L_delta tends to 1 WITHOUT an
assumption that S_n remains exactly maximizing at each delta.

The singleton outside halos have side sqrt(t)delta, fit strictly inside
their own cells and are mutually disjoint. Counting only cells fully in
C_n minus J gives residual-halo area liminf at least t|C_n minus J|.
Every deletion witness contains an old source and has radius below
sqrt(L_delta)/2, so deletion is contained in the enlarged C_delta tending
to C_n. This proves the canonical sequence bound

    limsup |Delta_atomic,delta(S_n)|
      <=U_n=(1+1/n²)²+(3-6/n²)/(n+1),
    U_5=1927/1250=1.5416.

The n5 independently certified Borel lower area therefore exceeds this
atomic limsup bound by at least

    11721967851/10485760000 >1.11789.

This is a strict deletion-area discontinuity for THIS canonical transfer,
even with L_delta tending to 1. It is not a bound for arbitrary atomic
outside measures and not a finite-atomic counterexample to 3.5. Additional
outside contraction factors tending to 1 can only improve the upper
mass bound and retain the same limiting singleton-halo lower density.

For the maximum-tie reduction, keep any specified maximum S unchanged
and multiply all outside masses by 0<lambda<1. No ineligible label is
activated. Any other eligible label either is a proper positive subset
of S, or has decreased outside mass, and has mass strictly below L.
Thus S becomes uniquely maximizing. As lambda increases to 1, the full
and residual strict halos separately increase to their original halos;
each fixed positive witness margin persists. Their deletion indicators
converge, and all halos lie in one bounded set for finite sources. Area
convergence follows even though deletion sets need not be monotone.
Finite eligible margins preserve connectivity for lambda close to 1;
unchanged support and decreased total mass preserve strict nonterminality.

It follows that the UNIVERSAL statement "every input in a perturbation-
stable class has some maximum S with Delta<=K L" is equivalent to the
corresponding statement for EVERY maximum S. This is not a pointwise
equivalence on one unchanged input. The sign <=K matters: an equality
limit does not establish a strict <K version. Extra domain conditions
must themselves survive the perturbation. In particular a strict
violation at one uniform maximum grid could be uniquely maximized, but
this reduction does not produce such a violation.

## 7. Source versions, computation ownership and reproduction

Final author proof versions reviewed are respectively
`7b49ddee130e5deee6bbf17641554fbfe122dae9925f73c865574e7073eddf53`
(Borel/internal/atomization) and
`673bf41ad02a5656e123f28db369e9a9a2fbcca63855cf0831678e16ab6fa981`
(tie reduction). The author scalar code, result and receipt were read:
4552 cap pairs, 4460 forward differences, 1119 selected interval overlaps
and six exact canonical-cell counts are author-run finite controls, not
my independent rerun and not an all-n proof. The analytic arguments
above and in those proofs carry the general conclusions.

The core, cell and complete-geometry/Borel component manifests are frozen
at SHA-256
`daf016d0f89634bf91a7553fffa7bcd5a1aa69ad69c710f453b3a8f198e47025`,
`a3f550f6a2a6afde8f5743382eeedf051c110ba4d631b7e24a97f5bfd39088f5`,
and `9953334275b4a5a69609ade52031afa5aa21b629136aad2470eb28a34532276b`.
The final round report is pinned at
`ee6f1849cbb489c25fbc31ae32345b0226461b8f36d4e5dc4ffc0e28224720db`.
Their presence/byte binding is distinct from proof reading and the eight
actual independent mathematical executions described here. Root owns
MATLAB/Wolfram evidence and publication; no such process was launched by
this workspace. The producer's 18 MIP solves were not independently
repeated, and the coarse144 exploratory family was not independently
area rechecked. No floating MIP objective is an exact upper certificate.
The final archive-only byte/provenance checker passed three component
manifests, 71 payload/dependency/source rows (including repetitions), the
eight independent mathematical result/receipt chains and the frozen
report. It also exactly compared the saved unique-input area lower
endpoint against 3/2 and the R71 chosen-input upper endpoint. PID 23027
exited zero and was absent after wait. This step did not rerun any
mathematical computation, LP or adaptive search.

Portable complete-geometry audit from the integrated topic directory:

    python3 -B computation/review_r072_actual_halo_y_slabs.py --input computation/inputs/r072/r071-positive33.json --saved computation/results/round072-actual/r071-positive33.json --producer computation/round072/actual_halo_geometry.py --bits 72 --output /Volumes/data/CodexTemp/T2/r072-replay/positive33.json

Use `s25-only.json` in the same two input/result directories for that control.
For the two new38-atom cases use inputs `aligned225-actual-input.json`
and `contracted-unique-actual-input.json`, with saved outputs `aligned225.json`
and `contracted-unique.json` in `computation/results/round072-actual`.
The k,l scalar controls are independently replayable:

    python3 -B computation/review_r072_uniform_grid_internal.py --input computation/inputs/r072-uniform-grid-internal.json --output /Volumes/data/CodexTemp/T2/r072-replay/grid-internal.json

Whole-cell certificate, including its contracted carryover:

    python3 -B computation/review_r072_whole_cell_certificate.py --model computation/results/round072-cells/aligned225-model-exact.json --result computation/results/round072-cells/aligned225.json --actual-input computation/results/round072-cells/aligned225-actual-input.json --contracted-input computation/inputs/r072/contracted-unique-actual-input.json --output /Volumes/data/CodexTemp/T2/r072-replay/whole-cells.json

The n5 complete saved-leaf audit:

    python3 -B computation/review_r072_borel_cells.py --summary computation/results/round072-borel/n5.json --certificate computation/results/round072-borel/n5.cells.json.gz --producer computation/round072/borel_grid_adaptive.py --mode all --output /Volumes/data/CodexTemp/T2/r072-replay/borel-n5.json

For the actual n17 review scope, change both n5 filenames to n17 and use
`--mode sample --sample-per-state 16`; `--mode all` would be a NEW broader
execution, not the review recorded above.

Archive-only final source/receipt check:

    python3 -B computation/review_r072_source_bindings.py --packet computation/results/round072-core-integration-manifest.json --packet computation/results/round072-cells/component-manifest.json --packet computation/results/r072-geometry-component-manifest.json --archive-root .. --report rounds/072.md --output /Volumes/data/CodexTemp/T2/r072-replay/source-bindings.json

No new general finite-atomic deletion or C2 bound has been established by
these computations or reductions. The remaining obligation is an actual
arbitrary-outside atomic budget or a certified strict atomic violation.
This workspace performs no commit or push; root owns publication.
