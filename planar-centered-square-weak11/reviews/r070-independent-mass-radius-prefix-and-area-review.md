# R70 independent review: mass, capture radius, prefixes and area certificates

Original review/runs 2026-11-15; recovered provenance and finite certificate
finalized 2026-11-16. Owner: pcw-sol-cover-structure. Prior rounds are frozen. This
review changes no producer or author file and launches no MATLAB process.
The accompanying manifest fixes the final source versions and actual runs.

## Mathematical conclusions

The necessary mixed-label inequality is correct, with a strict sign:
for maximum eligible S of mass L, mixed T of internal mass a>0 and outside
mass b>0, total m=a+b<=L, every strict receiver x satisfies

    c_S(x)>sqrt(L+b)-sqrt(a+b)/2.

The independent proof uses S union T of mass L+b>L. Within-S and within-T
spans are below sqrt(L+b); if the reverse weak inequality held, strict
receiver capture makes every cross span strictly below it as well. Thus
the union would be eligible. No bbox closure of T, eligibility of U=T∩S,
or independent source placement for each pair is needed.

For fixed a, f_a(m)=sqrt(L+m-a)-sqrt(m)/2 has derivative sign
3m-(L-a). For 0<a<L/4 its minimum on a<m<=L occurs at (L-a)/3. At
a=L/4 and above, b=0 is excluded: the value at m=a is an infimum, not an
attained minimum. In both cases the necessary inequality remains strict.
At c/s=3/4 the inverse gives a>L/4; at c/s=1/2 no mixed receiver exists;
at c/s=sqrt(3)/2 the inverse adds only a>0. The two inverse branches and
the L-dimensional scaling are correct.

For an internal U, the exact scalar condition uses
q=max(a,4r_U^2,(L-a)/3)<L and c_S>f_a(q). The strict q<L condition is
essential; q=L cannot satisfy the strict capture radius. If q is an
excluded endpoint below L, the strict gap permits a small upward mass
perturbation. The author's rational sign split isolates the nonnegative
radical before squaring, so it preserves strictness and equality cases.
Outside I=Halo(S), the entry a in the clamp can indeed be omitted because
otherwise x would lie in the internal receiver E_U(a).

The all-H strengthening is also correct: for every nonempty H subset S\U,
mass M with M+m>L, the same union argument for H union T gives
c_H>sqrt(M+m)-sqrt(m)/2. H need not itself be eligible. Equality in
M+m=L supplies no constraint from maximality. At fixed x,U,m, checking
the complete distance prefixes including all ties is equivalent to all
H: enlarging H to its complete prefix preserves c_H and increases M.
The same m must satisfy every prefix. A zero-distance prefix is retained.
The active prefix inequality is exactly
M<(R-r)(R+3r), r=sqrt(m)/2, with the alternative M+m<=L kept.

The author's inclusion W subset I union M_prefix subset I union M uses
actual closed mixed labels. Their internal U is closed within S. A code
enumerating all proper internal bbox closures must include internally
ineligible U. The claimed outer inclusion needs only actual closed mixed
labels, whose U is closed within S. For an arbitrary proposed U, bbox
closure preserves r_U. If its mass reaches the proposed m, the strict
capture radius already puts x in I; otherwise closure has mass below m,
decreases f_a(m), and leaves fewer H constraints. Thus this argument is
for the union with I, not an unqualified pointwise equivalence of the
mixed carrier alone. Raw all-subset formulations also contain the needed
actual captures. No one of these necessary
conditions asserts that an outside source state exists.

## Positive balanced four-corner theorem

The kernel proof covers arbitrary finite outside sources, provided the
equal-quarter square S retains maximum mass and its normalized side t is
in [sqrt(3)/2,9/10]. I read the full proof and the exact-checker source.
Its stronger union J of all eligible touching-label outer receivers
satisfies W subset J, and the captured-U carrier is valid for every such
label. For this t range the only pure eligible internal labels are the
four singletons and S; triples at the lower endpoint are excluded by the
strict threshold. The exact coarse-carrier union decomposes into:

* the singleton union C\K_(3/4), area 5t-5/4;
* five disjoint internal squares inside K_(3/4), area 5(1-t)^2;
* four disjoint new adjacent-pair strips, area (1-t)(sqrt(2)-1).

Diagonal-pair and triple receivers are contained in the internal halo.
The strip locations, mutual disjointness, endpoint direction and scale
give the stated exact formula
F(t)=15/4-5t+5t^2+(sqrt(2)-1)(1-t).
It increases over the interval, so F(t)<67/20<7/2. This is a general
sufficient class, not a three-fixture inference. The author's three
fixed-geometry exact runs and receipt are source-bound; I did not rerun
that checker or an independent complete kernel-fixture area calculation.

## Source-only obstructions and the stronger limiting result

The author's four-corner d=99/100 obstruction is valid. The four side-one
singleton carriers cover C up to its boundary. For any x in C\K_rho,
rho=(sqrt(7)-1)/2, choose a capturing singleton with r_U<1/2. Its clamp
lies below 1 and f is increasing from its minimum at mass 1/4, so the
strict c_S>rho implies membership. The area is
(199/50)sqrt(7)-7>7/2; the positive square-comparison margin is 791/1250.
The five-source ambient embedding is unique-max, connected, fat and
nonterminal by the stated strict rational margins. Its added source at
height d/2 connects through span d/2; its actual halo was not the S-only
production control. I read its 31-subset author checker and exact rows;
I did not rerun that small author checker.

The independent proof
`proofs/r070-grid-clamped-envelope-near-four.md` proves more: for positive
uniform grids S_n, both the full internal-capture clamp and its common-m
all-prefix strengthening have area tending to 4. It supplies a legal
single-outside-source ambient embedding by a general column-count bound,
not an inference from eligibility of S_n alone. That bound rules out
EVERY q-containing subset of mass at least 1. It also proves connectivity,
nonterminality, fatness and failure of the earlier support-tree criterion.
It makes no assertion that all heavy-core certificates fail.

The all-prefix extension is analytic: uniform open/closed rectangle
discrepancy, a moving common radius r_n=sqrt(1-n^(-1/2))/2, and the active
prefix condition produce one strict margin simultaneously for all
remaining prefixes at each fixed receiver outside a null set. The moving
radius overlap formula explicitly requires r_n>|delta_i-1/2|, which
eventually holds almost everywhere. Fatou's lemma then gives area 4.
This limiting argument alone is not a finite-grid area certificate for
M_prefix. The six finite exact areas below concern a fixed-total
sub-envelope of M only. The separately added proof section 6 does give
one explicit finite M_prefix/common-m area certificate: n=400001,
m=99/100, and Omega=((-49/100,49/100) union (51/100,149/100))^2.
I checked each analytic bound, including Omega subset C_n. Its area is
2401/625, the capture mass is strictly between 0 and m, and every active
prefix has cap-minus-mass greater than
1/20000-402/(25n)=78401/8000020000>0. The same m works for all H.
Root separately read and checked this complete analytic derivation.
The input's 160000800001 sources are specified by formula; none were
enumerated by this certificate. It makes no actual halo claim.

## Actual independent finite computations

All three entry points use the Python standard library, with Fraction
and integer arithmetic. No producer or previous geometry module is imported.
The three original run results record Python 3.13.18. Their original
contemporaneous execution reports stated exit zero, but no separate
persisted PID/exit receipts were found on recovery. On 2026-11-16 I
checked the existing output bytes, their source/input hashes, and result
contents; I did not repeat these runs or re-observe their old exit states.
The recovery-provenance JSON explicitly records this receipt gap.

1. `review_r070_grid_clamp_obstruction.py`: six odd grids n=3,5,9,17,33,65.
   It partitions at every source coordinate plus/minus 1/2, counts exact
   mass A, and subtracts the central square for cutoff sqrt(2-A)-1/2.
   Integer-isqrt at 80 bits bounds the radical outward. The resulting
   area lower bounds are 3.242330569895, 3.707650300417, 3.874437996709,
   3.952014089891, 3.981712671134, 3.993218376951 (decimal displays only).
   Exact rational endpoints and exact ambient margins are saved. It does
   not enumerate 2^(n^2+1) subsets; section 3 of the proof handles them.

2. `review_r070_saved_cells.py`: all three producer controls were replayed
   from a byte-staged archive-relative fixture. Literal source-subset
   enumeration checked 255,15,65535 subsets; all 8,8,72 proper internal
   bbox closures were reconstructed independently. All 12193,12148,12073
   saved quadtree leaves were proved to partition their exact domains,
   including unknown leaves. All classified leaves passed independent
   three-branch radical elimination, with 166551 label tests total.
   Unknown leaves remain exclusively outward. Coordinate-compressed 2D
   difference arrays checked 2745792 area cells and matched the saved
   exact endpoints for M, I, the additive inner/outer I-union-M budget,
   M-minus-(I-union-R), and balanced8 true Delta. Historical true P/Delta
   were only source-bound to their R69 output, not recomputed. Uniform4
   has no ambient outside input. Balanced8 actual Delta was independently
   rebuilt; its producer's true P was source/code reviewed but was not
   separately matched by this replay. Auxiliary M-minus-I and old-credit
   output fields were source/code reviewed, not independently area-matched.

3. `review_r070_grid_prefix_points.py`: only n=17 was independently
   replayed. It lists all 289 old sources and groups their actual 2D
   L-infinity distances at all 256 saved receivers. Every complete row,
   including captures, active prefixes, strict slack, violations and
   first violation, matches root's event-count implementation. There are
   5186 actual grouped distance events and 179 passing centers. Root's
   larger n=65,257,1025 runs were source/code/receipt bound only. Their
   188,255,256 passing centers are root's results, not independent replays.
   These finite-point checks use m=(1-1/floor(sqrt(n)))^2, which differs
   from the common m_n used in the analytic limiting proof. No finite
   sample fraction is interpreted as an area certificate.

4. `review_r070_finite_prefix_area.py` was newly run on 2026-11-16 using
   Python 3.15.0 and Fraction arithmetic. Its actual saved result checks
   every rational margin and the positive-radical square comparison in
   proof section 6. The process exited zero and its owned PID was absent
   after wait, with a separate persisted receipt. This scalar calculation
   complements the analytic rectangle-discrepancy proof; it does not
   enumerate the formula-specified grid, prefixes or source subsets.

The saved results preserve exact source/input/output hashes. The source
binding addendum also checks all producer, author, kernel and root packet
payload hashes and the root receipts. Presence and matching hashes are
reported separately from the mathematical reading and actual replays.
Root's 14 Wolfram scalar checks were read and bound, not rerun here.
The final source-binding check pins the captured-mass author proof at
SHA-256 43054eb9d7ec35b8bf8cd27067fb2360ca4ba13c2de012b7bdb8cc890030989d,
the four-corner obstruction proof at
94226fade8d92d444318d3e76e4c24d6e46b5465279af10546b58275bc6f50c8,
and the frozen round report at
c37e6fe36ee153c5b3d7a0f754135d654f3b416abc4fc849c32a1a3c22d328fd.
These hashes identify the versions reviewed; byte identity by itself is
not a proof. The final root prefix-grid manifest is
05e780669ce7e5cab88e5d61e9d4d18e7efbd155d0fe5f818eff67f71b0c7d6b.
Its receipt checks the parent and eight workers that returned results;
twelve were requested, and the process pool was closed. No claim that
all twelve worker PIDs were recorded is made.

## Portable reproduction

Run from the topic directory, using an explicit fresh dataTemp output:

    python3 -B computation/review_r070_grid_clamp_obstruction.py --input computation/inputs/r070-independent-grid-clamp.json --output /Volumes/data/CodexTemp/T2/r070-replay/grid.json

For the saved-leaf replay, `--archive-root` denotes the REPOSITORY root
containing `planar-centered-square-weak11`, with the compute integration
manifest and its payload/dependencies at their declared relative paths:

    python3 -B computation/review_r070_saved_cells.py --archive-root .. --output /Volumes/data/CodexTemp/T2/r070-replay/leaves.json

The small prefix replay uses only explicit supplied archive files:

    python3 -B computation/review_r070_grid_prefix_points.py --input computation/round070/r070-grid-prefix-input.json --saved computation/results/round070-root-prefix-grid/r070-grid-prefix-result.json --producer computation/round070/check_grid_prefix.py --output /Volumes/data/CodexTemp/T2/r070-replay/prefix-n17.json --n 17

The finite analytic area's scalar certificate is independently replayable:

    python3 -B computation/review_r070_finite_prefix_area.py --input computation/inputs/r070-finite-prefix-area.json --output /Volumes/data/CodexTemp/T2/r070-replay/finite-prefix-area.json

These entries have no fallback to `.local`, production workspaces or
unarchived runtimes. The original production sources remain unchanged.
No general W, Wprime, deletion selection or C2 conclusion is established
by the necessary-envelope obstructions. Publication and push verification
belong to root.
