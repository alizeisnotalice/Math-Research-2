# R72: complete actual halo and deletion enclosures

2026-11-16. Finite positive atomic input, threshold 1, R² Lebesgue area,
centered axis-parallel squares with a common radius. These computations concern
the true full and residual halos, with no receiver window or radius truncation.
They do not evaluate the R70 virtual outside envelopes and do not prove a
general 3.5L bound.

## Complete source labels and every maximum block

For a finite source subset T of mass m(T) and bbox [a,b]x[c,d], its exact open
receiver rectangle is

    E_T=(b-sqrt(m(T))/2,a+sqrt(m(T))/2)
           x(d-sqrt(m(T))/2,c+sqrt(m(T))/2).

It is nonempty precisely when max(b-a,d-c)²<m(T). The full halo equals the
union of these rectangles over all nonempty subsets. Indeed, a strict square
mass witness with captured source set T has radius less than sqrt(m(T))/2,
and its center lies in E_T. Conversely a point of E_T has farthest distance
strictly below sqrt(m(T))/2 and therefore admits a strict common-radius square
witness. Boundaries affect atomic source membership in the argument; the
resulting receiver rectangles are open, and their finite boundaries have zero
Lebesgue area.

Complete T to all sources in its closed bbox. This preserves its bbox and
increases its mass, so E_T is included in the completed rectangle. Consequently
it is sufficient, and exact, to enumerate all nonempty tight bounding-box
closures, including their strict eligibility test. This is not enumeration
of only the closures internally eligible for a preselected S.

The production code constructs every coordinate-interval source mask on each
axis and intersects every x-mask with every y-mask. Deduplicating these
nonempty masks gives exactly all bbox closures. Positivity also means that
every maximum-mass eligible subset is bbox closed: an omitted source in the
same bbox would increase its eligible mass. Thus maximizing over the complete
closure catalog enumerates **every** maximum block, including ties. All are
evaluated; the code does not select the worst or a convenient block.

## Directed rational roots and complete areas

Let B be a common denominator of all source coordinates and masses. Scale
coordinates by B and masses by B² to obtain integers. For precision b=72 put
D=B*2^b. For each eligible closure with integer mass M, calculate

    k=floor(sqrt(M*2^(2b))/2),
    4k² <= M*2^(2b) < 4(k+1)².

Both inequalities are checked with integers. The exact physical radius lies
between k/D and (k+1)/D. For an exact square root both endpoints are k/D.
Substituting these radii in the receiver rectangle gives rational inner and
outer rectangles, hence finite unions H^- subset H subset H^+ up to harmless
receiver boundaries. No floating point enters these enclosures.

For each maximum block S the same construction is performed on the complete
residual input nu=mu|_(A\S), and on S itself. True deletion satisfies

    H_full^- \ H_residual^+ subset Delta
       subset H_full^+ \ H_residual^-.

The code evaluates both bounding differences directly, so it does not assume
that its individually rounded inner and outer unions retain nesting.

The area engine is a new vertical x-event sweep. It maintains three independent
interval-union cover counts on a common y segment tree: the left rectangles,
the right rectangles, and their combined union. On each x slab the left-minus-
right length is length(combined union)-length(right union). Integer slab widths
and lengths give exact areas in units D^-2. At equal x all entering/leaving
events are processed before the next slab. This differs from the independent
reviewer's horizontal slab and merged-interval algorithm.

The closure enumeration and root principle follow the previously used
`computation/round064/mass_polytope_variation.py` (SHA-256
`4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246`),
specifically its `closures` and `exact_rectangles` methods. The R72 implementation
is self-contained and imports no old geometry module. Its output pins this
method source and saves the complete closure catalog and every directed
rectangle. No completeness claim is inferred from an old subdirectory name.

## Initial actual baselines

`computation/inputs/r072/r071-positive33.json` is a byte-identical copy of the
published R71 actual input (SHA-256
`4a2708f364a6198fbd774861eb2259441840461918fa961d6c3bd0c856677169`).
`s25-only.json` restricts it to its 25 declared S atoms, preserving all positions
and masses and recording the original hash.

The actual R71 input has 456 bbox closures, 177 strictly eligible closures,
maximum eligible mass 1, and two maximum masks: 33554431 (the original S25)
and 4059003870. Its eligible-source graph has one component. Its actual bbox
span squared minus (35/33) times total mass equals -371/15000, so it is a
terminal input for that criterion. It must not be presented as a nonterminal
counterexample or candidate. The S25-only control is also terminal.

Initial production enclosures, shown below only as decimal displays of the
exact rational endpoints saved in JSON, are:

| Input / chosen maximum block | H(full) area | Actual deletion area |
| --- | ---: | ---: |
| R71 33-positive input, mask 33554431 | 1.9719531518213766 | 1.4419810899674625 |
| Same input, mask 4059003870 | same H(full) | 1.3941316996330793 |
| S25-only, its unique maximum | 1.373993652089341 | 1.373993652089341 |

All three intervals have width less than 10^-21. These are complete finite-
input area enclosures; the short decimal displays are not independent proofs
or claimed exact constants. The candidate outcomes and final independent
review are recorded separately in the integrated round report.

## One finite candidate from the shared-cell optimizer

`computation/inputs/r072/aligned225-actual-input.json` is the byte-identical
38-positive-atom candidate supplied by the shared-cell worker: the same S25
and 13 positive outside atoms. The complete geometry calculation enumerates
816 bbox closures, 193 eligible closures, one eligible-source component, and
four maximum eligible blocks of mass 1. Its actual span is 48/25 and total
mass 2685811/1000000. The nonterminal margin is exactly 5529563/6600000>0,
so unlike the R71 baseline this candidate lies in the nonterminal domain.

The full halo has certified area about 3.3050127813748214. All four maximum
blocks are evaluated in `computation/results/round072-actual/aligned225.json`:

| Maximum mask | Actual deletion area, decimal display |
| --- | ---: |
| 33554431, original S25 | 1.5551882180744014 |
| 120292637696 | 1.3390731136614293 |
| 131012738048 | 1.3109816138066264 |
| 244845606912 | 1.2963772190195557 |

The chosen S25 deletion exceeds the R71 S25 deletion, but the smallest deletion
among this input's maximum blocks is smaller than both R71 baseline choices.
Thus comparing only the originally marked S25 would misrepresent the maximum-
block selection obligation. The optimizer's 89 certified receiver cells are a
lower witness, while this component evaluates complete actual halos and their
difference; it does not replace full geometry with their count or union area.
The exact endpoints, all source masks, residual and internal halo values are
saved for every tie. This finite candidate does not challenge a 3.5 selection
bound. Independent full geometry review is recorded by the integrated audit.

The single additional candidate `contracted-unique-actual-input.json` multiplies
all 13 outside masses of that same input by 99999999/100000000, leaving S25
fixed. Its exact input SHA-256 is
`a972d3890206cc347d67a2267473798802a1f7284cecc21af9bc6db2c68894a5`.
Complete enumeration finds a unique maximum block, the original S25. Its
full halo area is about 3.305012762967794 and true deletion area about
1.5551882171496079. It remains connected and nonterminal. The exact result
is `computation/results/round072-actual/contracted-unique.json`; no numeric
tolerance is used to declare uniqueness.

Thus a small common outside-mass contraction removes the lower-cost maximum
ties while essentially preserving the original S25 deletion area. The original
four-tie selection minimum is about 1.29637722; the contracted input's unique
selection cost is about 1.55518822. This is a meaningful finite comparison, not
a proof that the latter is extremal or close to 3.5. The general uniqueness-
by-contraction lemma and its scope belong to the core proof.

## Code and reproduction

The API is `evaluate(data,bits=72)` in
`computation/round072/actual_halo_geometry.py`. It accepts exact `positions`,
`masses`, threshold `1`, and optional `S_indices` for marking a declared block.
Zero-mass positions are explicitly removed; positive positions must be
distinct. The output preserves the mapping to original indices, evaluates all
maximum blocks, and records complete full/internal/residual rectangle data,
their area bounds, component count and terminal status. It does not rely on
producer-supplied maximum-block or eligibility claims.

From the topic directory:

    python3 -B computation/round072/actual_halo_geometry.py \
      --input computation/inputs/r072/r071-positive33.json \
      --output /Volumes/data/CodexTemp/T2/r072-positive33-recheck.json --bits 72

Repeat with `s25-only.json` for the control. The actual initial runs used
Python 3.15.0 standard library, exit 0, approximately 0.015 and 0.011 seconds
internally. These tiny jobs were run sequentially; no worker pool was created
to imply useful parallelism. Each output contains input/script hashes, exact
precision and integer coordinate scale, Python/platform and process id.
Only independent audits identified in the round's final review count as
independent reruns. All four archived finite inputs above have now been
independently replayed by direct coordinate-box enumeration with literal
source membership, followed by horizontal slab merging and direct set
difference. That code imported no production geometry module. Every complete
closure catalog, eligible label, directed rectangle and exact area endpoint
agreed, for all two/four/one maximum-block choices as applicable. No additional
finite geometry reruns are needed after those source bindings. The separate
audit package pins its code, results and receipts. This delegated component
makes no commit or push.
