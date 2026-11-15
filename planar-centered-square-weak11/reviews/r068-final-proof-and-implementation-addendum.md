# R68 final independent proof and implementation addendum

2026-11-15. PASS for the final four author proof notes, the frozen N128 specialized source/geometry interface, and the specified independent computations. General W/W′3.5 maximum-block selection and the original Δ/C2 problem remain open. The constructive strict-refinement/reachable-surplus results are proved separately in the accompanying independent review. R67 files remain frozen; no producer source, MATLAB process, commit or push was initiated by this reviewer.

## Central residual exclusion

The maximum-mass property gives d(S∪B)²≥L+m(B) for every positive residual subset B, without requiring a residual witness to be closed in the full input. For an actual residual witness, b≤L and2r<sqrt(b). If the center is in K_b, the four exact span expressions max{w,v−u,v,w−u} in each coordinate are strictly below sqrt(L+b). This contradicts the mass cap. The strict slack covers the boundary of the closed central rectangle; no almost-everywhere relaxation is needed.

The inequality f_L(b)=sqrt(L+b)−sqrt(b)/2≥sqrt(3L)/2 has equality only at b=L/3. The two sides are nonnegative before each squaring, and the final difference is (sqrt(L)−sqrt(3b))². Thus R avoids K and D,W′,W,H agree exactly inside K. The coefficient's singleton sharpness example respects the strict eligibility convention: the two-source whole set has squared span equal to its mass and is ineligible, while L=1.

The peripheral expression is valid only after checking containment of the two same-center rectangles. In the fat branch, K⊂C1 and the difference area is (1+sqrt(3))sqrt(L)(w+h)−2L<2sqrt(3)L. The general expression uses the product of the coordinatewise minimum side lengths. The proof correctly does not add this surplus estimate to an independently maximized occupied-area estimate to claim3.5. Root's saved nine-check Wolfram output/receipt is source-bound here as an external actual execution; it was not rerun by this reviewer. The witness geometry is proved in the note and was checked directly.

## Reachable-surplus note, final four-corner version

The final input changes the earlier tiny center source to the missing corner(3/10,6/5). I reread the updated proof and independently repeated all63-subset/whole-patch rational checks for this final version. Only the final version is in the delivery. It has24 closed bboxes/13 eligible closures, unique S15, both fat axes, rank4, one connected nonterminal component, and actual N_full49/N_remaining32 on the positive patch. The index subtraction N_full\S is mask48 and is not the complete residual state. All radius/cost signs and H(S) exclusion survive with a stronger far-source bound112338801/25000000.

Each support side edge costs9/10 and each diagonal99/100 at L=1. Every tree has three edges, hence minimum27/10, attained by a side path. The tree bound37/10 and exact support-union area361/100 exceed3.5. This keeps the example outside the newly paid support classes. It still makes no full W′-area or3.5-failure assertion.

The new general observation H(S)⊂W′ is valid: an outside point inside bbox(S) could be adjoined without increasing spans, contradicting maximum mass. Thus every S-only closed label is a full-input closure with zero residual core. The displayed patch is disjoint from H(S), so it proves old-halo surplus beyond internal-halo overlap even for the actual full source label at that receiver. This does not contradict central exactness because the patch is peripheral.

## All-N six-carrier family closure

The family classification is exhaustive: internal closed boxes are contiguous old intervals in one or both rows; crossing closed boxes are external-point suffixes in one or both rows. The maximality proof excludes crossing mass≥1 by a strict-eligibility equality check at z=0. All internal subsets are below the full S mass unless equal to S.

Each of the six carriers encloses the complete full receiver rectangles, hence all corresponding erosion pieces. The one-row internal radius is at most1/(2sqrt(2)), the internal two-row central height is at most1/8, crossing one-row radius is at mostsqrt(3)/4, and crossing two-row radius at most1/2. Expanding the six area expressions gives399/128+7sqrt(2)/8−25sqrt(3)/32. The rational substitutions sqrt(2)<10/7 and sqrt(3)>17/10 have the correct signs and give a strict upper bound389/128<3.5. The carrier count is fixed; this is not the invalid sum of all repeated source capacities. The argument applies to every positive N multiple of32, not merely the saved N128 run.

For every internal receiver, the leftmost source coordinate a gives m(T)≤1−8a/7 and x_1<a+sqrt(m(T))/2≤53/56, using the exact completed square in the note. Since the residual singleton halo starts at31/32, its gap from H(S) is5/224>0. Every crossing tight residual core is this same singleton halo. Distributing the set difference over the crossing union now proves W′=H(full)\R exactly for the entire family. W may still overpay. The result closes this input family's relevance to a3.5 counterexample; it supplies no statement about several-residual-source inputs.

## N128 source and independent area replay

I read the frozen specialized producer SHA b0c9bc02ea6219bb6c20a6a8bb0b3273fc0fb942859ea2f46318a3a2952e3efc. Its mass-unit/radius scaling is correct: coordinates have unit1/(32N), masses1/(4N), and squared scaled half-radii are64N·units·mesh². Internal labels have zero residual core; crossing labels use exact mass1/4 and exact same-bbox cores or the actual external singleton core. The segment-tree sweep, strict eligibility, unique maximum checks and separately diagnosed old-halo overlap have the intended direction and scope. No producer routine was imported for independent verification.

`review_r068_two_row_horizontal.py` first builds **generic coordinate-boundary bitmaps**. Intersecting x/y ranges and checking all four tight faces enumerates every closed bbox, independently of the producer's interval/suffix generator. All25,153 boxes, every stored descriptor field,15,205 eligibility decisions and139 independently binary-searched radius floors match the saved compressed catalog. There is no2^257 subset enumeration.

Area uses horizontal y strips and ordinary merged x intervals, not the producer's vertical segment tree or its four-piece rectangle subtraction. On each strip the label's full interval is cut by its same-bbox residual interval for W, or the external singleton interval for W′. Old, full and cross unions and overlaps are then measured exactly. The final actual run took1.792763seconds, Python3.13.18, exit0. All seven area interval pairs match every saved Fraction endpoint exactly:

    W upper ≈1.01165761049465,
    W′ upper = Δ upper ≈0.98466171446514,
    W old-halo surplus lower ≈0.02699589602951,
    W′ surplus exactly0 at both directed levels.

The initial replay already passed the geometry; the final replay added checks of the catalog's kind/index metadata without changing the mathematical input or search scope. Only the final entry pin/result is the formal delivery. The three source JSON/gzip files were copied bytewise into an archive-relative fixture. The entry's pinned spec accepts a repository or topic root and requires a new dataTemp output, so a clean checkout needs no historical source location.

## General support-tree and chain payment class

The final kernel note is correct for **every rooted spanning tree** on any maximum block support. A child adds at most the portion of its square outside its already-present parent's square, giving L+Σc_s. The proof does not require a monotone order, a selected witness graph or an externally supplied flow. Edge costs depend on support differences and the shared L, and are nonnegative because all coordinate spans are less than sqrt(L).

For any proposed rational-data tree, sqrt(L)P−Q≤5L/2 is equivalent to LP²≤(5L/2+Q)²; both squared sides are nonnegative. Global MST optimality is unnecessary for a positive certificate. To establish failure for every tree, the two-row example instead supplies an analytic lower bound and a tree attaining it.

For a coordinatewise monotone chain after reflections, intersections of the current square with all earlier squares are nested. The latest predecessor gives the whole prior intersection, so the chain's increment formula is exact. This yields L+sqrt(L)(w+h)−Σdx_i dy_i<3L for arbitrary rank and fat bboxes. The centered-star corollary requires an actual support point at the bbox center and at most four points; each edge is strictly below3L/4 and the support union below13L/4. It does not cover four arbitrary corner points.

The tree-method limitation uses0<epsilon≤h<1 in normalized units. Horizontal edges cost at least epsilon, cross-row edges at least h; a tree with k≥1 cross edges has cost≥2w+h+(k−1)(h−epsilon), with equality2w+h attained. Thus the eight-point example has minimum21/8>5/2 and actual support area225/64>3.5. It refutes a universal tree/support-only selection rule, not W/W′ selection.

For the three author fixtures, `review_r068_support_tree_examples.py` independently used inclusion-exclusion for support area and Prim for the MST, replacing the author's compressed-cell/Kruskal methods. All31/31/511 source-subset domain tests passed. Actual runtime0.006850seconds. Areas253/100,1063/400,225/64 and tree costs153/100,731/400,21/8 match exactly. No W or W′ was evaluated for these support-only fixtures. The source program and saved receipt were read and hash-bound separately.

## Reproduction and final scope

The accompanying constructive review gives the two small-example commands. For N128, in a clean archive checkout:

```sh
python3 -B planar-centered-square-weak11/computation/review_r068_two_row_horizontal.py \
  --output /Volumes/data/CodexTemp/T2/r068-replay-new/N128.json
python3 -B planar-centered-square-weak11/computation/review_r068_support_tree_examples.py \
  --input planar-centered-square-weak11/computation/inputs/r068-support-tree-payment.json \
  --expected planar-centered-square-weak11/computation/results/round068-support-tree-payment.json \
  --output /Volumes/data/CodexTemp/T2/r068-replay-new/support-tree.json
```

All entries use the Python standard library; the mathematical checks require assertions. The final manifest binds the four proof pins, producer/kernel/root-CAS saved receipts and source chains, portable entries/specs, actual independent outputs and process-absence checks. Saved author/producer/CAS executions are not claimed as independently rerun. No extra N, full six-source halo computation or canceled thin toy is included. No blocking issue remains in the reviewed final material. The new results restrict or pay structural classes and exhibit genuine state differences; they do not settle general3.5 maximum-block selection or C2.
