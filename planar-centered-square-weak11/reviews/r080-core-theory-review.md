# R80 core theory review: actual balanced block and complete-plane label accounting

Date: 2026-11-16. This note records a separate mathematical/source review by the core author. It is not a new numerical run or a second independent area-engine execution. The audit component records its own actual replay.

## 1. Actual single-corner theorem

The final [balanced-block proof](../proofs/r080-balanced-corner-halfmass-obstruction.md), SHA256 `9ce460cd9c16979e870ba46171a5bc071fde05cd0835109cd1fd5b1a081d4022`, was read in full. Its six original source positions and masses, full closed-square qualification, unique maximum old block, eligible-hyperedge connectivity and nonterminal gap are consistent.

In particular, the source-bbox necessity argument may use a smallest containing square which captures additional sources: positivity strengthens necessity, and does not invalidate the cap. Sufficiency follows by bounding every captured subset by its bbox closure and the monotone square cap. Passing from eligible subsets to closures preserves strict eligibility. Thus the unique maximum-closure claim, together with positive masses, covers every source subset.

On the actual SW corner, every effective radius is below 1/2 and only the original old atom at the origin occurs in the local source window. This is exact localization; no old atom is moved. The seven local full labels and three complete residual labels give the two disjoint pieces displayed in the author proof. The pair residual is eligible and is explicitly subtracted. The extra thin old-only strip is positive and is retained. The radical formula and its directed rational lower bound are correct:

    D_SW=244823/500000−(171/200)a−(857/1000)(U−V)
        >129915191/10^9>1/8,
    a=√(1/8), U=√(999/4000), V=√(749/4000).

The conclusion is an actual counterexample to a general balanced-block single-corner half-local-mass bound. It is not a prescribed uniform S_n input, nor a counterexample to a four-corner aggregate or total7/2 bound.

## 2. Source review of the complete computation

The production source `computation/round080-balanced/certify_balanced.py` was read in full, together with its README and saved 63-subset catalog. The two union routines are inherited read-only from the SHA-pinned R75 core. The source correctly cancels all full outside-only labels against the identical residual before applying directed bounds. Lower deletion uses old-full inner minus residual outer rectangles; upper deletion reverses these choices. Cropping and exact noncorner subtraction preserve the full-plane decomposition.

This review did not rerun the 63-label catalog, the complete bbox caps, the union routines, or any new input. The final saved complete catalog includes the narrow cross-old label in Section 3; the computation never suffered the omission discussed there.

## 3. Complete-plane analytic decomposition and the necessary narrow label

Write d=71/200, D=171/200, w=999/1000, e=1/1000 and v=1/4. Reading the complete saved catalog shows that the old-containing eligible masks are

    1, 2, 4, 8, 15, 17, 33, 35, 49.

The three pure-outside eligible masks are16,32,48. In source order (origin, right old, top old, diagonal old, left outside, bottom outside), masks17,33,49 are precisely the three mixed origin-cluster labels. Mask15 is the old four-corner block.

Crucially, mask35 consists of origin, the right old atom and the bottom outside atom. Its mass999/1000 exceeds its squared span w²=998001/1000000. Its nonempty receiver rectangle is

    (w−U,U)×(−U,U−D),

with area

    W=(2U−w)(2U−D)>0.

The direct inequalities 4U²=999/1000>w² and 2U>D prove positivity. This rectangle is disjoint from the original origin cluster, all three other old singleton halos, the all-old central box, and the complete residual: its horizontal coordinate lies near .4995, while the origin-cluster rightmost endpoint is V<.433 and the right singleton starts at .749. Its vertical coordinates are negative, separating it from the central box and the northern singleton halos. Neither outside singleton nor outside pair reaches its horizontal interval.

The first proposed hand-derived full-plane decomposition omitted this label. The omission was caught in this review before freezing a root formula. It changes only the draft symbolic decomposition; the author's local SW proof and the producer's complete area certificate were already correct.

For completeness, the remaining analytic pieces can be checked directly. The origin-only and origin-left receivers, after their shared left-singleton residual is removed, have area

    D−2Da.

Indeed their union area is1/4+(2a−d)2a−2v(a−d+v), and its residual overlap is the same2v(a−d+v). Here v=1/4 and D=d+1/2. The origin-bottom receiver and the triple receiver form a lower component, disjoint from the upper component because D−U>a. Its area after the complete pair residual is

    (2U−d)(2U−D)+(2d−U−V)(2V−D).

The pair residual is contained in the triple rectangle. The extra part of the origin-bottom rectangle is its right strip of width V−U+d; this strip lies beyond the pair residual. The bottom singleton residual lies strictly below the entire component, because its upper vertical endpoint is −D+√(499/4000)<−U. The left singleton residual lies above it.

The other three old singleton halos contribute3/4 and the all-old receiver contributes e²; these four pieces and the narrow mask35 rectangle are mutually disjoint and disjoint from both origin components. Hence the corrected full-plane identity is

    |D|=3/4+e²+D−2Da
        +(2U−d)(2U−D)+(2d−U−V)(2V−D)
        +(2U−w)(2U−D)
       =3779121/1000000−(171/100)a−(5273/1000)U
          +(91/40)V−2UV.

The root component separately supplies its exact symbolic/interval execution and geometric note. This review supplies the label-completeness and geometric check, not that software run. In particular, an exact local formula alone must not be promoted to a complete-plane formula without accounting for mask35.

## 4. Scope

The qualified example disproves a local half-mass payment even for equal old corner masses and a unique maximum old block. Its certified complete-plane deletion remains below7/2, and its four-corner total remains below the aggregate corner allowance. No general C2 bound changes. Further rounds must keep the distinction between failure of an individual corner payment and failure of a joint block budget.
