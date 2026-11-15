# R63: maximal compact blocks, an empty source collar, and a 3.5 deletion target

2026-11-15, Asia/Shanghai. This note proves the structural reduction below for **all finite positive atomic inputs**. Its proposed 3.5 deletion inequality remains unproved. The elementary general bound 4 obtained in the proof is weaker than the inherited general bound 3.590865312933, and is not a numerical improvement for C2.

## 1. Objects and the maximal block

Write A for a finite set of distinct source positions, with positive masses m_i; merge coincident positions first. For T⊂A put m(T)=Σ_(i∈T)m_i and let d(T) be its largest coordinate span. The true halo H(A) uses all centered, axis-parallel squares with one common radius, all r>0, and strict threshold 1. Call T eligible if d(T)^2<m(T). Define

    L(A)=max{m(T): ∅≠T⊂A, d(T)^2<m(T)}.

Singletons make this finite maximum positive. A maximizing block S has mass L=L(A), and each of its two coordinate spans is strictly below sqrt(L). Let its bbox be [a,b]×[c,d], and write w=b−a, h=d−c. Positivity implies that S contains every source of A in its bbox: adding any omitted such source preserves the bbox and eligibility and increases the mass.

**Witness and square caps.** Every actual strict witness captures an eligible subset, so its captured mass is at most L and its radius is strictly less than sqrt(L)/2. Moreover every closed axis-parallel square of side sqrt(L) contains at most L mass. Otherwise its captured sources have span squared ≤L and mass >L, and form an eligible block larger than L. These are consequences of maximizing over *all* eligible blocks; maximizing only over blocks of bounded cardinality does not establish either conclusion.

More generally, every square of side s≥sqrt(L) has mass at most s²; otherwise its captured sources form an eligible block of mass >s²≥L. A useful quantitative exterior constraint is

    m(T)≤d(S∪T)^2−L       for every nonempty T⊂A\S.       (1a)

Indeed S∪T has mass >L and therefore cannot be eligible. In particular any square of side s containing the whole bbox of S has outside-source mass at most (s²−L)_+: for s≤sqrt(L) that outside mass is zero by the unit-square cap, and for s≥sqrt(L) use the preceding larger-square cap. The large-square condition is essential; a hypothetical critical background cannot be inserted outside the collar without checking it.

**Empty source collar.** No source of A\S lies in

    J_S=[b−sqrt(L),a+sqrt(L)] × [d−sqrt(L),c+sqrt(L)].       (1)

Indeed a point p in (1) can be included together with the whole bbox of S in a square of side sqrt(L): independently choose its left coordinate in [b−sqrt(L),a]∩[p_x−sqrt(L),p_x], and its lower coordinate analogously. Both intervals are nonempty by (1). The square then has mass at least L+m_p>L, contradicting the square cap. Closed endpoints cause no exception because masses are positive and the inequality is in the mass.

## 2. A true deletion cost with the old source area retained

For any nonempty S⊂A define

    Delta_A(S)=|H(A)\H(A\S)|.

This is the full halo area that disappears after removal. It is neither the pure mixed area nor the area of a point proxy. Monotonicity gives the exact identity

    |H(A)|=Delta_A(S)+|H(A\S)|.                            (2)

If S maximizes L, every receiver in the difference has an actual witness meeting S; otherwise it would remain in H(A\S). The witness cap implies the receiver belongs to

    I_S=(a−sqrt(L)/2,b+sqrt(L)/2)
          ×(c−sqrt(L)/2,d+sqrt(L)/2).

Consequently

    Delta_A(S) ≤ (w+sqrt(L))(h+sqrt(L)) < 4L.              (3)

The last inequality is strict because both source spans are strictly below sqrt(L). Repeatedly remove a maximum-mass eligible block. The process terminates, the removed source blocks partition A, and the actual deleted receiver sets are disjoint by nested halos. Thus (2)–(3) prove |H(A)|<4m(A) for every nonempty finite input. This is a proved baseline for this mechanism, **not** a replacement for the stronger inherited upper bound.

There is already a useful exact geometric branch at target K=7/2:

    min(w,h) ≤ 3sqrt(L)/4  ==>  Delta_A(S)<(7/2)L.         (4)

For example if w≤3sqrt(L)/4, the two factors in (3) are at most (7/4)sqrt(L) and strictly less than 2sqrt(L). The test is rational: 16 min(w²,h²)≤9L. A less conservative test uses the full product in (3), with an outward root bound. Only maximum blocks whose two coordinate spans both exceed 3sqrt(L)/4 can obstruct this elementary 3.5 branch.

The exact deficit formulation, which retains the already-paid old halo, is

    Delta_A(S)=|I_S|−|I_S\H(A)|−|I_S∩H(A\S)|.           (5)

Thus the hard local target is to pay the excess |I_S|−(7/2)L by uncovered area and remaining-source halo inside I_S. The two subtracted sets in (5) are disjoint. Subtracting the mass of A\S in place of either area would be invalid.

**Exact spatial localization.** Put

    R_S=[a−sqrt(L),b+sqrt(L)] × [c−sqrt(L),d+sqrt(L)],
    A_loc={i∈A:p_i∈R_S}.

Then L(A_loc)=L because S remains and no eligible mass can increase on taking a restriction. Both A and A_loc have no strict witnesses of radius at least sqrt(L)/2. Every smaller-radius square centered in I_S lies in R_S, so the two original halos agree on I_S; their residual halos after removal of S agree there as well. Since both deletion differences lie in I_S, this proves

    Delta_A(S)=Delta_(A_loc)(S).                           (5a)

Furthermore m(A_loc)≤9L: each coordinate width of R_S is strictly below 3sqrt(L), partition it into three intervals of length at most sqrt(L), assign boundary points uniquely, and cover each of the nine rectangular cells by a closed square of side sqrt(L). Each cell has mass at most L by the square cap. Thus the target for a specified maximum block is scale-free and spatially local: after normalization L=1 it involves at most mass 9 in a square of side below 3, with the exact empty collar and all large-square caps above. It still permits arbitrarily many sources and arbitrarily small masses. It does not assert that a fixed spatial grid or fixed source count is complete, or that halo area is uniformly continuous under such compression.

## 3. The one selected next lemma

R53 supplies B=1994727599733432151936139/604462909807314587353088<33/10 and the proved local bound

    |H(C)|≤B max(m(C),d(C)^2).

Hence each source component with d(C)^2≤(35/33)m(C) is already paid at a cost strictly below (7/2)m(C). Components here are those of the hypergraph of all eligible source subsets. Their true halos are disjoint: two strict witnesses at the same receiver can be nested by radius, and the larger gives an eligible edge joining their source components.

**Selected conjecture (maximum eligible block).** For every eligible-connected finite A satisfying d(A)^2>(35/33)m(A), there is a proper maximum-mass eligible block S for which

    Delta_A(S)≤(7/2)m(S).                                  (6)

Every maximum block is proper in this domain since A is not eligible. The collar (1), the global witness cap, and (5) are concrete geometric constraints available for proving (6). The conjecture asks for at least one of the maximum blocks if several tie. Failing a test on one arbitrarily chosen tie does not refute it.

**Conditional general-bound theorem.** If (6) holds for every such finite input, then C2≤7/2 for arbitrary nonnegative L1 inputs.

Proof for finite atoms is by induction on source count. Split into eligible components. Stop on the components paid by R53. For any remaining connected component choose S from (6); its residual has fewer sources, so use induction on the residual and (2). The actual removed mass is m(S); no outside reserve is duplicated. This yields |H(A)|≤(7/2)m(A).

For completeness, a finite-atomic bound transfers to every finite positive Borel measure. First restrict the measure to a bounded square. Partition that square into a finite grid of cells of l-infinity diameter at most epsilon and move each cell's entire mass to a representative point in its closure. Every original captured source point moves by at most epsilon. Given a strict witness Q(x,r), the discretized measure on Q(x,r+epsilon) is at least the original captured mass, and exceeds 4(r+epsilon)^2 for all sufficiently small epsilon. Thus every point of the original halo lies in the liminf of the discretized halos. Fatou's lemma for receiver indicators proves the bound with the same preserved mass. Increasing bounded restrictions recover the full measure: a strict witness has mass strictly above its threshold and hence some finite restriction still has a strict witness at the same radius. Threshold lambda follows by replacing the measure by mu/lambda. In particular the result applies to all L1 densities. This argument uses strict witness slack and finite source displacement, not convergence of halo areas without justification.

The weaker conjecture permitting *any* compact proper S would also suffice. Bounded-cardinality probes test stronger restrictions of that weaker conjecture; they do not enumerate all possible compact blocks. The maximum-mass version above is the selected next target because it has the proved cap and collar, rather than only an existential reformulation of a good deletion tree.

## 4. A convergent finite certificate for a prescribed deletion

Use rational inner/outer approximations H_X^-⊂H(X)⊂H_X^+ to the full finite bbox rectangle union, for X=A and X=A\S. Their construction must include every bbox closure and use its actual captured mass. Then

    ell(S)=|H_A^-\H_(A\S)^+| ≤ Delta_A(S)
          ≤ |H_A^+\H_(A\S)^-|=u(S).                    (7)

If Gamma_X is the sum, over all retained bbox rectangles for X, of their outer area minus inner area, then

    0≤u(S)−ell(S)≤Gamma_A+Gamma_(A\S).                   (8)

Indeed any point counted by the upper difference and not the lower is either in H_A^+\H_A^- or in H_(A\S)^+\H_(A\S)^−. Each such union uncertainty is bounded by the corresponding sum of rectangle bands. For a fixed finite input these quantities tend to zero as the rational root brackets shrink. Dividing (8) by m(S) gives the ratio error; this need not be uniform when the deleted mass tends to zero.

Thus u(S)≤(7/2)m(S) is a rigorous successful deletion certificate. To refute (6) for an input requires ell(S)>(7/2)m(S) for **every** maximum-mass eligible block S, together with exact connectedness and the nonterminal check. A floating minimum above 3.5 is only a request to refine the arithmetic. R61 source-union labels do not by themselves determine a residual halo: after deleting S, witnesses with smaller actual remaining mass must be recomputed.

## 5. Relation to earlier work and next computation

R58 preserves old coverage by a global mask and true leaf cores. The exact increments (2) instead retain the full residual halo; telescoping across a chain prevents repeated receiver payment automatically. This gives a stronger input-specific certificate at a higher geometry-computation cost, not an automatic tree-choice theorem.

R59 establishes a positive costed local boundary insertion, while earlier platform examples allow arbitrary specified insertions to have large cost relative to their tiny added mass. Neither proves nor refutes (6): (6) chooses a maximum-mass compact block, and pays only its *full true deletion difference*. No arbitrary small insertion is asserted cheap. The old sources remain in the subtraction in (5) and (7).

The bounded probe assigned to next_compute uses exact rational positions/masses, complete residual bbox geometries, and (7), recording all maximum-mass ties as well as a separately labelled small-rank candidate menu. Structural inputs are wide connected configurations, two-scale clusters with bridges, frames/thin shells, and light boundary additions. Already-paid components are identified before judging failures. Only a fully certified counterexample can decide the conjecture negatively; successful finite inputs do not prove it. Root owns the actual large executions and their receipts. No frontier result is asserted in this proof before that receipt exists.

A separate small exact structural check actually ran in this workspace using `computation/check_r063_maximal_block_structure.py`, Python 3.15.0, standard integer/Fraction arithmetic, 32-bit strict root brackets, no randomness and no floating decisions. The final run checked nine inputs (singletons, separated points, a chain, a bridge, three square grids, an asymmetric frame, and the bridge with distant extra sources), all 1,981 nonempty source subsets, and all 25 maximum-mass eligible blocks. Square caps, bbox closure, empty collars, exterior growth constraints (1a), true residual union-difference bounds and the Gamma width bound passed. Local restrictions retain maximal eligible mass L and satisfy mass≤9L; independently recomputed localized deletion intervals intersect the original intervals. Final runtime was 0.1116 seconds, exit 0. This small calculation supports the implementation of the stated finite identities and does not certify (6) beyond those inputs. The raw exact output is `computation/results/round063/maximal-block-structure-exact.json`; decimal displays are not used in decisions. Earlier two development runs had eight cases; their counts are not added to this final-run ledger.

Read inherited sources and SHA-256 pins:

- R53 independent Borel/local proof, `reviews/r053-covering-core-independent-review.md`: `f466f85312c98c116f8276a10b14fbde3f8d9870be1ac5f970c4e196d76ece66`.
- R58 author proof, `proofs/r058-bbox-terminal-masked-tree.md`: `00eb0f728faa9738e934715b45e52403a3521b0552edc7fc30c7976217f6264f`.
- R59 author proof, `proofs/r059-mass-costed-boundary-reallocation.md`: `d9c548c2b58307e29168d2c687386bea1467e88ce320784899bc839527400169`.
- R61 author proof, `proofs/finite-source-union-cut-certificate.md`: `05d92b73f5b8bb077f8fbbb9b56e68b29c164758aeb79d7d1e062b83caf23e1e`.

These dependencies were read, not rerun here. Exact C2, matching legal extremizers, and the conjecture (6) remain open. No commit or push is performed by this delegated workspace.
