# R67: support rank and a mass-sensitive full erosion-union candidate

2026-11-15. Finite positive atomic inputs, distinct source positions, threshold1, strict mass greater than the area of a centered axis-parallel common-radius square. Let L be the largest mass of an eligible source subset and S any eligible block of massL. This note proves elementary support and erosion containments and states a new, unproved, stronger sufficient selection claim. It does not claim a new general weak-type bound.

## 1. Maximum blocks with at most three sources are already paid

Every actual strict witness in A captures an eligible source set of mass at mostL. Its half-side r therefore satisfies4r²<L. If x lies in the true deletion region D=H(A)\H(A\S), every strict witness at x must meet S; otherwise it survives deletion. In particular some p∈S lies within distance r<sqrt(L)/2 of x. Consequently

    D ⊂ Z_1(S):=union_(p∈S) (p+(-sqrt(L)/2,sqrt(L)/2)²),
    Delta_A(S)≤|Z_1(S)|≤#S·L.                           (1)

Thus #S≤3 implies Delta≤3L<3.5L. Coincident positions have already been merged, so #S counts distinct positive sources. This is a simple necessary-condition reduction, not a sharp general theorem.

Combined with R53 and R63, an unhandled component must be eligible-connected and nonterminal, and **every** maximum block must have both bbox widths greater than3sqrt(L)/4 and contain at least four sources. Existence of even one thin or rank-at-most-three maximum block pays that component by an already proved estimate. All these tests concern the actual component, not a disconnected whole input.

## 2. The full erosion union uses the captured mass at every label

For any eligible full-input bbox-closed set T meeting S, write

    m_T=m(T),  b_T=m(T\S),  a_T=m(T∩S).

For t≥0 define E_T(t) from the **same full-T bbox**, using radius sqrt(t)/2; set it empty when either coordinate receiver interval is empty, and set E_T(0)=empty. At t=m_T this is the usual true receiver rectangle. As proved in R66,

    E_T(b_T)⊂H(A\S),
    |E_T(m_T)\E_T(b_T)|≤a_T.                           (2)

Define the finite region

    W_A(S)=union_{eligible closed T, T∩S nonempty}
                   [E_T(m_T)\E_T(b_T)].                (3)

All eligible labels participate, including labels with small captured mass. There is no threshold theta and no filled bbox or equal-radius support container in this definition. It does not subtract the complete residual halo: each label uses only its own conservative same-bbox residual core.

The true deletion obeys

    D⊂W_A(S)⊂H(A)∩Z_1(S).                              (4)

For the first inclusion, a true receiver x has an eligible closed complete source state N(x), attained by a single strict witness. At x∈D this state meets S. Thus x∈E_N(m_N), while x cannot belong to E_N(b_N)⊂H(A\S). For the second inclusion, every constituent ring lies in its true E_T(m_T)⊂H(A). Moreover x∈E_T(m_T) implies distance(x,p)<sqrt(m_T)/2≤sqrt(L)/2 for every p∈T∩S, so x∈Z_1(S).

Since D=H\R with R=H(A\S), (4) also gives the exact overpayment identity

    |W_A(S)|=Delta_A(S)+|W_A(S)∩R|.                     (5)

Unlike the R66 filled-container budget, W introduces no points outside the full halo. Its only surplus is actual old-halo area lying in at least one label's erosion ring and paid by another residual witness. Formula(5) is a diagnosis, not the proposed estimate: computing W in(3) does not require R or Delta. Summing all a_T instead of taking the union is forbidden by the R66 fat repeated-capacity family.

A tighter optional version uses the actual bbox of the remaining subset B=T\S, without taking the complete residual halo. Write E_B(b_T) for its receiver rectangle and set it empty when B is empty. Because bbox(B) is contained in bbox(T),

    E_T(b_T)⊂E_B(b_T)⊂R.

Thus the alternative finite union

    W'_A(S)=union_T [E_T(m_T)\E_(T\S)(b_T)]

obeys D⊂W'_A(S)⊂W_A(S), and each of its ring differences still has area at most a_T. Its residual core can have a different center and need not be contained in E_T(m_T); one must compute the actual rectangle difference. The identity |W'|=Delta+|W'∩R| also holds. It remains an outer bound, since an unrelated residual source set can cover part of a ring. Neither construction is assumed stable under tiny-mass padding, and eligible closures that are not separately reachable complete N-states can contribute to their overpayment.

The containments imply the unconditional but insufficient bounds

    |W_A(S)|≤min(#S·L, (w+sqrt(L))(h+sqrt(L)))<4L

when the second expression is used. They also show that the rank-at-most-three reduction applies to W itself.

## 3. A precise next claim, deliberately stronger than true deletion

The following is a **candidate selection proposition**, not a result:

> For every finite positive atomic input A whose eligible hypergraph is connected and whose full bbox obeys d(A)²>(35/33)m(A), assume every maximum-mass eligible block S has both widths>3sqrt(L)/4 and rank at least4. Then at least one such S satisfies |W_A(S)|≤(7/2)L, where W is the explicit finite union(3).

This claim is falsifiable by complete finite geometry: to refute it one must certify its domain and a strict lower bound |W_A(S)|>3.5L for **every** maximum block, including all ties. One failed block or an outward upper bound exceeding3.5 is insufficient. Such a refutation would not refute the original Delta selection problem, because the surplus in(5) may be positive.

Replacing W by W' gives a weaker sufficient candidate with the same implication to the original Delta goal. Saving both regions in the same directed probe can identify whether a failure is only due to the old full-T bbox in the residual core. Neither candidate is asserted as a theorem.

If the claim is proved, the terminal, thin and low-rank reductions above, followed by the existing R63 induction and Borel/L1 transfer, would imply the general3.5 upper bound. This is why the claim is relevant. It does not merely rewrite Delta≤3.5: its left side is a stronger region computable from individual full-label masses and bboxes without the complete residual halo, and its possible overpayment is identified in(5).

The positive obstacle is now explicit. R66 already controls each ring by a_T, but the same source can be charged in arbitrarily many high labels even in the fat branch. A proof must control the spatial union, or a selected geometric subcover of those rings, with the source masses counted correctly. The present note supplies neither that global estimate nor a selection rule.

## 4. Bounded exact probe interface

For each label, directed rectangles yield

    E_T(m_T)^-\E_T(b_T)^+ ⊂ ring_T
                            ⊂ E_T(m_T)^+\E_T(b_T)^-.

Taking the union on each side and then measuring gives a valid interval for |W|. The actual remaining halo can be computed separately only as a diagnostic for Delta and the surplus in(5). It must not be silently used to shrink W before testing the candidate.

Two directed probes are sufficient for the initial interface: the positive ten-source R67 padding obstruction and the frozen final R64 input. No theta sweep or LP continuation is part of this test. A passing interval for either input is only an input-specific sufficient certificate. The main purpose is to verify that the construction responds to captured masses and removes the known filled-container pathology while preserving a separately measurable old-halo overpayment.

For W' the same directed difference rule applies with the genuine remaining-subset bbox. Full and residual rectangle locations are kept separately; no assumed concentric erosion is used. Any reported W'≤W comparison must allow directed enclosure slack unless a direct geometric subset calculation is made.
