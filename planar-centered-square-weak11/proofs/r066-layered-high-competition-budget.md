# R66: exact mass layers and the obstruction to a filled low container

2026-11-15. All objects are finite positive atomic inputs, with coincident atoms merged, for strict threshold1 and all centered axis-parallel squares with a common radius. The exact C2 goal remains open. This note gives a computable general theorem, identifies the precise relaxation loss, and states the remaining selection problem without turning finite probes into a general bound.

## 1. Source geometry and the two layers

Let S be a maximum-mass eligible block of A, with mass L>0, bbox [a,b]×[c,d], and widths w=b−a, h=d−c<sqrt(L). Write

    H=H(A),  R=H(A\S),  D=H\R,  Delta=|D|.

For theta in [0,1], let the low-mass container be

    C_theta=(a−sqrt(theta L)/2,b+sqrt(theta L)/2)
              ×(c−sqrt(theta L)/2,d+sqrt(theta L)/2).

At theta=0 this is the open bbox, possibly empty in area. For every eligible bbox-closed source set T, E_T denotes its complete receiver rectangle at its own true mass. Form the finite high-mass union

    P_theta=union{E_T: T eligible and bbox-closed,
                       T intersects S, m(T)>theta L}.       (1)

**S is included when theta<1.** It may be omitted only after separately proving E_S⊂C_theta. In general omission is false: for a singleton S, E_S has area L and C_theta has area theta L, so omission would miss a positive part of the true deletion.

R64 proves that the complete source union N(x) at any halo receiver is realized by one actual strict witness. It is eligible and bbox-closed. For x∈D it intersects S; otherwise this witness survives removal of S. If m(N(x))≤theta L, its strict radius is below sqrt(theta L)/2 and it meets S, hence x∈C_theta. Otherwise x∈E_(N(x))⊂P_theta. Conversely P_theta⊂H, since each E_T supplies an actual strict witness.

This proves the **exact exterior identity**

    D\C_theta = P_theta\(C_theta union R).                 (2)

It preserves all repeated source use and all receiver overlaps. No arbitrary eligible hypergraph replaces the geometric rectangles.

## 2. A complete finite budget and its exact loss

Define

    U_theta=|P_theta\(C_theta union R)|,
    K_theta=|C_theta intersect R|,
    B_theta=|C_theta|−K_theta+U_theta.                     (3)

The letter K here is an old-halo area credit, not a new weak-type constant. Equation(2) gives

    Delta ≤ B_theta,
    B_theta = Delta + |C_theta\H|.                        (4)

Indeed the exact low deletion area is |C_theta∩H|−|C_theta∩R|, while the remaining deletion area is U_theta. Replacing |C_theta∩H| by |C_theta| creates exactly the hole area |C_theta\H|. Equivalently B_theta=|(C_theta∪P_theta)\R|.

Thus every fixed input and maximum block has an explicit finite rectangle-union certificate. It is stronger than a sum of independent witness costs because the high rectangles are unioned before charging. It is still a sufficient certificate, and its failure is not automatically failure of the original cheap-deletion conjecture: holes in the low container can consume the entire excess.

There is an important limitation on threshold search. Since C_theta increases with theta while the actual H and Delta are fixed, (4) implies

    theta1≤theta2  ==>  B_theta1≤B_theta2.                 (5)

Consequently exhaustive optimization of this **exact** budget over theta merely favors the smallest container, with the minimum at theta=0. It does not create a new selection theorem. A useful mass layer must instead allow a uniform estimate of U_theta−K_theta, or a controlled simplification of the high source states. If the low container were a union of small source-centered squares rather than the bbox expansion, its theta0 limit would be empty and the corresponding exact budget would simply be Delta itself.

## 3. The fat component reduction and a refuted sufficient selection proposal

First split A into its actual eligible hypergraph components; their true halos are disjoint. R53 already pays every component with d(A)²≤(35/33)m(A) below3.5m(A). In a remaining connected component, any maximum block having w≤3sqrt(L)/4 or h≤3sqrt(L)/4 is paid by the R63 thin-axis bound. Therefore a genuinely unresolved component has **every** maximum block fat in both coordinates:

    w>3sqrt(L)/4,  h>3sqrt(L)/4.                          (6)

The following sufficient geometric selection statement was proposed during this round:

> In every such connected nonterminal component satisfying(6), some maximum block S has B_(3/4)(S)≤(7/2)L.

The three-source fat example proved in the separate attachment **refutes this statement**: its unique maximum block has B_(3/4)>3.5L. The implication to a general3.5 bound would follow from R63's conditional induction and Borel/L1 transfer, but its sufficient premise is false. The original selection question asking Delta≤3.5L remains open; it is strictly weaker because it does not include the artificial hole term in(4).

For a prescribed S, the exact shape-sensitive condition is

    U_(3/4)−K_(3/4)
       ≤ (7/2)L−(w+sqrt(3L)/2)(h+sqrt(3L)/2).             (7)

The cruder shape-free sufficient inequality is

    U_(3/4)−K_(3/4)≤(7/4−sqrt(3))L,                     (8)

because |C_(3/4)|<(7/4+sqrt(3))L. Both sufficient conditions(7) and(8) fail as selection assertions, even in the connected nonterminal all-fat branch, by the three-source example. The earlier separate five-source counterexample only refutes(8) without the fat restriction: its unique maximum S is thin and already paid by R63, and its actual B is below3.5L. Neither example refutes the original Delta selection conjecture or any general C2 bound. The two examples have different scopes and must not be conflated.

A diagnostic on a disconnected input must not use a fat maximum block from an independently paid compact component to claim an obstruction in this branch. Caps, maximum blocks, connectedness and terminal tests must be recomputed for the actual component under examination.

## 4. Exact finite arithmetic certificates

Let C^-⊂C_theta⊂C^+ be directed rational containers, P^-⊂P_theta⊂P^+ the unions of all selected high rectangles, and R^-⊂R⊂R^+ the complete residual halos. High-label membership uses the exact rational comparison m(T)>theta L; equality belongs to the low case. Residual source masses and bbox closures must be recomputed after removing S.

Valid intervals are

    |P^-\(C^+∪R^+)| ≤ U_theta ≤ |P^+\(C^-∪R^-)|,
    |C^-∩R^-| ≤ K_theta ≤ |C^+∩R^+|,
    |(C^-∪P^-)\R^+| ≤ B_theta ≤ |(C^+∪P^+)\R^-|.         (9)

The direct union-difference interval for B is usually tighter than adding separate interval endpoints in(3). The hole interval is

    |C^-\H^+| ≤ |C_theta\H| ≤ |C^+\H^-|,                (10)

where H^-/H^+ bound the full true halo. Independent interval compatibility with Delta plus the hole is a useful implementation check; the identity itself is proved by(2), not by overlapping decimal intervals.

If Gamma_C is |C^+|−|C^-| and Gamma_P,Gamma_R are sums of rectangle-band areas, the width of the B interval in(9) is at most Gamma_C+Gamma_P+Gamma_R. Any point in the upper difference but not its lower counterpart lies in a container/high uncertainty band or a residual uncertainty band. For fixed finite input these bands shrink to zero. This convergence is not uniform in source count or geometry.

To certify(8) with a rational upper bound g for (U−K)/L, one can use only rational comparisons: require q=7/4−g≥0 and q²≥3. Conversely a rational lower bound g violates(8) if q<0 or if q≥0 and q²<3. No unrounded decimal sqrt3 is used in a strict decision.

## 5. Interpretation of the new tests

The R64 final input has only five high closures at theta3/4, including S. Four crossing closures have mass close to L but bbox spans close to their critical square-root mass; their receiver rectangles are very small. A large L_star from R64 therefore does not imply a large unpaid high area. This observation motivates the separate erosion/eligibility-slack estimate, but finite examples do not prove a uniform payment rule.

The required general work must also retain enough of the low-container hole credit, adapt the low container, or estimate Delta directly. Bounding the *union* U_theta and retaining K_theta alone cannot establish the refuted fixed theta3/4 filled-container selection. Neither a sum of all source-intersection masses nor a threshold sweep of the exact identity(4) closes the original Delta-based obligation. The explicit obstacles and the actual finite diagnostics are archived separately from this theorem. No new general3.5 bound, exact C2 value or matching extremizing sequence is asserted.

## 6. Final status update: the fixed-theta selection proposal is false

The three-source construction in `r066-three-source-fat-envelope-obstruction.md` proves the refutation recorded in section3. It has a unique maximum block of mass1 with both widths1023/1024, and the complete input is connected and nonterminal. Its exact B_(3/4) is strictly above7/2: a rational lower certificate for the excess is246397/131072000. Thus the shape-sensitive selection proposal in section3, and consequently its stronger shape-free version(8), are false even in the all-fat branch. The earlier five-source example remains a separate thin obstruction with the limited scope stated there.

The three-source true deletion is at most1+1/1048576. Therefore the example does not refute the original Delta-based cheap-deletion selection conjecture. The exact identities(2)–(5), interval certificates(9)–(10), and component/thin reductions remain valid. What fails is the proposed uniform replacement of the true low deletion by the entire fixed theta3/4 low container. Subsequent positive work must retain its genuine occupied part or deduct enough of the explicit hole term, adapt the low container, or return to Delta itself. The final status is a proved finite-budget identity and a proved obstruction to this particular sufficient selection route, not a new general C2 upper bound.
