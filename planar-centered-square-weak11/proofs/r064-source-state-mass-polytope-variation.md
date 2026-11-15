# R64: true source states, an exact mass polytope, and legal local variation

2026-11-15, Asia/Shanghai. This proof fixes a finite-position optimization problem that preserves the true centered-square geometry. It gives no new general upper bound for C2. The completed computation starts from one exactly verified positive21-source seed and certifies nine successive finite mass improvements on the same positions. No chamber-gradient certificate or local/global optimality is claimed.

## 1. The true source state is itself one actual witness

Let A be finitely many distinct points with positive masses. For each receiver x in the strict threshold-1 halo H(A), define N(x) to be the union of all source sets captured by actual strict centered common-radius witnesses at x. For every source in N(x), choose one such witness containing it. There are finitely many sources. The largest of the finitely many selected radii is still one of these strict witnesses; all the selected squares are nested because they have the same center. Its captured source set is exactly N(x): it contains their union and cannot capture a source excluded by the definition of N(x).

Consequently N(x) is nonempty, eligible, and bbox-closed. In particular

    d(N(x))²<m(N(x))≤L(A),                                 (1)

where L is the maximum eligible block mass from R63. Every source in N(x) can use this **same** actual strict witness. This strengthens the existence of separate legal witnesses for source allocations; it does not change their costs or automatically bound a Hall cut.

For any nonempty eligible bbox-closed source set T, let

    E_T=(b_T−sqrt(m(T))/2,a_T+sqrt(m(T))/2)
          ×(d_T−sqrt(m(T))/2,c_T+sqrt(m(T))/2).

These are the true source-bbox receiver rectangles. The exact state region is

    V_T={x:N(x)=T}=E_T \ union_(U eligible, bbox-closed, U strictly contains T) E_U.  (2)

Indeed x in E_T has an actual strict witness containing T, hence T⊂N(x). If the inclusion is strict, N(x) is one of the U in (2), and x∈E_(N(x)). Conversely, if N(x)=T then its actual strict witness implies x∈E_T and no strict-superset rectangle can be active. Thus the V_T are disjoint and partition H(A), and

    D(S)=sum_(eligible bbox-closed T⊂S) |V_T|.              (3)

For any maximum-mass eligible S, there is no eligible strict superset, and V_S=E_S. The E_S for distinct maximum-mass ties are pairwise disjoint: otherwise the one source state at a common receiver would contain their union, whose positive mass is larger than L.

The nonempty certified inner source unions used in R61 are eligible as well, although they need not be bbox-closed. At x in several inner rectangles, form the union T of their source labels. Choose one of these labels I whose maximum source distance from x is largest. That distance is strictly below sqrt(m(I))/2≤sqrt(m(T))/2, and bounds the distance of every source of T from x. Hence d(T)²<m(T). The empty uncertainty label remains an exception and must still be paid in every nonempty upper cut as R61 requires.

These statements retain actual common-radius geometry. Replacing the V_T or their rectangle incidence by an arbitrary hypergraph discards constraints useful to the next step.

## 2. Exact finite mass constraints for L≤1

Fix distinct rational source positions p_1,...,p_n. The variable masses are nonnegative. Ignore zero-mass positions when describing the actual positive support, but retain their locations in the fixed finite constraint list; set L=0 for the zero measure. For every nonempty source-coordinate bbox closure B, let d_B be the largest span of its tight position bbox and let m(B)=Σ_(i∈B)m_i. Equal closure masks can be deduplicated.

**Theorem 1.** The following are equivalent:

1. No positive eligible source set has mass greater than 1, that is, L(A)≤1.
2. Every closed axis-parallel square Q of side s satisfies μ(Q)≤max(1,s²).
3. Every fixed complete source-bbox closure satisfies

       m(B)≤max(1,d_B²).                                  (4)

For 1⇒2, a violating square has captured mass greater than both 1 and its area, so its captured sources are an eligible block of mass >1. For 2⇒3, place the tight rectangle of B in a closed square of side d_B. The square may include more mass, which only strengthens the required inequality. For 3⇒1, complete any eligible positive subset T of mass >1 to its tight bbox closure B. The bbox is unchanged and the mass only increases; its mass is then greater than max(1,d_B²), contradicting (4). Zero-mass extra positions do not affect these arguments.

For fixed positions all the right sides in (4) are rational constants. Thus m≥0 together with (4) is a finite rational polytope (bounded, for example, by its singleton constraints). If a prescribed set S has d(S)<1 and we impose

    m(S)=1,                                               (5)

then S itself is eligible and (4)–(5) are exactly a finite linear description of L(A)=1 with S a maximum block. No threshold menu, floating source classification, or replacement of actual captured mass is used in this equivalence.

The R63 collar and exterior growth constraints follow automatically. In particular an outside point lying in the closed collar must receive zero mass: a side-1 square can contain it and S, and (5) already exhausts that square's allowed mass.

## 3. Deliberate design constraints and source-state separation

For a targeted adversarial experiment one may impose an additional rational eta>0 and, for each closure B different from S with d_B²<1, require

    m(B)≤1−eta.                                           (6)

This is a **design restriction**, not a necessary condition on general inputs. It makes S the unique maximum-mass eligible block. Every other eligible closed state has d_B²<m(B)≤1, hence d_B²<1 and mass at most 1−eta by (6). Together with (1)–(2), this gives

    x∈E_S ==> N(x)=S,
    x∈H(A)\E_S ==> m(N(x))≤1−eta.                        (7)

In particular every strict witness at a receiver outside E_S has radius less than sqrt(1−eta)/2. This is a rigorous smaller-scale state outside the full-S core for the restricted design. There is no claim of a uniform positive eta for arbitrary inputs or all maximizing ties.

The selected R63 conjecture concerns connected inputs outside the R53 terminal class. The optimization must preserve this domain. A positive floor m_i≥tau>0 fixes the actual support bbox and prevents using empty far positions to inflate d(A). For fixed actual span d_A, the linear cap

    m(A)≤(33/35)(d_A²−gamma_geom),  gamma_geom>0,           (8)

gives d_A²−(35/33)m(A)≥gamma_geom. To guarantee eligible connectivity, specify finitely many source sets C_j whose hyperedges, together with S, connect every source, and impose

    m(C_j)≥d(C_j)²+gamma_conn,  gamma_conn>0.               (9)

These are also experiment-design restrictions. Every proposed final input should still have its actual positive support, bbox, maximum blocks, and components checked exactly. If zero masses are allowed in a later variant, recompute the actual support before invoking (8); a fixed-position bbox is insufficient in that case.

## 4. The true deletion objective and a legal local variation

Let P be the rational polytope defined by (4)–(5) and any chosen additional linear design constraints. For m in P let R=A\S and define

    F(m)=Delta_A(S)=|H(A)|−|H(R)|.                        (10)

The area difference is exact because H(R)⊂H(A). It retains the old halo and pays mass m(S)=1, including the original denominator. It cannot be replaced by the area of source states merely satisfying N(x)∩S≠∅: many of those receivers still have strict residual witnesses after S is deleted.

On a relative open chamber where all relevant masses are positive, every source-bbox eligibility decision is fixed and the orders of all x/y rectangle endpoints are fixed, F is an explicitly finite real analytic root expression. Each endpoint is a position constant plus or minus sqrt(m(B))/2. A fixed rectangle arrangement expresses union area as a finite sum of horizontal slab widths times unions of vertical interval lengths; the latter are sums of fixed endpoint differences. Both full and residual areas therefore have such expressions. Persistent endpoint equalities can be grouped if they are identities on the affine face under consideration; a wall whose order can change is not silently treated as an open chamber.

For an allowed direction h, a moving endpoint satisfies

    d/dt [c ± sqrt(m(B)+t h(B))/2] at t=0
         = ±h(B)/(4sqrt(m(B))),                           (11)

so the chamber directional derivative F'(m;h) follows by ordinary product rules. If h is tangent to the equality constraints and points inward on each active inequality, then m+t h is feasible for all sufficiently small positive t. If it also remains in the chosen chamber and F'(m;h)>0, then

    F(m+t h)>F(m) for every sufficiently small t>0.         (12)

The proof is the first-order analytic expansion with a positive linear term. This gives a legal optimization target: choose h by a bounded linear program maximizing the chamber gradient paired with h over the polytope tangent cone. It is not a claim that one LP optimizes the nonlinear area globally.

Boundary and symmetric seeds can lie on endpoint-order walls. A generic floating gradient there has no certified status. One must choose and justify a one-sided chamber, use rigorous one-sided derivative bounds, or use finite steps certified directly as below. Activating zero masses also requires separate treatment; (11) is not invoked at a vanishing mass.

## 5. Finite rational steps require no unproved gradient sign

For implementation, a floating LP or finite-difference routine may propose a nearby mass vector. It becomes a legal candidate only after rationalization, exact verification of (4)–(9), and the equality m(S)=1. Rationalization can break an active inequality; its feasibility must be checked rather than assumed. An alternative is to recover a rational feasible LP basis or shrink toward a strictly feasible rational reference point and verify all constraints.

For any prescribed rational candidate, construct directed full-halo and residual-halo area bounds

    A^-≤|H(A)|≤A^+,  R^-≤|H(R)|≤R^+.

Then F lies in [max(0,A^-−R^+),A^+−R^-]. The sharper R63 actual set-difference bounds may also be used. A strict improvement from old m to new m' is certified by

    F_lower(m')>F_upper(m).                               (13)

This proves an actual admissible finite improvement even if the proposal crossed endpoint or eligibility walls. It does not require treating a solver's derivative or objective as exact. Failure of (13) means refine arithmetic or reject that proposed improvement; it does not establish local optimality. A candidate with a lower F>7/2 would refute the maximum-block conjecture **only** if uniqueness and the nonterminal connected domain have also been certified. An expensive arbitrarily selected block among cheaper ties would not suffice.

## 6. One exact connected positive seed

The pilot uses only one fixed position configuration. Let alpha=15/16, a=alpha/2=15/32, delta=1/64. The prescribed S is the 3×3 grid {-a,0,a}². Its four corner masses are 251/1024 and its five other masses are 1/256, so m(S)=1. Add eight near sources in the axis and diagonal directions at coordinate radius

    rho_near=1−a+delta=35/64,

each of mass 1/256. Add only the four far diagonal sources at radius

    rho_far=1+a−delta=93/64,

each of mass 3/4. There are 21 positive sources. Far axis points were not included: tiny masses there would be isolated, and claiming connectivity from unused positions would be invalid.

The builder `computation/build_r064_exact_mass_polytope_seed.py` actually ran with Python 3.15.0 using integer/Fraction arithmetic. It enumerated all 274 distinct complete bbox closures, verified every cap (4), found 54 strict eligible closures, and found only the prescribed S (mask 511) at mass L=1. It verified eligible connectivity, all source positions in the R63 local region, and no outside source in the closed collar. The actual positive bbox and mass are

    d_A²=8649/1024,    m(A)=129/32,
    d_A²−(35/33)m(A)=46979/11264>0.

The exact JSON specifies eta=1/2048, a positive floor 1/8192, gamma_geom=1/16, and twelve source-edge lower constraints with gamma_conn=1/8192. Each outside source has a prescribed eligible edge meeting S, so all feasible mass vectors obeying these constraints remain connected and nonterminal. The initial seed satisfies them all. These are a concrete feasible adversarial domain rather than a claim that every input can be moved into that domain.

The final input file is `computation/inputs/r064-fixed-position-positive21.json`, SHA-256 `46acda723e8590054912a18f71c67efc95b5a0bcfd449a9c60ce615413ac41dd`. The final builder's stdout and JSON record only the verified mass geometry. No full halo objective, optimizer success, derivative sign, new C2 bound, or extremizing sequence was computed by that seed-building run.

Separately, root actually evaluated the seed's true deletion with 40-bit directed root bounds. It scaled positions by 16 and masses by 256 to call the frozen R63 complete rectangle/deletion functions; this preserves threshold 1 and all area/mass ratios. It did not invoke the N≤20 CLI or allocate all-cut arrays, and it did not modify frozen R63. The result is

    78305295413464867120688079 / 77371252455336267181195264
        ≤Delta_A(S)/m(S)≤
    78305295413696011890096831 / 77371252455336267181195264,

approximately 1.01207222228823 at the upper endpoint. The run again found 274 closures, 54 eligible rectangles, one connected nonterminal component, and unique maximum mask511. Runtime was 0.00137 seconds, exit0. The root driver and output `seed-deletion.json` are separate provenance artifacts; driver SHA `1b14abb57d8707b117aa527b8f137019214b5d959091a77efebba7320a894435`, frozen geometry SHA `5e90e5f06926b6716f56c2194b0cedeadf15d761b85173e3829d707ccc8cdf36`. This seed is still very mild. Its role is a verified feasible starting point for mass variation, not evidence near the conjectured3.5 obstruction.

### Actual finite mass improvements

The root then ran the new mass-polytope proposal/certification module with Python3.13.18, NumPy2.5.4 and SciPy1.18.1. Finite differences and one trust-radius1/16 LP proposed a rational direction; the full exact constraints and 72-bit directed complete halo bounds decided feasibility and improvement. All eight pilot line points were exactly feasible and strictly improved on the seed. The best pilot step1 raised the certified deletion cost to approximately1.1699905166170264. The full eligible closure count changed54→66, so the accepted move crossed an activation wall; it is evidence for (13), not an instance of a certified fixed-chamber derivative.

The same frozen module was subsequently called at most eight more times on the same positions and polytope. Each accepted result had to satisfy F_lower(new)>F_upper(previous), not merely beat the original seed. The actual run used all eight calls, accepted all eight winners, and ended with status `FINITE_CALL_BUDGET_EXHAUSTED`. It did not stop for lack of improvement and did not certify convergence. The driver wall time was2.72347seconds, all invoked processes exited0, and the source/input/result hash chain was retained.

Across pilot plus continuation, nine LP calls evaluated72 exactly feasible rational line points. Only42 of those72 were certified strict improvements over their respective current baselines; nine winners were accepted. The final complete deletion interval is

    27563576327882167757252602439668073627139913
      /22300745198530623141535718272648361505980416
        ≤Delta_A(S)≤
    27563576327882167757598652066915336092487621
      /22300745198530623141535718272648361505980416,

with upper endpoint approximately1.2359935097459571. Direct Fraction comparison of the final lower bound against the original72-bit seed upper bound gives a strict increase greater than1/5. Every accepted state has21 positive sources, one eligible component, the true nonterminal bbox, m(S)=1 and unique maximum mask511. The final number of eligible bbox closures is70; this is not a count of independently verified positive-area receiver states.

The last saved overlapping eligible-closure competitor bound is268181345/268435456, approximately0.999053, above3/4. It is a conservative closure bound and is not labelled a positive-area reachable-state maximum. The separate competitor-envelope lemma therefore does not by itself certify a3.5 bound for this final input, although the directly computed deletion cost is already far below3.5.

Final input and execution ledger are archived as `computation/round064/final-certified-input.json` and `computation/round064/certified-continuation-ledger.json`. Production, receipts and scope are bound by the separate R64 computation integration manifest. The optimized quantity is Delta_A(S)/m(S) with m(S)=1; it is not |H(A)|/m(A), and these improvements do not give a new lower bound for C2. The finite position family is not exhausted or proved optimal.

## 7. Evidence boundary and dependencies

The true-state argument and mass-polytope equivalence above are direct proofs. R61 provides the source-union allocation interpretation and R63 provides the selected deletion target, maximal-block collar and exact localization. Existing finite piecewise-root area structure motivates the chamber organization, but the required finite formula is derived explicitly in Section 4. No unchanged legacy script is relabelled as an R64 run.

R63 remains frozen as a separate completed bounded probe. Its 72 structural inputs were too far below the target to justify further unstructured sampling. This R64 experiment changes the question: optimize the true full deletion area inside a precisely specified admissible mass polytope with near-full S bbox and old-halo subtraction. The local finite-position domain is not a complete discretization of all measures. Exact C2, a matching general bound, and legal extremizers remain unresolved.
