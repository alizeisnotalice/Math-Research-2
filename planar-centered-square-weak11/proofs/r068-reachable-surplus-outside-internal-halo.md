# R68: reachable full source states can overpay outside the internal halo

2026-11-15. This six-source example refutes two tempting simplifications of the R67 W′ surplus: it need not lie in H(S), and restricting to the genuinely reachable complete full source state N(x) does not remove it. The example has a unique maximum block of rank4, fat in both axes, in one connected nonterminal component. It does not refute the3.5 W′ selection claim or the original true-deletion claim; no such area comparison is made here.

## 1. Rational input and its exact domain

Let S have four sources

| Index | Position | Mass |
| --- | --- | --- |
|0|(3/10,3/10)|1/10|
|1|(6/5,6/5)|449/500|
|2|(3/10,6/5)|1/1000|
|3|(6/5,3/10)|1/1000|

Add the residual sources

    q=(-1/50,3/10),  m(q)=1/10,
    z=(7/50,7/50),   m(z)=1/1000.

All six masses are positive. S has mass1 and both spans9/10, so is eligible and fat. A source set omitting index1 has mass at most203/1000<1. Any set containing index1 and an outside source has squared span at least(53/50)²=2809/2500, greater than the total mass1101/1000. Therefore all such sets are ineligible. An eligible mass-at-least-one set must consequently be S itself. Thus L=1 and S is the unique maximum block, mask15.

S is an eligible edge. The pairs{p_0,q} and{p_0,z} are eligible and connect both outside sources to it. The full input is therefore eligible-connected. Its largest span is61/50 and

    d(A)²−(35/33)m(A)=17637/55000>0.                   (1)

The support of S is the four corners of a side9/10 square. It is not a monotone coordinate chain under either coordinate reflection. The support-tree edge costs at L=1 are9/10 on the four sides and99/100 on the diagonals. Every spanning tree has three edges, so its minimum cost is27/10, attained by three side edges. The resulting tree bound is37/10>7/2, and the exact uniform support-neighborhood union has area(19/10)²=361/100>7/2. Thus neither the new support-tree test nor the filled support bound pays this particular support. This observation is not a lower bound for W′ or Delta.

## 2. An explicit positive-area patch with complete source states

Consider the open square

    P={x: ||x−(7/50,7/50)||_infinity<1/10000},
    |P|=1/25000000.                                   (2)

For every x∈P, the three-source set

    N={p_0,q,z},  m(N)=201/1000                         (3)

is captured by the common half-side r=161/1000. The maximum source distance is below1601/10000<r, and

    4r²=25921/250000<201/1000.

This is one actual strict witness. Every other S source is at distance greater than10599/10000 from x. No full-input strict witness can capture such a source: its area would be greater than

    4(10599/10000)²=112338801/25000000>1101/1000=m(A).

Equivalently, the maximum-block radius bound r<1/2 already rules it out. Thus the union of **all** full strict witness source sets is exactly N throughout P, not merely a chosen eligible label. N is also bbox-closed in A; its bbox contains precisely the three displayed sources.

The complete residual source state at x is only{z}. Indeed z itself has distance less than1/10000 and mass1/1000, so supplies a residual strict witness. Any residual witness containing q would have half-side greater than1599/10000, and

    4(1599/10000)²=2556801/25000000
                   >101/1000=m({q,z}).                 (4)

No residual strict witness can capture q. Thus N_A(x)=N but N_(A\S)(x)={z} on this entire patch. Removing S from the full complete state does not produce the complete residual state.

## 3. Surplus survives both proposed simplifications

The actual tight-bbox residual core associated with N uses the set N\S={q,z}, of mass101/1000. Inequality(4) shows x is outside its receiver rectangle E_(N\S)(101/1000). Since x belongs to the actual full-state rectangle E_N(201/1000),

    P⊂{x:N_A(x)=N}
          ∩[E_N(201/1000)\E_(N\S)(101/1000)]∩R.       (5)

In particular P⊂W′∩R, and this persists if W′ is restricted pointwise to the genuinely reachable complete full label N_A(x). The source of overpayment is an actual surviving smaller residual witness, not a nonexistent full source state.

Furthermore P is disjoint from H(S). A witness for S that captures any of p_1,p_2,p_3 would have area greater than112338801/25000000>1=m(S). The only remaining possible nonempty captured set is{p_0}, whose mass1/10 is strictly below the lower area in(4). Thus no internal S witness is strict on P, proving

    P⊂W′∩R\H(S).                                     (6)

For every maximum block, H(S)⊂W′: an outside source cannot lie inside bbox(S), since adjoining it would retain S's eligible bbox while increasing its mass beyondL. Hence every S-only bbox closure is also a full-input closure with zero residual core. Applying this observation to the example, the positive patch in(6) even yields

    |W′∩R|≥|H(S)∩R|+1/25000000.

The analogous assertion that all old surplus can be charged solely to internal-halo overlap is therefore false. The R68 central residual-free rectangle is not contradicted: this patch lies in its peripheral complement.

## 4. Scope and reproducibility

The author check `computation/check_r068_reachable_surplus.py` actually enumerated all63 nonempty subsets,24 closed bboxes and13 eligible closed bboxes, and verified the maximum/domain conditions, the exact support-tree cost and all rational patch inequalities above. Python3.15.0, runtime0.001398seconds, exit0. The final minimal input `computation/inputs/r068-reachable-surplus-six-source.json` has SHA `98c4fdc86206906a2ccf7302c56746e340e0d7141d3eab288ce449170730a307`; its actual execution output is `computation/results/round068/reachable-surplus-arithmetic.json`. This is a source-domain and positive-patch certificate, not a complete-halo-area experiment. The final four-corner support replaces an unarchived preliminary variant with one tiny source at the center; that variant was already paid by the newly proved support-star bound.

The independent reviewer separately checks the same63-subset domain and patch by exact Fraction arithmetic. No additional full-halo producer run is requested for this construction, and the earlier thin four-source toy is discarded rather than archived as another result. The full W′ area might still satisfy the3.5 bound; equations(5)–(6) only rule out two proposed surplus simplifications. A successful general payment argument must account for the distinction between full and residual reachable states, or control the actual peripheral surplus without identifying those states.
