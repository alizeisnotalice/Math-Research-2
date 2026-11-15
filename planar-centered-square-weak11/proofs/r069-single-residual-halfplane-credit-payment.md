# R69: a single-residual halfplane credit pays a new maximum-block class

2026-11-15. Finite positive atomic inputs with distinct source positions, strict threshold1 and centered axis-parallel common-radius squares. This note gives a sufficient class for |W′|≤3.5L that is not covered by the R68 support-only tests. It does not establish unrestricted W′ selection or the general C2 bound.

## 1. One residual atom and an unoccupied half-square

Let S be an eligible block of maximum mass L>0 in A. Suppose A\S consists of exactly one atom q with mass b>0. Put s=sqrt(L), let bbox(S)=[a,a+w]×[c,c+h], and write

    I=H(S),  Q=q+(-sqrt(b)/2,sqrt(b)/2)²,
    C1=(a−s/2,a+w+s/2)×(c−s/2,c+h+s/2).

Use the actual-tight-residual-bbox erosion union W′ from R67. Assume

    Q⊂C1,   I⊂{x:x_1≤q_1}.                            (1)

Then

    |W′_A(S)| ≤ (s+w)(s+h)−b/2.                        (2)

Consequently S is paid whenever the right side is at most7L/2. Either coordinate and either halfplane direction may be used by reflection.

Proof. No residual source lies inside bbox(S): adjoining one would preserve the strictly eligible spans while increasing mass beyond L. Thus internal S-only closed labels are full-input closed labels, and their receiver union is exactly I. Their residual cores are empty. Every crossing label has the same actual residual subset{q}, so its W′ core is exactly Q, regardless of the full-label bbox. Therefore

    W′ = I ∪ union_crossT (E_T(m(T))\Q),
    W′∩Q = I∩Q.                                       (3)

The R67 support containment gives W′⊂C1. By(1), the right half Q+={x∈Q:x_1>q_1} is inside C1 and disjoint from W′. Its area is exactly b/2. Subtracting this **one geometrically excluded region** from C1 proves(2). No sum of repeated label masses and no complete-halo computation is used. Boundaries have zero area and do not change the estimate.

The central exclusion theorem from R68 gives Q∩K=empty because Q is the actual entire residual halo. Hence the credit in(2) is entirely peripheral; the proof does not separately maximize and add central/peripheral costs. This is not a claim that arbitrary W′ old-halo surplus disappears. The previously proved six-source surplus example remains valid.

## 2. Finite source certificates for the hypotheses

The second condition in(1) concerns the **complete internal halo**, not a selected internal witness. Its rightmost coordinate is the supremum of the receiver right endpoints

    min_(p∈T) p_1 + sqrt(m(T))/2

over every eligible S-only closed source label T. Thus it is sufficient, for every such label, to check q_1−min p_1≥0 and

    m(T) ≤ 4(q_1−min p_1)².                            (4)

These comparisons are exact rational tests when positions and masses are rational. They are source-mass endpoint caps, not an assumption that N_A\S equals the residual complete state.

For Q⊂C1, each signed distance from q to a face of C1 must be at least sqrt(b)/2. If square-root intervals are used, an outward bound on Q and an inward bound on C1 certify this direction. Reversing those bounds would not certify containment. In the example below L=1 and all four face distances are rational, so squared rational comparisons suffice.

For true deletion there is also the simpler bound

    Delta=|H(A)\Q|≤|C1\Q|=(s+w)(s+h)−b                 (5)

under Q⊂C1. Only D⊂C1 is used for this bound. The stronger statement of interest in(2) pays W′ despite the possible internal-halo overlap with Q.

## 3. A strict one-parameter family outside the previous support tests

Take four S sources at the corners of{0,9/10}², each of mass1/4. Add

    q=(23/20,9/10),
    11/50 < b < 3457/14000.                            (6)

Then S is the unique maximum block of mass1, fat in both axes, rank4, and the full input is one connected nonterminal component.

Indeed b<3457/14000<1/4. Every proper outside-containing subset has mass at most3/4+b<1, while every proper S subset has mass below1. S is strictly eligible since its squared span is81/100<1. The full input has squared span529/400 and

    b < (33/35)(529/400)−1 =3457/14000,

so it is nonterminal and in particular ineligible. The eligible full S and the pair consisting of q and the top-right S corner connect the full input.

Every S-only pair or triple has squared span81/100 but mass at most3/4, so is ineligible. The complete internal halo therefore uses only the four singleton labels and the full S label. Their right endpoints are at most9/10+1/4=23/20; the full S receiver has right endpoint1/2. Hence the halfplane hypothesis is proved for the whole I.

C1=(-1/2,7/5)². The smallest distance of q to a C1 face is1/4. Because b<1/4, the residual half-side sqrt(b)/2 is strictly below1/4, so Q⊂C1. The theorem gives

    |W′| ≤ 361/100−b/2 < 361/100−11/100 =7/2.          (7)

The four-corner support is not a reflected monotone chain and has no bbox-center support point. Every spanning tree has three edges, each of cost at least9/10, and a side path attains minimum27/10>5/2. Its exact unit-square support union is C1, of area361/100>7/2. Thus the R68 support-tree, chain, centered-star and filled-support sufficient tests do not already pay this family. Equation(7) adds a mass-sensitive peripheral credit class. It is not a lower bound for W or a universal rule for several residual sources.

## 4. One exact source-domain check, no new full-halo experiment

Choose b=6/25 in(6). The minimal input is `computation/inputs/r069-single-residual-halfplane-payment.json`. The independent script `computation/review_r069_single_residual_payment.py` actually enumerated all31 nonempty subsets, checked the unique maximum S15, connectivity, rank4/fat conditions, and the complete five-label internal halo right endpoint. It checked Q⊂C1 by b<4(1/4)², the peripheral location by9/10>sqrt(3)/2, the support/MST obstructions and the strict bounds

    d(A)²−(35/33)m(A)=97/13200>0,
    |W′|≤349/100<7/2,   Delta≤337/100.

Python3.13.18, exact Fraction arithmetic, runtime0.000641seconds, exit0. No full H(A), W or W′ area was measured; these are analytic upper bounds supported by an exact source-domain and endpoint certificate. No extra parameter sweep, production halo job or MATLAB run was used.

The one-parameter statement(6)–(7) is proved by the displayed strict rational inequalities, not inferred from the single computed input. The general remaining task is a payment rule for arbitrary residual source configurations and all maximum-block ties; that task remains open.
