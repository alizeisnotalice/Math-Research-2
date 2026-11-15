# R66: three sources refute the fixed three-quarter envelope selection

2026-11-15. This example has a unique maximum eligible block, both its coordinate spans exceed three quarters of the square-root mass, and the whole input is one connected nonterminal component. Nevertheless its exact low-container budget B_(3/4) exceeds3.5L. The true deletion is at most1.000001L. Thus the fixed-container selection proposal in the main R66 note is false even in its intended fat branch. The general cheap-deletion conjecture and the exact C2 problem remain open.

## 1. The complete input and its domain

Put alpha=1023/1024 and epsilon=1/1024. The three sources are

    p_0=(0,0),                         m_0=1/5,
    p_1=(alpha,alpha),                 m_1=4/5,
    p_2=(alpha+1/8,alpha),             m_2=epsilon.

Let S={p_0,p_1}. Its mass is1 and both bbox spans are alpha, with

    3/4<alpha<1.

The seven nonempty subsets can be checked directly. Both S and {p_1,p_2} are eligible: their squared largest spans are alpha²<1 and1/64<4/5+epsilon. The pair {p_0,p_2} and the full input are ineligible, because their common largest span squared (alpha+1/8)² exceeds their masses. Singletons are eligible. Consequently S is the unique maximum eligible block, L=1. At theta3/4 the complete list of high bbox-closed sets meeting S is

    {p_1},  S,  {p_1,p_2}.                              (1)

The eligible edges S and {p_1,p_2} connect all three sources. The whole input is nonterminal, with the exact margin

    d(A)²−(35/33)m(A)=6982433/34603008>0.                 (2)

There are no zero-mass points, additional components, alternative maximum blocks, or thin axes that could invalidate the domain of the proposed selection assertion.

## 2. A strict lower bound for the actual envelope budget

Write r=sqrt3/4. The low container is the square

    C=(−r,alpha+r)²,
    |C|=(alpha+sqrt3/2)².

After removing S the residual halo R is the singleton square about p_2 with half-side1/64. It lies wholly in C: its greatest horizontal displacement to the right of alpha is9/64<r, and the remaining endpoint inequalities are immediate. Therefore its exact old-halo credit is

    K=|C∩R|=|R|=epsilon.                               (3)

The high singleton {p_1} has receiver square Q centered at p_1 with half-side s=1/sqrt5. Since s>r and alpha>s−r, its intersection with C is a square of side s+r. Because R⊂C, all of Q\C is in the true high unpaid union. Hence

    U≥|Q\C|=4/5−(s+r)²
             =33/80−sqrt(3/20).                         (4)

Using (3)–(4) in the exact definition of B gives

    B_(3/4)≥alpha²+93/80−epsilon
                    +alpha sqrt3−sqrt(3/20).            (5)

This is strictly above7/2 by an entirely rational certificate. Indeed

    sqrt3>433/250,       since 187489<187500,
    sqrt(3/20)<97/250,  since 187500<188180.

Both substitutions are in the correct direction in (5), and exact arithmetic gives

    alpha²+93/80−epsilon+alpha(433/250)−97/250−7/2
       =246397/131072000>0.                             (6)

In particular, B_(3/4)>7/2 for the only maximum block. Thus even the exact shape-sensitive fixed-theta selection proposal is false. The stronger shape-free criterion fails as well; directly, its compensated-high excess is bounded below by

    (33/80−epsilon−97/250)−(7/4−433/250)
       =707/128000>0.                                  (7)

No interval overlap or floating root is used in these strict conclusions.

## 3. The actual deletion is small

Let Q_b be the square centered at p_1 with half-side sqrt(4/5+epsilon)/2. The receiver rectangle for {p_1,p_2} is a subset of Q_b, as is the singleton receiver square Q. The residual R is also a subset of Q_b: its furthest horizontal displacement from p_1 is9/64, whose square is below1/5. The receiver rectangle E_S has both side lengths1−alpha and area(1−alpha)².

The complete eligible list from section1 therefore yields

    H(A)⊂E_{p_0}∪Q_b∪E_S,
    |H(A)|≤1/5+(4/5+epsilon)+(1−alpha)².

Since R⊂H(A), exact removal gives the rigorous bound

    Delta=|H(A)|−epsilon≤1+1/1048576<3.5.               (8)

The identity B=Delta+|C\H(A)| from the main note consequently shows that the envelope failure is paid by a large artificial hole term. In particular this example disproves neither the original Delta-based maximum-block selection conjecture nor a general weak-type upper bound. It does show that restricting to fat blocks does not justify filling the whole low bbox expansion at theta3/4.

## 4. Reproducible checks and the remaining route

The author's `computation/check_r066_three_source_fat_envelope.py` actually performed the seven-subset eligibility/closure census, exact connectedness and domain tests, and the rational signs (2), (6) and (7): Python3.15.0,0.000621seconds, exit0. The minimal input `computation/inputs/r066-three-source-fat-envelope.json` has SHA `7a576758b824d97b0c1274620d69b78716b1161b3a00e2774dda024b5d498cd2`. The arithmetic output is `computation/results/round066/three-source-fat-arithmetic.json`.

The separate frozen production budget and erosion programs also actually ran this one deterministic input. Their directed full-geometry certificate gives B_lower approximately3.50309040474105>3.5 and Delta_upper approximately0.999688668942791<1; both child processes exited0 and were verified gone. Those decimal summaries describe saved rational intervals and are not the proof of the signs. The producer's adapted input is archived with its own hash, independently of the author minimal-input hash above. These finite calculations check this displayed example; the elementary formulas(1)–(8) prove the obstruction.

After this counterexample, the next valid target must retain genuine low-region occupancy or deduct hole area, use a smaller adaptive low container, or bound Delta directly. Merely sharpening the high-competition payment while insisting on this fixed theta3/4 filled container cannot prove the proposed universal selection statement. No search over many new inputs is needed to establish that failure.
