# R68: the repeated-capacity family is paid by six fixed carriers

2026-11-15. This is an analytic result for the explicit R66 two-row family at every positive multiple N of32. It is not a bound for arbitrary inputs. It explains why the divergent sum of per-label capacities in that family does not threaten the true spatial-union budget, and why increasing N in that family cannot refute the W′ selection candidate.

## 1. The family and the label classification

Use the R66 input

    p_i^±=(7i/(8N), ±7/16),  m(p_i^±)=1/(2N),
    i=0,...,N−1,
    q=(39/32,0),  m(q)=1/4,
    S={p_i^±},  A=S∪{q}.

Put w=7(N−1)/(8N)<7/8. The R66 proof establishes that L=1 and S is the unique maximum eligible block; it is fat, and the input is connected and nonterminal. For completeness, a bbox-closed crossing set has q and a suffix of one or both rows. One-row crossing mass is at most3/4. In the two-row case, with z=1/4−k/N, its mass is1+z and horizontal span1+7z/8; when z≥0 their squared-span excess is(3/4)z+(49/64)z²≥0. Thus no crossing mass-at-least-one set is strictly eligible. All proper internal subsets have mass below1, while S is eligible.

Every full bbox-closed label meeting S belongs to one of four types: internal one-row interval, internal two-row interval, crossing one-row suffix, or crossing two-row suffix. The two signs of a single row yield six carrier rectangles in total. The outside singleton does not meet S and is not a W label.

## 2. Six carriers cover the entire W union

Write D=39/32. Receiver rectangles have the open endpoint convention. Taking closed endpoints in the following outer carriers only adds null boundaries.

* An internal one-row label has mass≤1/2, hence receiver half-side≤1/(2sqrt2). For each sign it lies in

      [-1/(2sqrt2), w+1/(2sqrt2)]
         ×[±7/16−1/(2sqrt2), ±7/16+1/(2sqrt2)].

  The sum of these two carrier areas is at most1+7sqrt2/8.

* An internal two-row label has mass≤1 and vertical span7/8. Its receiver lies in

      [-1/2,w+1/2]×[-1/16,1/16],

  of area at most15/64.

* A crossing one-row label has mass≤3/4. Its horizontal receiver interval is contained in[D−sqrt3/4,w+sqrt3/4]. For the top row its vertical interval is contained in[7/16−sqrt3/4,sqrt3/4]; for the bottom row it is contained in[−sqrt3/4,−7/16+sqrt3/4]. The sum of these two carrier areas is at most

      2(sqrt3/2−11/32)(sqrt3/2−7/16)
        =461/256−25sqrt3/32.

* A crossing two-row label has mass≤L=1 and vertical span7/8. Its receiver lies in

      [D−1/2,w+1/2]×[-1/16,1/16],

  of area at most21/256.

Every ring in W lies inside its full receiver rectangle, so these six **N-independent** carriers cover the complete W union. We may sum their six areas as a union upper bound; no equality of their union with that sum is asserted. This operation is different from summing an unbounded list of repeated source-intersection masses.

Consequently, for every N under consideration,

    |W_A(S)|≤399/128+7sqrt2/8−25sqrt3/32
              <389/128<7/2.                            (1)

The strict rational comparison uses sqrt2<10/7 and sqrt3>17/10. Both follow by squaring positive rationals, and substitution gives399/128+160/128−170/128=389/128. This proves payment for the entire family without enumerating its receiver cells at every N.

## 3. The tighter state has no surplus anywhere in this family

The remaining halo R is the singleton square about q with half-side1/4. Its left edge is31/32. We show that H(S) lies wholly to the left of53/56.

Take any nonempty internal source set T and let a be its leftmost horizontal coordinate, necessarily a=7i/(8N). Its mass is at most the total suffix mass

    m(T)≤1−i/N=1−8a/7.

For any receiver x∈E_T(m(T)),

    x_1<a+sqrt(m(T))/2
        ≤a+(1/2)sqrt(1−8a/7)≤53/56.                    (2)

To verify the last bound put t=sqrt(1−8a/7). Then

    a+t/2=53/56−(7/8)(t−2/7)²≤53/56.

The gap31/32−53/56=5/224 is positive. Thus H(S)∩R is empty.

There is only one residual source. Every crossing label T therefore has the same actual tight-bbox residual core E_(T\S)(1/4)=R, while every internal label has an empty residual core. Accordingly

    W′=H(S) ∪ union_crossT (E_T(m(T))\R).

The full halo is H(S)∪R∪union_crossT E_T. Since H(S) is disjoint from R, these formulas give the exact set identity

    W′_A(S)=H(A)\R=D_A(S).                             (3)

Thus the entire family has no tight-state overpayment. The same-full-bbox state W can still overpay, because its residual cores are smaller than R. Neither conclusion is extrapolated to general inputs with several residual sources.

## 4. The one actual N128 check and the decision it supports

The authorized N128 instance has257 sources,25,153 complete bboxes and15,205 eligible bboxes. Its specialized exact geometry computation actually completed with wall about1.35seconds and sampled peak memory about52MiB. It gives W upper approximately1.01165761049465, W′ upper and true-deletion upper approximately0.984661714465140, and a positive W surplus approximately0.0269958960295. The W′ and deletion directed arrays agree, and the H(S)-to-R separation is separately checked. All these decimals summarize saved directed rational intervals; the exact family statements are proved by(1)–(3).

The finite run verifies implementation on one larger repeated-label catalog. It is not the proof for all N. The analytic carrier bound and exact W′ identity settle this family's relevance to the3.5 candidate: no enlargement of N within the same family can produce the desired counterexample. This round therefore ends its computational extension rather than scanning more sizes.
