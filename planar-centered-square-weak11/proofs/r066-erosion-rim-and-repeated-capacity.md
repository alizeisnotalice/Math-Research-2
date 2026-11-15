# R66: source erosion, the outer rim, and repeated capacity

2026-11-15. This attachment gives general analytic bounds and an explicit obstruction to summing their capacities. It does not prove a general cheap-deletion selection theorem. All source masses are positive, positions distinct, and halos use strict threshold1 and actual centered axis-parallel common-radius squares.

## 1. One receiver rectangle costs at most the source mass deleted from it

Fix any source set T intersecting the deleted set S. Put

    m=m(T),  a=m(T∩S),  b=m(T\S),  m=a+b,
    u=horizontal span(T),  v=vertical span(T).

Assume T is eligible, so E_T(m) is its true receiver rectangle with side lengths sqrt(m)−u and sqrt(m)−v. Using the **same full-T bbox**, form a smaller rectangle E_T(b) with radius sqrt(b)/2, discarding it if empty and setting it empty when b=0. For b>0 every point in E_T(b) has a strict witness capturing all of T\S with mass b. Hence

    E_T(b)⊂H(A\S),
    E_T(m)\H(A\S) ⊂ E_T(m)\E_T(b).                      (1)

This is a proved true residual-halo core. No proxy of uncertain containment is subtracted. Its use of the old full-T bbox is conservative because the actual residual bbox can only shrink.

Let F(t)=(sqrt(t)−u)_+(sqrt(t)−v)_+. Then

    |E_T(m)\H(A\S)|≤F(m)−F(b)≤m−b=a.                    (2)

To prove the last inequality, let d=max(u,v). F(t)=0 for t≤d²; for t>d² its derivative is

    F'(t)=1−(u+v)/(2sqrt(t))∈[0,1].                       (3)

F is continuous and locally absolutely continuous, including at d². Integrating its derivative from b to m gives(2). The u=v=0 case is simply F(t)=t, and b=0 causes no exception. This improves a loose perimeter estimate2a to the exact coefficient1.

For a low container C_theta, the stronger individual outer contribution is

    |E_T(m)\(C_theta∪H(A\S))|
       ≤ |(E_T(m)\E_T(b))\C_theta|
       ≤ min(a, |E_T(m)\C_theta|).                       (4)

Taking a **union of the erosion rings** gives a valid bound for the high unpaid union:

    U_theta≤|(union_highT [E_T(m(T))\E_T(m(T\S))])\C_theta|.  (5)

Overlaps are retained before measuring the right side. It may strictly exceed the true high unpaid area because another residual witness can cover a part of a ring. Equality of finite computed intervals does not prove a general set identity.

## 2. Near-critical bboxes can have tiny receiver area

The bound a may be much larger than the actual rectangle. Put d=max(u,v), e=min(u,v). Directly,

    |E_T(m)|=(sqrt(m)−d)(sqrt(m)−e)
             =(m−d²)(sqrt(m)−e)/(sqrt(m)+d).              (6)

If u=v=d, this becomes

    |E_T(m)|=(m−d²)²/(sqrt(m)+d)²≤(m−d²)²/m.             (7)

Thus mass close to L does not mean a large receiver region: a nearly critical two-dimensional bbox can have quadratically small area. Both(2) and(6) can be used, taking the smaller bound. Neither formula pays the union of many competitors without overlap accounting.

For the R64 final input, the four crossing high rectangles at theta3/4 are of this nearly square type. Their residual masses b are below the *full-T* bbox squared, so the conservative E_T(b) is empty. Their small unpaid areas are explained by(6)–(7), not by a large residual core inside each rectangle. The finite diagnostic does not establish this behavior for general competitors.

## 3. A universal rim is too large by itself

Let L be maximal eligible mass and S its maximum block. Every high E_T touching S lies in the radius sqrt(L)/2 expansion C_1 of bbox(S). Therefore

    U_theta≤|C_1\C_theta|
      =(1−sqrt(theta))sqrt(L)(w+h)+(1−theta)L
      ≤(3−2sqrt(theta)−theta)L.                           (8)

For theta3/4 this is at most (9/4−sqrt3)L. Combining it with the worst-case low container (7/4+sqrt3)L merely recovers4L. Compared with the target3.5L, these two *worst-case* bounds still require saving L/2 through actual uncovered rim, residual coverage or more precise shapes. It does not follow that every particular input has to exhibit a literal L/2 old-halo credit; its actual containers may already be smaller.

## 4. Crossing capacities repeat without bound even in the fat branch

The following family shows why summing the a-values in(2) cannot give a uniform source budget, even after the thin and terminal branches have been removed. It is an obstruction to that accounting rule, not a large deletion example or a failure of the actual union budget.

Let N be a positive multiple of32 and put

    S={p_i^±=(7i/(8N), ±7/16): i=0,...,N−1},
    m(p_i^±)=1/(2N),
    p_*=(39/32,0),  m(p_*)=1/4,  A=S∪{p_*}.

Then S has mass1 and bbox widths

    w=7(N−1)/(8N)>3/4,  h=7/8>3/4.

It is eligible, since its largest span squared is49/64<1. We claim that L=1 and S is the unique maximum eligible block. Any set of old sources has mass≤1, with equality only for S. Every bbox-closed set containing p_* and old sources contains either a suffix of one row or matching suffixes of both rows. The one-row type has mass≤3/4 and cannot compete with S. For the two-row type write

    T_k={p_i^±:i≥k}∪{p_*},  q=k/N,  z=1/4−q.

Its mass and coordinate spans are

    m(T_k)=5/4−q=1+z,
    d_x(T_k)=39/32−7q/8=1+7z/8,  d_y(T_k)=7/8.

When m(T_k)≥1, z≥0 and

    d_x(T_k)²−m(T_k)=(3/4)z+(49/64)z²≥0,              (9)

so T_k is not strictly eligible, including equality at z=0. The outside singleton has mass1/4. This exhausts all closed boxes: without p_* they contain a contiguous old-source interval in one or both rows; with p_* they have the forms just listed. Bbox completion preserves spans and increases mass, so an arbitrary eligible subset of mass≥1 would have an eligible closure of mass≥1. The preceding classification proves the claim over all subsets, not just the displayed tails.

For each

    k=N/4+1,...,15N/32−1

we have −7/32<z<0. The factor 3/4+(49/64)z is positive throughout this interval, so(9) is strictly negative. Also

    m(T_k)>25/32>49/64=d_y(T_k)².

Thus every selected T_k is an eligible closed **crossing** high competitor at theta3/4. Its removed and residual masses satisfy

    a_k=m(T_k∩S)=1−k/N>17/32,  b_k=1/4.

There are 7N/32−1 such tails. Hence

    sum_k a_k>(17/32)(7N/32−1) → infinity,             (10)

while m(S)=L=1. In fact the exact sum is (41/64)(7N/32−1), since the selected k are consecutive with average23N/64.

The entire input is eligible-connected: S is an eligible edge containing every old source, and any selected T_k connects p_* to S. At least one is available for every N≥32. It is nonterminal, since

    m(A)=5/4,
    d(A)²=(39/32)²=1521/1024>(35/33)(5/4).

Thus this family has a unique fat maximum block in a connected nonterminal component. The repetition obstruction genuinely persists in the unresolved branch. It excludes an estimate of the form sum_highT m(T∩S)≤cL for a fixed c, even if all internal competitors are removed. It does not exclude a geometry-sensitive or overlap-sensitive charge. The receiver rectangles overlap; equations(4)–(5) remain valid because they retain that overlap before measuring.

## 5. Actual bounded arithmetic

The author's `computation/check_r066_erosion_and_overlap.py` actually ran with Python3.15.0 and exact integer/Fraction arithmetic. It checked720 scalar boundary/erosion cases and an earlier one-row repetition family at N=8,16,32,64; runtime0.01898seconds, exit0. That auxiliary thin family is recorded in the script/output and is distinct from the fat family proved in section4.

Root's `computation/round066/check_fat_repeated_capacity.py` independently enumerated every closed bbox for N=32,64,...,384:12 inputs and1,009,644 boxes in total. It checked unique maximum mass1, both strict fat widths, all selected high tails and(10), using integer coordinates in units1/(32N) and masses in units1/(4N). Python3.13.18, wall0.144919seconds, exit0;12 workers were requested but only5 worker PIDs actually handled these short tasks, and all were verified gone. The script SHA is `524d2884d955d85f71237e5594610ee95beef4859e819f093650117120630b91`. These finite checks do not replace the proof for every N, and no full halo area was computed for this family.

The root also independently checked the erosion inequality and derivative range with the actual Wolfram v2 run archived in this round. Those computer-algebra checks supplement the elementary proof(3); they do not replace the source-to-residual inclusion in(1). The unresolved task is a shape- and overlap-sensitive budget in the fat component branch.
