# Covering core: independent review and nonlocal old-reserve obstruction

Date: 2026-11-15 (Asia/Shanghai). Worktree `pcw-astra-covering`, based on `f9502e0df212da80e489494087f67a27645ab356`. No commit or push; root owns integration. This is a strict obstruction to candidate general budgets, not a new numerical bound for C₂.

## 1. Fixed object and exact reviewed source

Lebesgue receiver area in R², Q(x,r)=x+[-r,r]² with one common radius, every r>0, closed source incidence and strict threshold 1. Write H(ρ)={x:ρ(Q(x,r))>4r² for some r>0} and A(μ,ν)=H(μ+ν)\(H(μ)∪H(ν)).

The independently reviewed source is `pcw-sol-cover-counterexample/planar-centered-square-weak11/proofs/finite-atomic-threshold-platform-insertion.md`, SHA256 `f5c8bbf546a16241f13ab920d268ed13fe16280d68d7680558b9554544eb5c6b`. Its author states a strict all-n finite-atomic counterexample. I checked the actual source locations, source mass, every radius branch, closed capture boundaries, positive mixed margin, and single-new-atom exclusion; the theorem is verified by the analytic audit below, with separate exact computation.

For n≥1 set L=16n, h=1/L, ε=1/n=16h and

    μ_n=4h Σ_(k=1)^L [δ_(1/2+kh,0)+δ_(-1/2-kh,0)],   ν_n=εδ_0.

The 32n old atoms are distinct, positive and rational, and μ_n(R²)=8.

## 2. Analytic audit, including all scales

For |u|≤1/4, source capture implies the positive-prefix condition 1/2+kh≤r+u and the negative-prefix condition 1/2+kh≤r−u. Because every source is outside this vertical strip, the other horizontal containment inequality is automatic. If r<|v| no source is captured. Thus at every center (u,v) in the strip and every radius,

    μ_n(Q) ≤4[(r+u−1/2)_+ +(r−u−1/2)_+].

Writing w=|u|, this bound is zero for r≤1/2−w; in the one-sided branch it is at most 4(r−1/4)≤4r²; in the two-sided branch it is 4(2r−1)≤4r². The last two differences are exactly 4(r−1/2)² and 4(r−1)². Truncating the source sequence at L only decreases these bounds. All endpoints satisfy weak inequalities, so the strict old halo misses the entire infinite strip. No finite scale truncation is used.

At r=1 and |u|<1/4, the two prefix counts are floor((1/2+u)/h) and floor((1/2−u)/h). Their sum is at least L−1, since their arguments sum to the integer L. Hence

    μ_n(Q((u,v),1)) ≥4−4h

whenever |v|<1. Adding ν_n raises this to at least 4+12h>4. In particular, the author's weaker margin ε/2 is valid. The old halo exclusion and H(ν_n)=(-√ε/2,√ε/2)² then prove the author's fixed rational receiver region E={|u|<1/4,1/2<|v|<3/4}, of area 1/4, is pure mixed.

A stronger consequence uses F={|u|<1/4,|v|<1}, of area 1:

    |A(μ_n,ν_n)| ≥1−min(1/2,√ε)√ε.

For n≥4 this is ≥1−1/n, and tends to 1 although ε tends to zero and old mass stays 8. In particular, no B(M,ε) with B(8,ε)→0 can bound all finite-atomic single-source insertion areas. This includes Kε[1+log(1+M/ε)] for every fixed K. For the author's simpler 1/4 estimate, any source-count Lipschitz coefficient K_N in |A|≤K_N ε must satisfy K_(32n)≥n/4=N/128. These statements permit the source count to grow; they do not deny continuity for a fixed support.

## 3. Exact geometric cause: empty old rectangles are born in one mass layer

This is a genuine same-measure obstruction, not merely an abstract family of rectangles. At r=1, for non-grid u the two captured prefix counts i,j satisfy i+j=L−1. Their source bounding box has horizontal span d=2−h, zero vertical span, and old mass m=4−4h. Its complete old receiver rectangle has horizontal width

    √m−d=2√(1−h)−(2−h)=−(1−√(1−h))²<0,

so it is empty. After insertion, the origin lies inside this same box, the mass becomes m+ε=4+12h, and the new width is

    √(m+ε)−d=2√(1+3h)−2+h>0.

The vertical side is √(m+ε). As i varies, order 1/h translated boxes are born, all with the same mass m and the same new radius. They cover a fixed portion of the central strip. There is no logarithmic number of distinct mass scales to pay for; every contributing box may already be in one scale layer.

For comparison, an individual fixed source box with coordinate spans d_x,d_y has receiver area

    a(m)=(√m−d_x)_+(√m−d_y)_+.

It is nondecreasing and 1-Lipschitz in m: on its active branch its derivative is 1−(d_x+d_y)/(2√m), between 0 and 1, while the other branch is zero and the transition is continuous. Therefore each box gains at most ε, but the number of simultaneously born boxes is essential. Summing this correct one-box estimate cannot produce an absolute one-layer Kε bound.

The independently developed thick-rectangle lemma by `sol_cover_structure` remains compatible: its hypothesis requires all old receiver rectangles to have both sides ≥ηR. Here the relevant old rectangles are empty. The failure pinpoints its omitted geometric branch.

## 4. Fixed L1 version with receivers away from all sources

A simpler fixed L1 construction clarifies where an old-mass reserve must live. Let

    f=8·1_{([−3/2,−1/2]∪[1/2,3/2])×[−1/4,1/4]},
    g_ε=4ε·1_[−1/4,1/4]²,   0<ε≤1/4.

Then ||f||₁=8 and ||g_ε||₁=ε. The horizontal marginal of f is density 4 on the two intervals. For any square its f integral is bounded by the same prefix-length expression in §2. Thus H(f) misses {|u|≤1/4}, for every v and every r. Also 0≤g_ε≤1, hence H(g_ε)=∅.

For each (u,v) in R={|u|<1/4,|v|<3/4}, the square Q((u,v),1) contains the full vertical thickness of f and all of g_ε. Its horizontal intersection with the two f intervals has total length 1, so

    ∫_Q f=4,   ∫_Q g_ε=ε.

Thus R⊂A(f,g_ε), giving |A(f,g_ε)|≥3/4 independently of ε. Moreover, the subset E={|u|<1/4,1/2<|v|<3/4}, of area 1/4, meets neither support. On E both f and g_ε vanish pointwise. Therefore a reserve that charges old mass only at the newly activated receivers, such as a bound on this E by K||g_ε||₁+L∫_E f, fails for any finite K,L. A successful reserve must allow geometric transport from old mass away from the receiver region. This exact statement does not claim that ∫ over the entire A is zero.

This strengthens the location information in R22's local plateau example. It does not contradict the ordinary weak (1,1) inequality: order-one old mass remains available to pay the order-one new halo.

## 5. Actual independent computation

`computation/astra_covering_platform_review.py` was actually run with Python 3.15.0, standard-library Fraction only, one process, no random seed or tolerance. Cases n=1,3,17 used 273 rational test centers and 32,126 true source-capture events r=max(|a−u|,|v|). Events were independently sorted and equal radii grouped before accumulating mass. Every event satisfied old_mass≤4r². Between successive events mass is constant and the denominator grows, so this verifies every radius for each listed center, not every center. The analytic proof supplies the latter.

It also checked the exact r=1 mixed counts and 50 rectangle-integral witnesses for the fixed L1 construction at ε=1/4,1/1024. The first actual run took about 0.188 seconds, process exit code 0, status PASS. For R34 integration, provenance paths in the same arithmetic script were changed to derive from its project directory; a fresh run of every case passed again in about 0.131 seconds with exit code 0. The latter is the archived result. Results and source hashes are saved in `computation/results/astra-covering/platform-review-exact.json`.

Reproduce from this clone:

    python3 -B planar-centered-square-weak11/computation/astra_covering_platform_review.py --out planar-centered-square-weak11/computation/results/astra-covering/platform-review-exact.json

The source author's separate event run, read as evidence rather than repeated here, used n=1,2,7,32 and 42,380 events. Its analytic source and script were inspected directly; the two programs use different event arrangements. Neither computation is claimed as formal verification of the all-n theorem.

## 6. Route decision and remaining general inequality

The proposed logarithmic old-total-mass insertion repair is strictly false, even with a single new atom, fixed old mass, collinear rational sources and a common witness radius. Any uniform near-critical allocation must reserve a nonvanishing amount of old mass for a family of almost saturated source boxes and allow that reserve to be released through remote receiver regions. An exact-threshold plateau reserve alone is insufficient on the finite-atomic inputs: their maximal functions have no positive-area exact level set, yet arbitrarily small insertions activate a macroscopic near-critical region. Indeed, away from the finite support, a finite atomic maximal value is attained at a source-distance event. Equality Mμ(x)=1 therefore implies |x−p_j|∞=√m(S)/2 for some source j and some nonempty source subset S. These finitely many square boundaries have zero area.

No matching C₂ upper bound, all-input extremal inequality, or complete reserve allocation is obtained here. A precise next target is a simultaneous allocation that charges the *union* of near-critical receiver boxes to old mass without spending a source repeatedly and without summing a separate ε cost for each born box. The thick-box branch can be separated, but the empty-to-positive branch must enter the same global accounting. The existing R5 capacity interface is not a solution to that union-area problem.

Read-only sources: R22 cross-scale proof SHA256 `8a47d929f40acf09e04dba13d4f07dfc217ef54e6517b979f367c3f4c30894c8`; R22 plateau proof SHA256 `8f313777fde68863c843923da6c5b8cbdb3b859a82f54e44bde5875b7333e6e2`; reviewed finite-atomic proof hash in §1. No old T files were modified, no core halo implementation was changed, and main was not written.


R34 integration supplement: the final thick-rectangle proof `proofs/sol-cover-structure.md`, SHA256 `4c3aa93e4987bc2975d9d5f8b8724d8ea19fcf6231b9aed96123345c14e7d1c7`, has been analytically reviewed. Its two sequential slice expansions give 64R²δ/a+16Rδ²/a, which with δ≤ε/(8R), a=ηR, ε≤R² is at most 33ε/(4η). Its final state correctly records the unrestricted logarithmic and degenerate-layer budgets as disproved. The abstract rectangle examples are auxiliary geometry checks, not claimed as shared-source realizations.
