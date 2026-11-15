# Point insertion: sharp old-mass square-shell budget

Date: 2026-11-15. Independent check of the covering core's argument. This file is new; the R34 platform proof, script and result remain unchanged. No commit or push.

Fix the strict threshold 1 for centered axis-parallel squares with a common radius, all r>0. Let μ be any finite nonnegative Borel measure of mass M, ν=εδ_p with ε>0, and

    A=H(μ+ν) \ (H(μ)∪H(ν)).

Define S_t(p)={x: ||x−p||∞<sqrt(t)/2}; its area is t. Then

    A ⊂ S_(M+ε)(p) \ S_ε(p),
    |A|≤M.                                                 (1)

This proves the stronger of the proposed bounds |A|≤M and |A|≤M+4ε. No assumption on source count or atomicity of μ is needed.

Proof. Take x∈A and a strict witness Q(x,r). If p∉Q, then (μ+ν)(Q)=μ(Q)>4r², contradicting x∉H(μ). Thus ||x−p||∞≤r. Total mass bounds give

    4||x−p||∞²≤4r²<(μ+ν)(Q)≤M+ε,

so x∈S_(M+ε)(p). Also x∉H(ν)=S_ε(p), proving the inclusion. The difference of the two concentric square areas is M. The inner boundary is allowed in the set difference and has area zero; the outer inequality is strict. □

For μ=Mδ_p and ε≥M, the old halo is contained in the new halo, and

    A=S_(M+ε)(p)\S_ε(p),   |A|=M.

Thus coefficient 1 cannot be reduced. For coincident points and arbitrary positive M,ε the exact area is min(M,ε).

A useful overlap-sensitive form retains the subtraction of the old halo:

    |A| ≤ M − |H(μ) ∩ (S_(M+ε)(p)\S_ε(p))|.              (2)

The discarded old-halo area in this shell is necessary if shell bounds are to be added across multiple insertions. For example, for k equal co-located atoms of mass m, the k-th insertion (k≥2) has actual pure mixed increment m, whereas (1) alone charges (k−1)m. The annuli in (1), all with inner square of area m and outer square of area km, overlap heavily. Subtracting the old halo via (2) gives exactly m. Simply summing the crude shell areas therefore does not give a uniform weak constant.

For varying insertion points, (2) identifies a concrete accounting term, but no argument here jointly controls it using source mass without repeated payments. Establishing such a control remains the next general-budget problem.

## Two explicit failures of simple cumulative accounting

First, N co-located atoms, each of mass 1/N, have final halo area 1. Regardless of insertion order, summing the crude old-mass bounds gives

    Σ_(k=1)^N M_(k−1)=(N−1)/2.

The overpayment can be unbounded. This is not merely an artifact of repeated labels at one location: take N distinct points p_j=(j/(10N),0), with masses 1/N. The same cumulative sum holds for every order, while every strict halo center is within 1/2 in each coordinate of a captured source, so the final halo lies in a rectangle of area at most 11/10. Thus merging equal-location atoms cannot fix this cumulative rule.

Second, retaining only the old-halo shell overlap in (2) can still be arbitrarily wasteful relative to actual activation. Take μ=δ_(2,0), ν=δ_(0,0). A square containing both points has r≥1, hence area at least 4>2, and cannot be a mixed witness; thus A=∅. The old halo at (2,0) is disjoint from S_2(0), so the right side of (2) is nevertheless 1. A useful cumulative allocation must retain both old coverage and actual common-capture geometry.

## A finite-source refinement retaining common-capture geometry

For finite atomic μ=Σ_j m_jδ_pj, and each nonempty source subset S, put m(S)=Σ_(j∈S)m_j. Let d(S,p) be the largest coordinate span of the bounding box of {p}∪{p_j:j∈S}. Define

    M_eff=max{m(S): d(S,p)²<m(S)+ε},

with M_eff=0 if there is no eligible subset. Then

    A ⊂ S_(M_eff+ε)(p)\S_ε(p),   |A|≤M_eff.             (3)

Indeed, any actual mixed witness captures a nonempty old subset S and p, so d(S,p)≤2r<sqrt(m(S)+ε). Thus S is eligible and ||x−p||∞≤r<sqrt(M_eff+ε)/2. This proves (3). Source atoms on a square boundary are captured; eligibility remains strict because the witness mass inequality is strict. The condition can be checked using rational arithmetic for rational positions and masses.

In the far-separated two-point example above there is no eligible subset, so (3) correctly yields zero, whereas the coarse shell bound yields 1. Formula (3) is still a single-insertion interface; it provides no uniform bound for cumulative charges through many insertions and no new value of C₂.

## A separate auxiliary obstruction

The claim |H(μ)|≥M is false: μ=δ_(−1/10,0)+δ_(1/10,0) has exact halo area 2−sqrt(2)/5<2=M, as proved in `joint-old-reserve-counterexample-notes.md`. Thus (2) cannot be simplified by replacing the old-halo intersection with total old mass.

## Actual bounded checks

Before the core's direct proof arrived, the explicitly requested 100-case seed-20261115 probe had already completed in 0.031 seconds, with no live background job. It used N=2,…,5 old atoms, M=1, ε/M∈{1/100,1/10,1,10}, rational inputs and double-precision complete subset-rectangle unions. No case exceeded |A|/M=1; the maximum observed was approximately 0.5533184645630573. This exploratory negative search is superseded by the proof of (1), not used as its justification. All 100 inputs and outputs are retained in `computation/results/joint-old-reserve/small-probe.json`; no broader search was launched.

The companion `computation/point_shell_budget_check.py` checks exact rational shell-area and coincident-source telescoping identities. It is a check of scalar arithmetic only; the geometric inclusion is proved above.
