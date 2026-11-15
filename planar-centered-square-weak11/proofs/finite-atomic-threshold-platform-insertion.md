# Finite atomic critical platform: fixed old mass and a single vanishing insertion

Date: 2026-11-15, Asia/Shanghai. Author workspace `pcw-sol-cover-counterexample`; draft for the covering core. No commit or push. This is a counterexample to insertion budgets, not a new bound for C₂.

## Precise statement

Work in R² with Lebesgue receiver area, closed centered axis-parallel squares
Q(x,r)=x+[-r,r]², a common radius in both coordinates, all r>0, and strict threshold 1. For a finite positive measure ρ write

    H(ρ)={x: some r>0 has ρ(Q(x,r))>4r²}.

For each integer n≥1, set

    ε=1/n,   L=16n,   h=1/(16n),   a_k=1/2+kh (1≤k≤L),
    μ_n=4h Σ_{k=1}^L (δ_(a_k,0)+δ_(-a_k,0)),
    ν_n=εδ_(0,0).

Thus μ_n has exactly 32n distinct positive rational atoms, all on the horizontal axis, total old mass 8, and ν_n consists of one new atom of mass ε. Then

    |H(μ_n+ν_n) \ (H(μ_n)∪H(ν_n))| ≥ 1/4.             (1)

Consequently there is no finite absolute K such that all finite positive atomic μ, of total mass M, and all ν=εδ_p satisfy

    |H(μ+ν) \ (H(μ)∪H(ν))|
       ≤ K ε[1+log(1+M/ε)].                            (2)

In fact (1) rules out any bound B(M,ε) that tends to zero as ε↓0 at M=8, uniformly in the number and geometry of the old atoms.

## Full-scale exclusion of the old halo

Fix x=(u,v) with |u|≤1/4. If r<|v|, Q contains no source atoms. Otherwise the positive atoms captured satisfy a_k≤r+u and the negative atoms captured satisfy a_k≤r−u. Put w=|u|≤1/4 and t_+=max(t,0). Since a_k=1/2+kh and every atom has mass 4h, the elementary prefix count gives, for every r>0,

    μ_n(Q(x,r))
       ≤ 4[(r+u−1/2)_+ +(r−u−1/2)_+].                (3)

The cap k≤L can only lower this upper bound. Formula (3) is used only in the branches in which it is at most 4r²:

* If r≤1/2−w, its right side is zero.
* If 1/2−w≤r≤1/2+w, at most one term is positive, and

      μ_n(Q)≤4(r+w−1/2)≤4(r−1/4)≤4r²,

  since r²−r+1/4=(r−1/2)²≥0.
* If r≥1/2+w, both terms are nonnegative, and

      μ_n(Q)≤4(2r−1)≤4r²,

  since r²−2r+1=(r−1)²≥0.

The endpoints are included in these weak inequalities. Thus every center in the entire infinite vertical strip {|u|≤1/4}, every v, and every r>0 satisfies μ_n(Q)≤4r². In particular,

    H(μ_n) ∩ {|u|≤1/4}=∅.                             (4)

No half-radius menu, truncation, generic-position assumption, or limiting halo-set argument occurs here.

## Fixed-area strict mixed witnesses

Let

    E={ (u,v): |u|<1/4, 1/2<|v|<3/4 }.

Its area is (1/2)·2(1/4)=1/4. For x∈E use the same radius r=1. Since |u|<1/4 and |v|<3/4, Q(x,1) contains the origin. The exact old source counts at this radius are

    N_+=floor((1/2+u)/h),   N_-=floor((1/2−u)/h).

Both counts lie strictly between 0 and L, so no support-end cap is needed. The elementary bounds floor(t)≥t−1 give

    μ_n(Q(x,1))=4h(N_++N_-)≥4−8h=4−ε/2.

Consequently

    (μ_n+ν_n)(Q(x,1))≥4+ε/2>4.                         (5)

These witnesses use a common radius and the threshold inequality is strict with an explicit positive margin. Indeed the integer identity ((1/2+u)/h)+((1/2−u)/h)=L gives the stronger lower bound 4−4h for the old mass, but (5) suffices.

By (4), no point of E lies in H(μ_n). The single-atom halo H(ν_n) is exactly the open square (−sqrt(ε)/2,sqrt(ε)/2)². Since ε≤1 and every x∈E has |v|>1/2≥sqrt(ε)/2, E also misses H(ν_n). Together with (5), this proves (1).

The receiver boundary |v|=1/2 is omitted, so the n=1 strict halo endpoint causes no issue. Source squares are closed throughout, and atoms at exact capture boundaries are included by the floor counts.

## Consequence for the proposed old-mass logarithm

Here M=8 and ε=1/n, so the right side of (2) equals

    (K/n)[1+log(1+8n)] → 0,

whereas the left side is at least 1/4. Thus for every finite K some integer n violates (2). The old atom count must grow with n in this construction; no claim is made for a uniform modulus at a fixed number of sources. The example does not imply an infinite weak (1,1) constant because the old mass remains available to a general weak bound.

Compared with R22, this example preserves finite old atoms and a single new atom but replaces separated-scale activation by a single fixed receiver-radius platform. It shows that a budget depending only on old total mass and new mass can require a nonvanishing payment from old mass, even when the newly inserted mass vanishes.

## Actual calculation record

The accompanying standard-library Fraction script `computation/finite_atomic_platform_insertion_check.py` independently enumerates every exact source-distance event for fixed test centers, including the vertical capture event, and checks the old inequality at all those radii. Between consecutive events the captured mass is constant and 4r² increases, so those event checks cover all radii for each tested center. It also checks exact source masses, exact r=1 mixed counts, halo exclusion of the new atom, area 1/4, and the three analytic polynomial identities above.

This finite calculation supplements the all-center, all-n proof. It is not a finite-grid substitute for that proof. Runtime, interpreter version, inputs, source hashes and actual output are saved with the result. No random numbers, tolerances, large search, MATLAB, or extra dependencies are used. Main repository and T were not modified.
