# R56: exact logarithmic obstruction to proportional captured-mass budgets

Date: 2026-11-15 (Asia/Shanghai). Strict analytic theorem for bounded, compactly supported L1 densities. It gives the exact cost of a specified broad transport class on an explicit receiver set, not the sharp weak-(1,1) constant C2. Arbitrary nonproportional source allocation remains allowed by the original problem.

## 1. The proposed budget and its quantifiers

At threshold 1 let dmu=f(p)dp, f>=0 integrable, Q(z,r)=z+[-r,r]^2, and H={z: sup_r mu(Q(z,r))/(4r^2)>1}. A proportional witness plan on a measurable receiver set Omega subset H is any measurable probability kernel pi_z(dr) supported, for almost every z, on **actual strict witnesses** r>0 with mu(Q(z,r))>4r^2. It transports one unit of receiver demand by

    K_z(dp)= integral [1_Q(z,r)(p)/mu(Q(z,r))] pi_z(dr) dmu(p).

The load density and its budget cost are

    ell_pi(p)= integral_Omega integral 1_Q(z,r)(p)/mu(Q(z,r)) pi_z(dr) dz,
    C_Omega(pi)=ess sup_mu ell_pi.

This includes deterministic choices, any receiver-dependent tie rule, all radii, all source overlaps, and arbitrary random mixtures. By Tonelli, integral ell_pi dmu=|Omega|. A universal C_Omega(pi)<=K would therefore provide a genuine L1 area budget. This proposal uses each full capturing square's actual total mass; it is not the earlier fixed-range or pair-endpoint template.

## 2. A compactly supported L1 family and common positive-mass source patch

Fix L>16 and put

    D=[-L,L] x [0,L], A=[-1,1] x [0,1], f_L=1_D+1_A,
    S=[-1,1] x [1,2], Omega_L={(s,t):4<=t<=L/4, |s|<=t/2}.

Then 0<=f_L<=2, mu_L(total)=2L^2+2, and mu_L(S)=2. All sets have finite area; there are no singular sources or periodic tails. Boundaries carry zero source and receiver measure.

For z=(s,t) in Omega_L the radius r=t gives Q subset D and A subset Q. Thus mu_L(Q)=4t^2+2>4t^2 and Omega_L subset H_L.

**Common-patch lemma.** Every strict witness Q(z,r) at every z in Omega_L contains all S.

Indeed, the background density 1_D is at most 1. Strictness requires |Q intersect A|>0, so r>t-1. If p=(u,v) in S, then

    |u-s|<=1+t/2<=t-1<r,  |v-t|<=t-1<r.

Hence S is even contained in the square's interior. In particular an arbitrarily small grazing intersection with the excess A cannot avoid the common background mass S. This is why smoothing a single shared atom is not needed.

For a witness with r>=t, its background area is at most 2r(r+t), while its excess area is at most 2. Therefore strictness implies

    4r^2<2r(r+t)+2,   r(r-t)<1,
    r<R(t):=(t+sqrt(t^2+4))/2<t+1/t.

If r<t, the bound r<t+1/t is immediate. Since t<=L/4 and t>=4,

    |s|+r<3L/8+1/4<L,  t+r<L/2+1/4<L.

Thus every such witness remains within the horizontal and upper boundaries of D; only its lower edge can leave D.

## 3. The exact minimax cost, with a matching legal sequence

For t<=r<R(t), the square contains A in full and has exactly

    mu_L(Q(z,r))=2r(r+t)+2.

This is increasing in r, its supremum is 4R(t)^2, and it is strictly less than that supremum at every strict witness. For r<t, mu_L(Q)<=4t^2+2<4R(t)^2 as well. Hence

    sup_{r strict at z} mu_L(Q(z,r))=4R(t)^2,

with the supremum approached but not attained.

For any plan pi, the common-patch lemma gives one scalar

    I_pi=integral_Omega_L integral [1/mu_L(Q(z,r))] pi_z(dr) dz

such that ell_pi(p)=I_pi for every p in S, whereas ell_pi(p)<=I_pi for all p. Since mu_L(S)>0,

    C_Omega_L(pi)=I_pi>J*(L):=integral_4^(L/4) t/[4R(t)^2] dt.

The strict inequality holds because each probability law is supported on strictly smaller denominators, so its integrand exceeds the displayed lower envelope almost everywhere on a positive-area receiver set. It also holds if I_pi is infinite.

Conversely, for each integer k>=1 choose the deterministic radius

    r_k(t)=t+(1-1/k)(R(t)-t).

These are actual strict common-radius witnesses. Their maximum source load is exactly

    J_k(L)=integral_4^(L/4) t/[2r_k(t)(r_k(t)+t)+2] dt,

again attained on all S. The denominators increase to 4R(t)^2; dominated convergence, using the integrable r=t denominator, gives J_k downarrow J*. Consequently

    inf_pi C_Omega_L(pi)=J*(L),

and this infimum is **not attained** by a strict plan when L>16. Mixtures cannot improve the infimum.

An explicit elementary formula follows from t=R-1/R and dt=(1+R^-2)dR:

    t/(4R^2)dt=(1/4)(R^-1-R^-5)dR,
    J*(L)=[(1/4)log R+1/(16R^4)] from R=2+sqrt(5) to R=R(L/4).

In particular J*(L)=(1/4)log L+O(1). This is an exact optimization result for the stated **restricted receiver set** Omega_L; no claim of an exact full-halo transport optimum is made.

## 4. Rational scale certificate and robust impossibility

For every t>=4, all witness masses satisfy

    mu_L(Q)<4(t+1/t)^2+2
             =t^2(4+10/t^2+4/t^4)
             <=(297/64)t^2.

Set L=16*2^n, n>=1. For T=4*2^j, j=0,...,n-1, the disjoint receiver slab

    T<=t<=2T, |s|<=t/2

has exact area 3T^2/2, and every witness denominator on it is at most (297/16)T^2. It contributes at least 8/99 to ell_pi(p) for each p in S. Thus

    C_Omega_L(pi)>=8n/99.

For every proposed finite uniform K, choosing n>99K/8 refutes the budget. The argument remains valid if a plan is extended to the full halo, because all additional charges are nonnegative.

It also defeats every hybrid scheme which requires a fixed fraction theta>0 of each receiver demand to be sent through this proportional rule: even if all remaining demand is allocated arbitrarily, its load on S is at least theta J*(L). More generally a source allocation bounded below by theta times the proportional measure for each chosen witness has the same obstruction. Thus a valid uniform repair cannot retain a fixed compulsory proportional component across all these scales.

For reference, the legal r=t choice has exact maximum

    J0(L)=(1/8)log[(L^2/4+2)/66].

The universal denominator 4t^2+41/4 gives

    J1(L)=(1/8)log[(L^2/4+41/4)/(64+41/4)] <= J*(L) <= J0(L),
    0<=J0(L)-J1(L)<(1/8)log(9/8).

The rational divergence and the closed-form minimax formula are analytic proofs; finite computational audits only check their displayed arithmetic and geometry.

## 5. What this changes, and what it leaves open

This construction uses the actual total mass and arbitrary physically correct common radius, retains every overlap, and allows all measurable choices and mixtures. It therefore blocks the natural attempt to repair R50 solely by increasing the scale range and replacing individual masses by total captured mass while keeping proportional redistribution.

It does not exclude source-dependent spatial allocation. In fact Omega_L lies in the region f_L=1, so the identity allocation dz -> dmu_L(z) pays all of this receiver demand with source load exactly 1, and with no recurrence at S. That identity is supported inside every centered witnessing square, but it is not proportional to that square's full source mass. The decisive information missing from the proposed budget is which part of the local background already pays its receiver demand.

The accompanying baseline/residual interface gives a precise general formulation of this cancellation. Proving a uniform residual allocation bound that improves the inherited general C2 upper bound remains open. No C2 numerical upper or lower bound changes in R56, and the exact sharp constant remains unresolved.

## 6. Evidence and reproducibility

The independent helper's theorem and actual exact computation are archived in:

- `reviews/shared-background-proportional-independent-review.md`;
- `computation/round056/shared_background_proportional_independent.py`;
- `computation/round056/shared-background-proportional-independent-output.json`.

That run used Python 3.15.0, standard-library Fraction arithmetic, no randomness and no rational-check tolerance: 5,855 geometry/constants checks and 36 minimax/antiderivative checks, exit 0. Floating log diagnostics are numerical evidence only. The lead's separate short symbolic check and exact minimum input are recorded in the R56 manifest; they do not rerun prior R50 artifacts.

R50 dependency, read but unchanged: `proofs/nonlinear-local-potential-range-obstruction.md` SHA256 `e0038121c05263103f64b421ee13f92c6ebaf7d9a4c00a64cb4338eb6319c781`, and `proofs/pair-cumulative-triple-obstruction.md` SHA256 `79c0768da77029d38af3254cfe72466bb27c228a788fe9e1b953da4fea8490c2`. These explain the preceding failed fixed-range/pair routes; the present construction and proof do not require their conclusions as mathematical lemmas.

No commits or pushes are performed by this delegated branch. Root integrates the final whitelist and supplies the round's Git receipt.
