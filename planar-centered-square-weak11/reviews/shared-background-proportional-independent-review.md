# R56 independent audit: shared-background proportional transport

## Scope and status

**Strict analytic obstruction, genuine compactly supported bounded L1 density.**
Independently checked in this helper workspace on 2026-11-15. This concerns the rule that sends each receiver's unit demand proportionally to *all source mass in an actual strict capturing square*, allowing arbitrary adaptive witness choices and probability mixtures. It does not rule out arbitrary nonproportional source allocation, and is not an upper or lower bound on the weak-(1,1) constant C2.

The construction was supplied by the recursive lead; this document independently checks its quantifiers, geometry, and exact dyadic constants. The sharper exact restricted-receiver minimax formula below is an additional derivation.

## Object and all-witness geometry

Fix L>16, threshold lambda=1, and define

- D=[-L,L] x [0,L], A=[-1,1] x [0,1];
- f=1_D+1_A, dmu=f dx, so mu(R2)=2L^2+2 and 0<=f<=2;
- S=[-1,1] x [1,2], with mu(S)=2;
- Omega={(s,t):4<=t<=L/4, |s|<=t/2}.

For every z=(s,t) in Omega, r=t gives Q(z,r) contained in D, contains A, and has mass 4t^2+2>4t^2. Thus Omega lies in the strict halo.

Let any r>0 be a strict witness at z: mu(Q)>4r^2. Since 1_D<=1, Q must meet A in positive area, implying r>t-1. This forces every such Q to contain all of S, up to null boundaries:

    s-r < 1-t/2 <= -1,
    s+r > t/2-1 >= 1,
    t-r < 1,   t+r > 2t-1 >= 7.

This is the key point: grazing the extra density does not avoid the common positive-mass background patch.

If r>=t, the background mass is at most 2r(r+t) and the excess mass is at most 2. Strictness therefore implies

    4r^2 < 2r(r+t)+2,
    r(r-t)<1,
    r<t+1/t.

If r<t the same last upper bound is immediate. For *every* witness,

    mu(Q) <= 4r^2+2 < 4(t+1/t)^2+2
           =4t^2+10+4/t^2 <=4t^2+41/4.

For t>=4, the normalized bound is also

    mu(Q)/t^2 <=4+10/t^2+4/t^4 <=297/64 <5.

All inequalities use the total actual source measure, not a single-atom mass, a column mass, or a fixed geometric reach.

## Adaptive witnesses and probability mixtures

For almost every z choose an arbitrary probability law pi_z supported on its actual strict witness radii. The unit receiver charge is

    K_z(dp)= integral [1_Q(z,r)(p)/mu(Q(z,r))] pi_z(dr) dmu(p).

The induced source load measure eta_Omega is integral_Omega K_z dz. For every p in S, all witness indicators are one, so Tonelli gives the uniform lower bound

    d eta_Omega / dmu (p)
      >= integral_Omega dz/(4t^2+41/4)
       = (1/8) log[(L^2/4+41/4)/(64+41/4)] = J1(L).

The same lower bound applies to any extension that transports the full halo. Thus no absolute C can make eta<=C mu for every bounded compactly supported L1 input by this proportional rule, regardless of adaptive radii or probability mixtures. It already fails on one fixed source patch of area two.

## Exact rational dyadic certificate

For L=16*2^n, partition Omega into disjoint trapezoid slabs T<=t<=2T, |s|<=t/2, T=4*2^k, k=0,...,n-1. Each slab has area 3T^2/2, while every witness denominator there is at most (297/16)T^2. It therefore adds at least

    (3/2)/(297/16)=8/99

to every p in S. Hence the exact rational bound is load>=8n/99. The simpler constant 3n/40 follows using denominator 5t^2. These are analytic bounds for every n, not conclusions extrapolated from the finite script audit.

## Matching upper choice and exact restricted-Omega minimax

The canonical r=t rule has maximum source load exactly

    J0(L)=(1/8) log[(L^2/4+2)/66],

attained on all of S. Every other source indicator is at most one. Moreover

    0<=J0(L)-J1(L)<(1/8)log(9/8).

Thus the least possible maximum source load on Omega grows as (1/4)log L with a bounded, explicitly controlled error.

In fact one can identify the restricted-Omega infimum exactly. Define

    R(t)=(t+sqrt(t^2+4))/2,   R(t)(R(t)-t)=1.

The preceding bounds imply every actual witness lies horizontally and above within D. For r>=t it contains all of A and has exact mass 2r(r+t)+2; strictness is exactly r<R(t). The mass supremum over such witnesses is 4R(t)^2. Witnesses with r<t have mass at most 4t^2+2<4R(t)^2. Consequently every proportional plan has source load at least

    J*(L)=integral_4^(L/4) t/[4R(t)^2] dt

on S. Conversely choose r_m(t)=t+(1-1/m)(R(t)-t), m>=1. These are actual strict witnesses; their denominators depend only on t. Every source load is at most the integral of the reciprocal denominator over Omega, with equality on S, and these maxima decrease to J*(L) by dominated convergence. Hence J* is exactly the infimum of the restricted-Omega maximum load. For nonempty Omega the infimum is not attained by a strict witness plan: every actual witness has a strictly smaller denominator, so every probability mixture has a strictly larger reciprocal denominator at each receiver; integrating the positive gap gives a strict excess on S.

An exact elementary antiderivative is

    J*(L)=[(1/4)log R +1/(16R^4)]_(R(4))^(R(L/4)).

Indeed t=R-1/R and dt=(1+1/R^2)dR transform the integrand to (1/4)(1/R-1/R^5)dR. This exact matching statement is only for the prescribed receiver set Omega, not for the whole halo.

## Actual computation and remaining issue

Python 3.15.0 ran `computation/round056/shared_background_proportional_independent.py` with exit 0: 5,855 exact Fraction geometry/constant checks, n=1,...,12, and 36 exact minimax/chain-rule checks, plus floating J0/J1/J* diagnostics. Randomness: none. Exact arithmetic assertions have no tolerance. The proof above supplies the universal quantifiers; sample checks alone do not.

The research consequence is a required change of allocation mechanism: even actual total capturing mass and arbitrary scale-adaptive witness mixtures cannot support a uniformly bounded *proportional* source load. A successful transport must use source-dependent nonproportional allocation or a different budget. The general weak-(1,1) upper bound and exact C2 remain unresolved.
