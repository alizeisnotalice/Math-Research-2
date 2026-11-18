# A spatial modulus forced at the maximum of every late time window

Use the normalized graph and the common extinction time σ of a component C of {f<1} from rounds 164 and 166. The argument applies to any real continuous solution of the heat equation on this cylinder, smooth in its interior, with q(σ,·)=0 and activity arbitrarily close to σ in every interior spatial interval. Fix J=[a,b] compactly contained in C, of length L>0, and let 0<h<σ−1. Write

    M(h) = max{|q(t,x)| : σ−h ≤ t ≤ σ, x∈J},
    Ω(h,r) = max{|q(t,x)−q(t,z)| : σ−h ≤ t ≤ σ, x,z∈J, |x−z|≤r}.

Activity gives M(h)>0. For any point (t*,x*) attaining M(h), put u=σ−t*∈(0,h] and ε=sign q(t*,x*). For 0<r≤L set p=6 barΦ(r/(2√u)), where barΦ denotes the standard normal upper tail. Whenever p<1 there exists z∈J with |z−x*|≤r and

    ε q(t*,z) ≤ [p/(1−p)] M(h).                                        (1)

Consequently for every h≤L²/16,

    Ω(h,√h) ≥ max_{x,z∈J, |x−z|≤√u}|q(t*,x)−q(t*,z)| > M(h)/5.         (2)

In particular for each fixed J the slices q(t*,·)/M(h) have no common modulus as h decreases to zero, regardless of the rate at which M(h) vanishes. For r(h) tending to zero with r(h)/√h tending to infinity, (1) guarantees an opposite sign or a value at most O(exp(−r(h)²/(8h))) times the peak in the peak's signed direction within distance r(h) at that very time. Each conclusion holds for every choice of late-window maximizer, which may lie at a spatial endpoint or at the first time of the window.

## Proof

The boundary traces used in R166 determine their full Dirichlet heat evolution. In particular, let W be Brownian motion started from y∈(a,b) and let τ be its first exit from (a,b). Since q solves q_t=q_xx/2 on a neighborhood of [t*,σ]×J, the stopped process q(σ−s,W_s), s≤u∧τ, is a bounded martingale. Hence

    0 = E_y[q(t*,W_u) 1_{τ>u} + q(σ−τ,W_τ) 1_{τ≤u}].                  (3)

At an extinction endpoint σ=H the same identity follows by continuity from σ−e and e↓0 with starting time t* fixed. Neither the boundary signs nor the values of q outside J enter (3).

Choose y as the projection of x* onto [a+r/2,b−r/2]. Thus |y−x*|≤r/2. The reflection principle and a union bound give P_y(τ≤u)≤4 barΦ(r/(2√u)). Also P_y(|W_u−y|>r/2)≤2 barΦ(r/(2√u)). On the complementary good event in (3) the time is t*, the spatial point lies in J within r of x*, and the probability of this event is at least 1−p. If m is the minimum of εq(t*,·) on this part of J, (3) implies 0≥m(1−p)−p M(h) whenever m≥0. If m<0 the claimed inequality holds already. This proves (1).

Take r=4√u≤L. Mills' bound gives 6 barΦ(2)<3 exp(−2)/√(2π)<1/6 (for instance e²>7.3 and π>3.1). Thus εq(t*,z)<M(h)/5 for some such z. Divide the interval from x* to z into four equal parts. The triangle inequality yields (2).

## Interior bound and higher derivatives

If K⊂J has distance d>0 from the endpoints, applying (3) at each x∈K yields

    M_K(h) ≤ E ceil(|Z|) Ω(h,√h) + 8 barΦ(d/√h) M(h),  Z~N(0,1).

Indeed on paths surviving to u the displacement has modulus cost bounded by ceil(|W_u−x|/√h) Ω(h,√h); on paths hitting the boundary the two values being compared each have absolute value at most M(h). In particular, a proposed common modulus for the normalized late-window slices forces strict relative interior decay. Estimate (2) excludes such a common modulus even on the closed interval J itself.

The mean value theorem in (2) gives ||q_x(t*,·)||_{∞,J}>M(h)/(5√u). More generally, for any integer n≥2 set ℓ=20(n−1)²√u. If ℓ≤L/2, then

    ||∂_x^n q(t*,·)||_{∞,J} ≥ n! M(h) / ℓ^n.                           (4)

At a point of the preceding derivative bound choose a one-sided interval of length ℓ inside J and let P be the degree n−1 Taylor polynomial there about that point. Taylor's formula gives ||P||≤M(h)+ℓ^n||∂_x^n q||/n!. The classical Markov polynomial inequality on that interval gives ||P'||≤2(n−1)²||P||/ℓ. Since P' at the Taylor center equals q_x, (4) follows. By the heat equation even orders also quantify higher temporal derivatives attained at some point of J at the very same t*. The spatial points can depend on the derivative order.

## Scope and check

The certificate concerns the signed density on compact subintervals of a low component, normalized by its own actual window amplitude. The graph f itself has no normalization by M(h), and these estimates alone do not prove positivity of the signed flow or discontinuity of f. Absolute spatial increments of q may be arbitrarily small while the normalized bound holds. The proof of the constant is analytic; a standard-library Python 3.15.0 evaluation gives 6 barΦ(2)=0.1365007916890752 and [1−p/(1−p)]/4=0.21048033675789873.
