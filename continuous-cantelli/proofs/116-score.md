# One-sided score concentration and partial escape

Let γ be standard Gaussian measure, let 0≤f≤H be measurable, and assume that every real tilt ν_t(dx)=exp(tx−t²f(x)/2)γ(dx) has total mass one. Set W_t=X−tf(X). Continuity and reciprocal-center ordering enter only in the second theorem. Write p_v for the centered Gaussian density of variance v. All assertions have reflected versions with x and t replaced by −x and −t.

## Concentration on a half-line

For t,R>0 and u>0 one has

    ν_t(X≤R, W_t≥u) ≤ exp[−tu²/{2(R−u)}]       (0<u<R),
    ν_t(X≤R, W_t≤−u) ≤ exp[−tu²/{2(R+u)}].

The first event has zero probability for u≥R since W_t≤X and ν_t has no atoms. Substitution of f=(X−W_t)/t into the exact probability change ν_t→ν_{t+s} gives, for all real s,

    E_{ν_t} exp[(s+s²/(2t))W_t−s²X/(2t)] = 1.

For the first inequality take s=tu/(R−u). The exponent on the event is at least su−s²(R−u)/(2t). For the second take s=−tu/(R+u), for which the coefficient of W_t is negative. The lower bound on the exponent is −su−s²(R+u)/(2t). These choices give the stated exponents exactly. No comparison between inputs at different positions is needed.

For every integer k≥1 and R,t>0 this implies

    E_{ν_t}[1_{X≤R}|W_t|^k] ≤ C_k [(R/t)^{k/2}+t^{−k}],
    E_{ν_t}[1_{X≤R} f^k]    ≤ C_k [(R/t)^k+t^{−2k}].

Indeed the two-sided bound is at most 2 exp[−tu²/(2(R+u))]. Split the moment integral at u=R, using Gaussian decay for u≤R and exponential decay for u≥R. For f use tf≤R+|W_t|, and absorb (R/t³)^{k/2} into (R/t)^k+t^{−2k} by AM–GM. In particular, for every deterministic R(t)=o(t) with R(t)>0, the unnormalized score on X≤R(t) tends to zero in every fixed positive integer moment. If the mass of this event stays above a positive constant, the same holds for its conditional law. More quantitatively, any fixed collection of sublinear sources carrying appreciable mass has tf(X)−X=O_P(√(R/t)+1/t) there.

For a fixed measurable set E⊂(−∞,R] these inequalities imply

    |(d/dt)^k ν_t(E)| ≤ C_k (√(R/t)+1/t)^k.

Differentiating the tilt k times gives a finite sum of monomials W_t^{k−2j} f^j of total weighted degree k. Hölder's inequality with the preceding moment bounds gives the conclusion, since the norm scale for f is bounded by C(√(R/t)+1/t)^2. Differentiation on compact parameter intervals is justified by f≤H and Gaussian tails. Integrating the absolute first derivative of the density also gives the total variation bound for restrictions

    ∫_{x≤R}|dν_{t+h}−dν_t| ≤ C |h|(√(R/t)+1/t),       |h|≤t/2.

R may depend on the base t in this comparison, holding its value fixed along the integration segment. In particular every fixed compact-source mass changes by o(1) over bounded parameter shifts.

## An inverse-square law under partial escape

Assume now that f is continuous, eventually positive on the positive ray, and m(x)=x/f(x) is eventually nondecreasing. Suppose some θ∈(0,1] satisfies

    limsup_{t→∞} ν_t([0,R]) ≤ 1−θ                 for every fixed R>0.    (P)

Then f(x)=1+O_θ(x^{−2}) as x→+∞. More precisely limsup x²|f(x)−1|≤C_θ. Fix R₀ sufficiently large that f(x)≥15/16 for x≥R₀, which the argument below ensures is possible. Write c(t)=ν_t([0,R₀]) and a(x)=x²(f(x)−1). For every compactly supported integrable ψ,

    ∫ ψ(u)[exp(−a(t+u)/2)−(1−c(t))]du → 0.

If a is uniformly continuous on the positive tail, this implies a(t)+2 log(1−c(t))→0. Thus a compact mass limit c forces a translated average 1−c and, with uniform continuity, pointwise limit −2 log(1−c). Under (P), c(t)≤1−θ/2 eventually. The choice of R₀ among all sufficiently large fixed cutoffs changes c(t) by o(1).

Here are the precise modifications to the companion proof `r101-monotone-centers.md`. Its fixed-fraction low-point contradiction uses only normalization and ordering, so liminf f≥1. At a fixed-fraction high point its macroscopic interval around t=(x+m(x))/2 contributes o(1). So does the region outside that interval on X≥R for a sufficiently large fixed R, by the shifted Gaussian estimate using liminf f≥1. The negative half-line has mass O(1/t) by f≥0. The rest has limsup ≤1−θ by (P), a contradiction to total mass one. Hence f→1. For any sufficiently large fixed R₀ the interval [R₀,t/2] has exponentially small mass since the multiplier exp[−t²(f−1)/2]≤exp(t²/32) and P_{N(t,1)}(X≤t/2)≤exp(−t²/8).

Use the companion's Gaussian and bandlimited smoothing at large t, with K_t=[t/2,2t], η=sup_{K_t}|f−1|, and r=√(1+t²η)/t. Before removing retained exterior sources, the smoothed density converges uniformly to p₁. This uses the companion's exponentially accurate complex-tilt Fourier estimate and convergence of the reference convolution to p₁. Fix small ε>0. Retained sources on [0,R₀] with f≥ε have mass at most exp(tR₀−εt²/2) and their smoothed density is at most Ct times that mass. Compact sources with f<ε have independent additional Gaussian variance at least 1−ε, and their entire normalized density is bounded everywhere by c(t)/√(2π(1−ε)). For X<0, split at f=ε; the low set has mass O(1/t) and variance at least 1−ε, the high set mass ≤exp(−εt²/2). For X∈[R₀,t/2] and X>2t, the above shifted Gaussian estimates give exponentially small mass and the bandlimited noise bounds density by Ct times this mass. Select ε and a fixed ℓ>0 small enough that

    (1−θ/2)/√(2π(1−ε)) < inf_{|w|≤ℓ} p₁(w).

After removing the exterior sources the local smoothed density still has a strictly positive lower bound on [−ℓ,ℓ] and a universal upper bound everywhere. The companion's strip argument therefore gives μ_t(J)≤C(|J|+r) globally and c_θ|J|≤μ_t(J)+Cr for J⊂[−ℓ,ℓ], where μ_t is the unsmoothed W_t law restricted to K_t. Its high-point test uses s=min((x−m(x))/8,ℓ/2); the low-point test uses only the global upper bound. The same neighboring-shell recurrence then proves the asserted inverse-square law, with limiting bound depending only on θ.

Finally c(t+s)−c(t)→0 for fixed s by the first theorem. Let t_j→∞ with c(t_j)→c. On each bounded u-interval any local weak-star subsequential limit H of exp(−a(t_j+u)/2) exists after extraction, because a is uniformly bounded on the positive tail. One can extend the multiplier boundedly to the rest of the line for this extraction. At t=t_j+s, remove X≤R₀ with mass c(t)+o(1) (the negative part is O(1/t)); remove R₀<X<t/2 with exponentially small mass. On X≥t/2 the multiplier exp[−t²(f−1)/2] is uniformly bounded by the inverse-square law, and on |X−t|≤K it differs uniformly from exp(−a(X)/2) by o(1) for each fixed K. The shifted Gaussian formula and then K→∞ give (p₁*H)(s)=1−c for every s. Gaussian convolution is injective on bounded distributions, by Fourier testing with compactly supported smooth functions and dividing locally by the nonvanishing multiplier. Thus H≡1−c. Compactness gives the asserted varying-center integral statement, and local uniform compactness of translates for uniformly continuous a gives the pointwise consequence. On both wings these tail results imply an attained maximum >1 for every nonconstant graph; as in the companion, round 23 rules out such graphs if some attained maximum has a nonconstant analytic germ. A continuous rough maximum is still possible under these necessary conditions.
