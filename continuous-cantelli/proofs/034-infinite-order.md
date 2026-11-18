# Bilateral infinite-order flatness forces a zero neighborhood

Suppose the exact variance identity of round 31 holds, so the positive potential of the means |u|<R, in coordinates centered at an attained zero, is

    S_R(y)=∫_{|u|<R} ρ(x₀+u) H_{v(x₀+u)}(y−u)du,
    H_s(z)=½∫₀ˢ p_t(z)dt.

Assume a local positive-order Hölder bound, and that U(b)=sup_{|u|≤b}v(x₀+u)=O_N(b^N) as b↓0 for every fixed N. Then v vanishes throughout some neighborhood of x₀. In particular this holds for every bilateral C∞ zero with vanishing full Taylor jet. No prescribed stretched-exponential envelope is needed.

Write v(u) for the translated field, g(u)=ρ(x₀+u), and q(u)=u²/v(u) whenever v(u)>0. All integrands below are taken as zero when v(u)=0. Uniform analytic continuation and the bound B of S_R on |z|≤3R/8 come from the round-31 exterior argument; this never requires pointwise continuation of H across its mean. We prove a new estimate on the Taylor coefficients of that continuation.

**Off-diagonal formula for the complete jet.** For u≠0, H_{v(u)}(z−u) is holomorphic when |z|≤|u|/4. In that disc Re((z−u)²)≥u²/2, so its maximum is ≤C √v(u) exp(−q(u)/4). By Cauchy's formula its coefficient h_n(u) of z^n satisfies

    |h_n(u)| ≤ C (4/|u|)^n exp(−q(u)/4).

Shrink the initial neighborhood so that v(u)≤K u⁴ and R≤1, giving q(u)≥1/(Ku²). For each n this estimate is integrable down to zero. To verify that a_n=[z^n]S_R(z) equals ∫g(u)h_n(u)du, fix any Taylor cutoff k independently of y. For real y→0, on |u|>8|y| the remainder through degree k is ≤2C(4|y|/|u|)^(k+1)exp(−q/4), whose integral is O_k(|y|^(k+1)). On |u|≤8|y| the original nonnegative potential integral is ≤||g||∞ U(8|y|)/2, since H increases with its variance and ∫_R H_s=s/2. It is O_k(|y|^(k+1)). The integrals of each of the first k+1 Taylor terms over this small set are also smaller than any fixed power of |y| by the above Gaussian-in-1/u bound. This proves the full off-diagonal formula by uniqueness of the real asymptotic Taylor coefficients of a holomorphic function. The argument uses infinite-order flatness only to justify this identity for every finite k; it does not need uniform control of its k-dependent remainders.

**A quantitative long initial segment of coefficients.** At the hitting level m=m_j of round 31 put R=r_{j−1}; then v≤2m throughout |u|<R. Hence

    q(u)≥max{u²/(2m),1/(Ku²)},
    q(u)≥c₀ m^(−1/2),       q(u)≥c₁ m^(−1/4)/|u|.

The first derived estimate is the geometric mean of the original two bounds, and the second is the geometric mean of that estimate and 1/(Ku²). Since R≤1, log(4R/|u|)≤4/|u|. Fix c₂>0 small enough that n≤N=floor(c₂m^(−1/4)) entails n log(4R/|u|)≤q(u)/8. For 0≤n≤N one has uniformly

    R^n |a_n| ≤ C∫_{|u|<R} exp(−q(u)/8)du ≤ C exp(−c m^(−1/2)).

Here the full coefficient formula was proved before n was chosen as a function of the fixed core. For n>N use the independent holomorphic bound |a_n|≤B(8/(3R))^n. Consequently summing the power series gives

    sup_{|y|≤R/32}|S_R(y)| ≤ C exp(−c m^(−1/4)).                  (1)

The low coefficients form a bounded geometric sum after scaling by R/32, and the high terms have ratio 1/12.

**Conclusion by first-hitting propagation.** If no zero neighborhood exists, define m_j=2^(−j)m₀, with 0<r_j<r_{j−1} the first bilateral envelope hits at those levels, so r_j↓0. At a hitting point x_j=±r_j the Hölder lower estimate of round 31 gives S_R(x_j)≥c m_j^(1/α+1/2). The round-31 propagation chain from |y|≤R/32 to x_j gives

    c m_j^(1/α+1/2) ≤ C exp[−c (1−r_j/R)^p m_j^(−1/4)].

Thus 1−r_j/r_{j−1} ≤ C[(j+1)m_j^(1/4)]^(1/p), a summable series for fixed p>0. The infinite product of the ratios is therefore strictly positive, whereas r_j→0. This contradiction proves the claim. The C∞ case supplies local Lipschitz continuity and all flat bounds by Taylor's theorem. If the graph only obeys every flat upper bound without local Hölder or a comparable mass bound at its hitting points, the last lower-estimate step needs another hypothesis.
