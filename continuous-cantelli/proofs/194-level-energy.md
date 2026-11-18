# Integer heat derivatives averaged over distance levels

Let g=p₁ and let E be a compact Ahlfors d-regular subset of the line, 0<d<1, Q=1−d. Suppose an exact graph has terminal variance deficit v≥0 equal to zero on E, with ∫p_{v(x)+s}(a−x)g(x)dx=p_{H+s}(a) for all a and s>0, H>0. Other zero sets and vanishing deficits at infinity are allowed. The following theorem supplies all inverse power budgets when the level counting measures become uniform on the heat length scale. In particular it applies to equal-contraction strongly separated Cantor sets with an asymptotically radial deficit.

For a.e. sufficiently small r>0, let ν_r count the finitely many points with distance r to E, and let N(r) be their number. In the line the measure dx disintegrates as ν_r dr on the complement of E, by splitting into the two halves of each complementary gap and the exterior rays. Gap midpoints form a countable exceptional set of radii. Assume N(r) is comparable to r^(−d). Let μ be a d-regular probability supported on E. On a fixed collar suppose there is a positive V(r) tending to zero with V(r)/r²→∞. Assume

    sup_{||φ||∞≤1, Lip(φ)≤1/√V(r)} |∫φ d(ν_r/N(r)−μ)| → 0.                 (1)

Suppose v(x)=V(r)(1+o(1)) uniformly on levels as r↓0. More generally it suffices that v≥cV on the collar and

    N(r)⁻¹∫ min(1,|v(x)/V(r)−1|)dν_r(x) → 0.                              (2)

Under these hypotheses ∫₀^ε r^(−d)V(r)^(−Q/2−n)dr is finite on a sufficiently small interval for every fixed nonnegative integer n. Hence no power law V comparable to r^p for 0<p<2 is possible. For the last implication choose n with p(Q/2+n)≥Q. A different ε for each n suffices.

**Positive mean on each level.** Put σ=gμ, τ=Q/2+n and T_b=(−∂_b)^n P_bσ, where P_b convolves with p_b. Gaussian derivative bounds and the Ahlfors upper estimate on dyadic annuli give

    ||T_b||∞ ≤ C_n b^(−τ),  ||∂_xT_b||∞ ≤ C_n b^(−τ−1/2),
    ||∂_bT_b||∞ ≤ C_n b^(−τ−1)                                               (0<b<1).

The mass J(b)=∫P_bσ dσ is between positive constant multiples of b^(−Q/2), by lower regularity on a √b-ball and the same annulus estimate. By Plancherel for Gaussian regularization, J(b)=(2π)⁻¹∫exp(−bξ²/2)|σ̂(ξ)|²dξ. Fix a large B so e^(1/2)J(Bb)≤J(b)/2 for small b, and discard |ξ|≤1/√(Bb) at cost at most e^(1/2)J(Bb). On the remainder (ξ²/2)^n≥(2Bb)^(−n). Consequently

    ∫T_b dσ=(2π)⁻¹∫(ξ²/2)^n exp(−bξ²/2)|σ̂(ξ)|²dξ ≥ c_n b^(−τ).             (3)

With b=V(r)+s small, the test φ=b^τ gT_b divided by a constant meets (1), because √V(r)≤√b. Thus replacing ν_r/N(r) by μ costs o(b^(−τ)) in ∫gT_b, uniformly in s. To replace V(r) by the actual v(x), use the time derivative bound along the segment from b to v(x)+s whenever |v/V−1|≤1/2. At the other points the size estimate suffices using v≥cV; if v+s is bounded away from zero its kernel is bounded. The resulting loss is at most C b^(−τ) times the mean in (2). Together with (3) this proves

    ∫g(x) T_{v(x)+s}(x)dν_r(x) ≥ c'_n N(r)(V(r)+s)^(−τ)

uniformly for all sufficiently small r,s>0. Integrate in r. The full signed integral on the line equals ∫(−∂_s)^n p_{H+s}dσ by exactness, bounded uniformly as s↓0. For fixed s the identity can be differentiated using domination by a constant depending on s. The region at fixed positive distance from E contributes a uniformly bounded absolute amount, since short-age heat derivatives there decay exponentially in inverse age and all positive ages away from zero give bounded kernels. E has zero Lebesgue measure. Fatou therefore gives the claimed budget.

**Verification of (1) for homogeneous Cantor sets.** Suppose E comes from m≥2 similarities of common ratio ρ, with strictly separated first-level hulls, and use the equal-weight probability μ. Orientations may reverse. Every complementary bounded gap is an image of one of the m−1 first-level gaps under a finite word. The number whose length exceeds 2r is comparable to r^(−d) by summing m^j for j up to log(1/r)/log(1/ρ)+O(1), since there are finitely many positive basic gap lengths. Hence N(r) is comparable to r^(−d). Choose a depth k with hull diameter ℓ comparable to √(r√V(r)), so r≪ℓ≪√V(r). Each depth-k hull includes precisely the same number of level-r points coming from gaps at depth k and deeper. The remaining points come from the O(m^k) ancestral gaps and exterior, costing O((r/ℓ)^d) of the normalized count. The equal per-hull portions can be coupled to μ, matching corresponding hulls, at distance at most ℓ. Thus the discrepancy in (1) is at most C[(r/ℓ)^d+ℓ/√V(r)], which tends to zero. The same proof gives (1) for any finite hull length after changing constants. The argument applies at all radii except the countable gap-midpoint radii, enough for disintegration.
