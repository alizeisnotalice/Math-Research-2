# Independent check and improved OU sign-cut constant

Let b=a−m, ε=∥b∥∞>0, u=ε. Cancelling or adding the independent variance u−m shows that X+√(u+b)G has law N(0,1+u), since u+b≥0 even when u<m. For a bounded Borel k with ∥k∥∞≤1 the endpoint Taylor identity and ∥∂² P_T k∥∞≤E|Z²−1|/T imply

 |E[b P_u k(X)]| ≤ κ E[(u+b)log(1+b/u)−b],    κ=E|Z²−1|/2=2g(1).

The integrand is b²∫₀¹(1−s)/(u+sb)ds with continuous endpoint convention 0log0=0. The identity uses a twice primitive Φ with Φ″=k, by first smoothing; bounded measurable k may be obtained by convolution because the Gaussian baseline has positive variance. Multiplication by 1−s in the Taylor remainder cancels the possible denominator singularity at b=−u. The convex integrand lies below the chord on each of [−ε,0],[0,ε], hence is ≤|b| on the negative side and ≤(2log2−1)|b| on the positive side. By centering, E b_+=E b_-=E|b|/2. Thus |E[bP_uk]|≤κ log(2)E|b|, sharpening the direct κ bound.

Let S denote any ±1 completion of the sign of b on {b=0}. For ρ=(1+ε)^(−1/2), X′=ρ(X+√ε G) is standard normal correlated by ρ with X. Take k(y)=S(ρy). Then
E[|b(X)|1_{S(X)≠S(X′)}] ≥ (1−κ log2) E|b|/2.
The ±1 convention matters when the zero set has positive Gaussian mass. If S switches at finitely many cuts z_i, every mismatch crosses one cut. The probability of crossing z for two ρ-correlated standard normals is π^−1∫_ρ¹ exp(−z²/(1+r))/√(1−r²)dr. This follows, for instance, by integrating the correlation derivative of the bivariate Gaussian lower quadrant at (z,z); ∂_r p_r(x,y)=∂_x∂_y p_r(x,y). The crossing probability is ≤g(z)/g(0) arccos(ρ)/π. Consequently

 Σ_i exp(−z_i²/2) ≥ π(1−κ log2)E|b| / [2 ε arctan(√ε)].

If S is only a proposed sign classifier rather than a completion of sign b, mismatches with that classifier must be included separately. The inequality quoted here concerns an actual sign completion.

**Hermite energy and the Gaussian Dirichlet form.** The supremum over k first gives ∥T_ρ b∥₁≤κ log(2)∥b∥₁ by reversibility of the OU kernel. Independently, k(y)=b(ρy)/ε and the pointwise F(b)≤b²/ε yield
0≤〈b,T_ρ b〉≤κ Var(b).
For a free t≥1 take u=tε, ρ_t=(1+tε)^(−1/2). Maximizing ∫₀¹(1−s)/(tε+sb)ds over |b|≤ε at b=−ε gives
Σ_{j≥2} ρ_t^j ∥Π_j b∥² ≤ q(t) Var(b),
q(t)=κ[1+(t−1)log(1−1/t)], q(1)=κ.
Here Π_j denotes the Gaussian Hermite projection of degree j; the constant and linear coefficients vanish. This follows from the usual expansion of OU-correlated Gaussian characters, or by applying their exponential generating function to the normalized Hermite basis. Therefore the energy fraction in degrees at most J is ≤q(t)(1+tε)^{J/2}. With probability weights ∥Π_j b∥²/Var(b), Jensen applied to the strictly convex exponential gives
Σ j∥Π_j b∥²/Var(b) ≥ 2 log(1/q(t))/log(1+tε).
The left hand side equals E|b′|²/Var(b) in Gaussian Sobolev space, and is interpreted as +∞ for an infinite-energy graph. Thus small-amplitude candidates must carry a fixed positive proportion of their L² energy in Hermite degrees of order at least 1/ε. R001's fixed-finite-dimensional isolation is logically different and does not give this scale-free spectral conclusion. The estimates are derived only for exact candidates; true nonconstant continuous candidates are not assumed to exist.
