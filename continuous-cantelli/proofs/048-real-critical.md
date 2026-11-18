# Sextic germs with a critical real error, with no analyticity of that error

Let V(u)=c u⁶ H₀(u)² with H₀(0)=1 real analytic, c>0. Suppose v≥0 is merely continuous and |v(x₀+u)−V(u)|≤C exp(−d/|u|⁴). Put d_H=(2·2^{2/3}/3−1)/(4c). If d>d_H there is no exact graph. For H₀≡1 the same conclusion holds at d=d_H if the error is o(|u|^{12}exp(−d_H/|u|⁴)). This result supplements the sector theorem; it puts no assumptions on the sign, variation or analytic extension of the real remainder.

1. Define B_N=binom(3N/2,N/2)c^{N/4}E|Z|^{N/2} for N divisible by four. The analytic-reference proof r043-general-analytic-sextic-sector.md supplies on an infinite such subsequence local density coefficients A_N≥B_N exp[−o(N)]. Indeed the negative-b central real fold integrated on |b|∈[√(cN/2),√(cN/2)+1] contributes ≥ const C_N ρ^−N E|b|^{N/2}exp(−C N^{3/4}), where ρ=2/(3√3), C_N=4^{-N}binom(2N,N), and Stirling gives C_N ρ^−N comparable to binom(3N/2,N/2). With H₀'≠0 the positive-b possible negative part is smaller than the favorable real-fold term, as shown there; with H₀'=0 its integral is positive along the selected subsequence. For H₀≡1 the pure central Lagrange formula gives A_N≥c₀ B_N without a subexponential loss.

2. The Nth asymptotic coefficient of the actual-reference local density difference equals (1/N!)∫g(u)[p_v^{(N)}(−u)−p_V^{(N)}(−u)] du. For each fixed Taylor cutoff k use |u|>8|w| to dominate the Cauchy Taylor remainder by a constant times |w|^{k+1}|u|^{−k−4}exp(−κ/u⁴). On |u|≤8|w| the difference of the actual kernels is bounded by C|v−V||u|^{−9} using ∂_s p_s=p_s''/2, hence contributes o(|w|^k). The inner integrals of the fixed off-diagonal coefficient polynomials are negligible to every order. This needs no smoothness of v and avoids branch-conditioning the rough graph.

3. For ε>0 shrink the core so every intermediate variance is c u⁶ times a factor in [1−ε,1+ε]. The Hermite–Sobolev bound for k=N+2 gives
|D_N|≤ C_ε N^{5/4}/√(N!) [c(1−ε)]^{−(N+3)/2} ∫₀∞u^{−3N−9} exp[−(1/(4c(1+ε))+d)/u⁴]du.
An exterior annulus, if added back, contributes only geometrically, by the positive distance estimate for the fixed core. With A=1/(4c(1+ε))+d the integral is A^{−(3N/4+2)}Γ(3N/4+2)/4. Exact Stirling followed by ε↓0 gives
limsup_{4|N} (|D_N|/B_N)^{1/N} ≤ [ (32/27)/(1+4cd)^3 ]^{1/4}.
The right side is <1 at d>d_H, which contradicts exact analyticity by the reference lower bound.

4. At d=d_H, H₀≡1 and |δ|≤C|u|^{12}e^{−d_H/u⁴}, split first at |u|=N^{−1/6}. In the outer region the common lower variance is c'/N, giving a geometric error via the original Hermite bound. Inside, |δ|/V is exponentially small compared to N^{-2}; one can replace 1±ε in the last bound by 1±C/N² without exponential or polynomial loss. Including |u|^K replaces Γ(3N/4+2) by Γ(3N/4+2−K/4). At equality the result after division by B_N is O(N^{3−K/4}); in particular K>12 suffices. At K=12 the small-o envelope is uniform in the shrinking interval and supplies o(1), proving the endpoint statement.
