# Global two-sided convex order and Stein layer flux

Let X,G be independent standard normals, A=a(X)≥0 bounded, m=EA, b=A−m, Y_s=X+√(m+s b)G, 0≤s≤1. Suppose Y_1 has the Gaussian law of Y_0. Let P_v be the heat semigroup of variance v (P_0 the identity) and g the standard Gaussian density. If Var A>0 then m>0.

For each real strike k set ψ_k(y)=(y−k)_+. There are two strict sign identities

  −E[bψ_k(Y_0)] = (1/2) E[b²∫₀¹(1−s)p_{m+sb}(k−X)ds] > 0,
   E[bψ_k(Y_1)] = (1/2) E[b²∫₀¹s p_{m+sb}(k−X)ds] > 0.                 (1)

There is no extremum or differentiability assumption. Both right sides have integral Var(A)/4 with respect to k. At any s<1 the interpolated variance is ≥(1−s)m>0; the 1/√(1−s) bound permits dominated convergence even at a zero terminal variance. To prove (1) apply the direct and reverse integral Taylor formulas in s to E[P_{m+sb}Φ(X)] with Φ(y)=(y−k)_+³/6. Its endpoint values agree, Φ″=ψ_k, Φ⁽⁴⁾=δ_k, and ∂_v P_v=∂²_y P_v/2. Gaussian regularization of Φ, followed by dominated convergence, justifies all endpoint formulas. The same argument displays the full positive interpolation gap, since the second s derivative equals E[b²p_{m+sb}(k−X)]/4>0 for s<1.

The marked finite measures of outputs with weights b_+ and b_- have the same mass and first moment at either endpoint: E b=E[Xb]=0 (the second equality follows by the third Gaussian moment). After dividing by E b_+, at the independent-noise endpoint the positive mark is smaller in convex order than the negative one; at the actual endpoint it is larger. Formula (1) is strict for every strike, hence also for nonaffine convex tests with finite required expectations by the call representation. Fixed variance inf a may first be subtracted from a by cancelling the common Gaussian characteristic factor, even if the infimum is not attained. This gives the same tests at the smallest permissible constant-noise endpoint m−inf a.

**Derivative-free Stein cut.** Set H_k=1_{(k,∞)}, F_k(x,t)=∫₀ᵗ P_s H_k(x) ds. Then

 E[X F_k(X,A) − ∫₀ᴬ p_s(k−X)ds]
       = −E[b²∫₀¹ s p_{m+sb}(k−X)ds] < 0.                              (2)

Indeed for φ=(y−k)_+²/2, the Gaussian target Stein identity and conditional Gaussian integration by parts give
0=E[X P_Aφ(X)+(A−1−m)P_Aφ′(X)].
The heat differences P_Aφ−φ=F_k/2 and P_Aφ′−φ′=½∫₀ᴬ p_s(k−X)ds, combined with E[Xφ(X)−φ′(X)]=0, imply that the left side of (2) equals −2E[bψ_k(Y_1)]. Both terms are finite since |F_k|≤sup a and the second integral is bounded by C√sup a. In particular the formula applies without assuming bounded variation. Its left side equals the integral in s of
∫_{a(x)>s} g(x)[xP_sH_k(x)−p_s(k−x)]dx
over 0<s<sup a.

If a is locally absolutely continuous and ∫g|a′|<∞ this layer flux is ∫g(x) P_{a(x)}H_k(x) a′(x) dx. This follows by integration by parts for F_k(x,a(x)); one can first smooth H_k with a small independent Gaussian, for which the chain rule applies without a degeneracy at a=0. Dominated convergence removes the smoothing since a′=0 almost everywhere on {a=0}. The same conclusion for a continuous locally BV graph holds with its Stieltjes variation da, assuming ∫g|da|<∞. In this case smoothing converges pointwise away from the single x=k whenever a=0; the continuous graph has a nonatomic Stieltjes measure. Positive and negative derivative measures g(da)_± have equal mass D by Gaussian integration by parts and E[Xa]=0. The output from the positive derivative measure is strictly smaller in first-order stochastic order than the one from the negative derivative measure, by (2) for all k. If ∫(1+|x|)g|da|<∞, the W₁ distance of their normalized outputs equals the integral of their tail differences, exactly Var(A)/(2D) by (1). No unconditional W₁ assertion is needed when either first moment is infinite.

**Scope and archived overlap.** The R001 quartic covariance is the integral-in-strike version of either side of (1). The R002 occupation identity uses the unmarked signed time residual, whereas (1) keeps the actual variance marks and (2) orders derivative-weighted outputs. Neither sign implies rigidity for all continuous graphs: even the exact non-graph three-layer coupling in R005 satisfies (1), whose proof uses only the two Gaussian endpoint laws. The BV derivative interpretation additionally requires a graph.
