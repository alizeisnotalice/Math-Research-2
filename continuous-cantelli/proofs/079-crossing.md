# Noise sensitivity of the sign of the variance deviation

Let X,G be independent standard normal variables. Write a≥0 for the conditional variance, m=E a(X), h=a−m and ε=‖h‖∞, taking essential suprema under standard Gaussian measure γ. For a nonconstant exact graph the project's boundedness theorem gives 0<ε<∞. Put L=E|h(X)|, θ=arctan √ε, ρ=cos θ and X'=ρX+sin(θ)G. Write φ for the standard Gaussian density. Define c=2φ(1)log 2 and w(dx)=|h(x)|γ(dx)/L.

**Theorem.** For every measurable B with values ±1 that agrees with sign h where h≠0, the weighted probability of B(X')≠B(X), with X drawn from w and independent G, is at least (1−c)/2. More generally, for *any* two-valued step classifier B with cuts z₁,…,z_J, let e_B=w{B≠sign h}. Then

    e_B + ε θ/(π L) ∑_{j=1}^J exp(−z_j²/2) ≥ (1−c)/2.                 (1)

One can take the infimum over all finite classifiers. In particular an exact graph whose mean-level sign is described by J cuts satisfies J≥π(1−c)L/(2ε arctan√ε). In a family with L/ε≥q>0 and ε↓0, at least a constant times q/√ε crossings are necessary. Formula (1) also counts the price of disregarding low-amplitude sign oscillations or placing all cuts in the Gaussian tails.

**Global estimate.** For u=ε, the two variables X+√(u+h(X))G and X+√u G have the same N(0,1+u) law. This follows directly by multiplying their characteristic functions by the scalar that shifts m to u, using u+h≥0. For any bounded real k with ‖k‖∞≤1, the heat semigroup P_v gives

    E[h(X)P_u k(X)] = −(1/2) E[h(X)² ∫₀¹(1−s)(P_{u+sh(X)}k)″(X)ds].   (2)

For smooth k, compare the expectations of a twice-antiderivative of k at the two endpoints and apply the integral Taylor formula in the variance, using ∂_v P_v=∂²P_v/2. The antiderivative grows at most quadratically. For bounded measurable k first convolve with a Gaussian and send its variance to zero. The left side converges because u>0. For s<1 all heat times are positive, and the second derivative has absolute value at most E|G²−1|/v=4φ(1)/v. Since v≥u(1−s), domination suffices even if the final variance vanishes.

On {h>0} the ratio h/(u+sh)≤1/(1+s), and on {h<0}, |h|/(u+sh)≤1/(1−s). Since both signed wings have absolute mean L/2, (2) yields

    |E h P_u k| ≤ [4φ(1)/2] (L/2) [∫₀¹(1−s)/(1+s)ds +1] = c L.     (3)

Choose k(y)=B(ρy). Then E[h(X)B(X')]≤c L. Since E[h(X)B(X)]=L(1−2e_B), their difference is at most 2E[|h(X)|1_{B(X')≠B(X)}]. This proves the stated noise sensitivity for a perfect sign classifier. For an imperfect one it gives e_B plus its weighted flip probability at least (1−c)/2.

**Crossing count.** Along X(t)=X cos t+G sin t, 0≤t≤θ, a classifier flip entails crossing some cut z. Conditional on radius R=(X²+G²)^(1/2), the initial angle is uniform. The expected number of hits of z is θ/π if R>|z| and zero if R<|z|. The equality event is null. Averaging the Rayleigh tail gives (θ/π)exp(−z²/2). Sum over cuts and use |h|≤ε to obtain (1). The same argument works for countably many cuts if every flip crosses a cut and the displayed sum is finite.

For a continuous h there is a corresponding geometric estimate requiring no regularity of the zero set Z={h=0}. The very same weighted coupling has Pr{segment [X,X'] meets Z}≥(1−c)/2, since opposite signs (or arrival at a zero) require an intermediate zero when starting with h≠0. All statements are necessary conditions. A number of tiny oscillations, or sufficiently rapid oscillations in the main part of the input distribution, is consistent with them.

**Whole-line variation.** In fact (3), after replacing k(y) by any test of ρy, and Gaussian reversibility give the global contraction ‖T_ρ h‖_{L¹(γ)}≤c‖h‖_{L¹(γ)}, where T_ρ is the stationary OU operator with correlation ρ. Consequently any exact graph of locally bounded variation satisfies

    ∫ φ(x) d|Dh|(x) ≥ √(π/2) (1−c) L / arctan√ε.

For smooth h integrate the absolute derivative along the same rotating Gaussian line. At each angle its velocity is an independent standard normal, so E|h(X')−h(X)|≤θ√(2/π)∫|h'|dγ. For a bounded locally BV function directly integrate expected level hits against the ordinary variation measure |Dh|: except for a null event at the endpoints, the increment between the two Gaussian values is bounded by the variation on the segment, and every level strictly between the values is crossed by the rotation. The hit bound above and Tonelli give the stated inequality, including the case of infinite weighted variation. The opposite inequality E|h(X')−h(X)|≥(1−c)L follows by conditioning on X and taking the L¹ norm.
