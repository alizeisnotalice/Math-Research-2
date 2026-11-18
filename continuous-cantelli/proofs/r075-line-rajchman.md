# A local Rajchman test for equal-ratio Cantor variance zeros on the line

Let K⊂R be the attractor of m≥2 similarities with a common absolute contraction ratio ρ∈(0,1). Assume the compact first-generation images are pairwise disjoint, let d=log(m)/log(1/ρ)<1, and let ν be the equal-weight self-similar probability. Reflections and arbitrary real translations and ratios are permitted. Let Y=X+√v(X)G, G independent standard normal. On a neighborhood of K the law of X has a continuous strictly positive density g, and v is continuous, vanishes on K and is positive off K there. Assume uniformly as D(x)=dist(x,K)↓0 that

    v(x)=D(x)^α[a(x)+o(1)],       0<α<2,

where a is continuous and positive near K. Outside this neighborhood the source can be any finite measure with any nonnegative measurable variance, including other variance zeros. Set β=2(1−d)/α. If the output is locally smooth near K, ν is Rajchman, meaning ∫exp(iξx)dν(x)→0 as |ξ|→∞ through arbitrary real frequencies. It suffices that every smooth test ψ compactly supported in this neighborhood satisfy

    |ξ|^β ∫ ψ(y)e^{iξy}dLaw(Y)(y) →0.                             (1)

It also suffices to have the assumptions on a neighborhood of a single point of K: pass to a fixed sufficiently deep cylinder wholly inside it. That compact cylinder is separated from the remaining cylinders at the same depth, and its equal-weight measure is a similarity image of ν, which preserves the Rajchman property.

**Heat localization.** Fix a bounded open U compactly inside the prescribed neighborhood, containing K, with g,a bounded above and below there and v comparable there to D^α. It may be chosen to be a sufficiently small fixed tube around K. Put σ=|ξ|^(−2/α) and Z=σ^(1−d)=|ξ|^(−β). The tube length satisfies |K_r|≍r^(1−d) near zero. In fact at level j with ρ^j≍r, the m^j sets of diameter O(r) cover K, and one point from each is separated from the others by ≳r: sets first differing at level h≤j are at least Δρ^(h−1) apart, where Δ is the minimum first-generation separation. Intervals of radius c r around these points then give the lower bound. Layer-cake integration gives ∫_U exp(−c (D/σ)^α)dx≍Z for every fixed c>0.

For ψ smooth compactly supported in U, the part of the output from input outside U is smooth on supp ψ, by positive separation and sup_{s>0}|∂_y^j p_s(y−x)|≤C_j |y−x|^(−j−1). Deterministic sources outside U vanish there. For the remainder Fourier inversion of ψ gives

    ∫ ψ(y)e^{iξy} dLaw(Y)(y)
        = ∫_U ψ(x)e^{iξx}e^(−ξ²v(x)/2) g(x)dx + O_ψ(Z/|ξ|)+O_{ψ,J}(|ξ|^(−J))

for every J. For a shift s with |s|≤|ξ|/2 the two Gaussian multipliers at ξ+s and ξ differ by at most C |s|/|ξ| exp(−ξ²v/16), by the mean value theorem. Integrate against the rapidly decreasing Fourier transform of ψ and use the tube bound. The complementary frequencies and the exterior input have rapid decay. Consequently (1) makes the right main integral o(Z). The same normalized cancellation holds for any family h_ξ of uniformly bounded equicontinuous tests supported in a common compact subset of U: use a finite uniform ε-net of fixed smooth tests and the local total mass bound O(Z), then ε↓0. No regularity of a or g beyond continuity is needed.

**Congruent tubes and continuous inverse weights.** Choose j=j(ξ) with ρ^j comparable to (σ/|ξ|)^(1/2); thus σ≪ρ^j≪|ξ|^(−1), by powers. For each level-j word w let K_w be its cylinder and T_w={dist(x,K_w)<δρ^j}, with δ>0 smaller than Δ/(3ρ). These tubes are disjoint and on each one D(x)=dist(x,K_w). They are congruent up to reflection, including when convex hulls of different cylinders overlap. On U outside their union the tilt mass is O(exp(−c(ρ^j/σ)^α))=o(Z). The union lies in an arbitrarily small fixed neighborhood of K eventually.

In any one tube define W_ξ(t)=∫_{T_w} exp(−ξ² t D(x)^α/2)dx. It is identical for every word w. For t on a compact positive interval including 1 and the range of a, W_ξ(t/2)/W_ξ(t)≤2^(1/α). Indeed |(K_w)_{qr}|≤q |(K_w)_r| for q≥1 on the real line: internal gaps of length l contribute min(l,2r), the two exterior ends contribute 2r, and the set itself contributes its own length. The inequality remains true after truncation at the fixed tube radius δρ^j. Layer-cake then gives the desired doubling. Convexity of log W_ξ implies

    0 ≤ −∂_t log W_ξ(t) ≤ (2/t) log[W_ξ(t/2)/W_ξ(t)] ≤ 2 log(2)/(α t).

Therefore ω_ξ(t)=W_ξ(t)/W_ξ(1) and its reciprocal are uniformly bounded and equicontinuous on that compact range. The frozen tube sum m^j W_ξ(1) is comparable to Z by the tube bound and the negligible exterior for the reference field D^α.

Take χ smooth compactly supported in U and identically one near K. The tests h_ξ(x)=χ(x)/[g(x)ω_ξ(a(x))] are admissible in the localized cancellation estimate. For arbitrary x_w∈K_w, the total actual unweighted exponential mass on T_w is [1+o(1)]W_ξ(a(x_w)), uniformly in w: throughout this shrinking tube v is bracketed between D^α[a(x_w)±ε_ξ] with ε_ξ→0, and the log derivative bound controls the mass ratios. The phase varies there by o(1) because |ξ|ρ^j→0, and so does 1/ω_ξ(a(x)). Hence the normalized cancellation gives

    o(Z) = W_ξ(1)∑_{|w|=j} exp(iξx_w) + o(Z)
         = m^j W_ξ(1) [∫exp(iξx)dν(x)+o(1)].

The last equality uses ν(K_w)=m^(−j) and the same phase variation. The comparability m^j W_ξ(1)≍Z proves the assertion for all real frequencies.

**Required output regularity.** Sobolev regularity H^β locally implies (1), including at the endpoint. For a fixed smooth compact ψ choose η smooth, equal to one on supp ψ. The weighted Fourier transform of ψ times the output is bounded by the convolution of a Schwartz-decaying kernel (1+|s|)^β|ψ̂(s)| with the L² function (1+|s|)^β|widehat{η Law(Y)}(s)|, using the weight product inequality. A convolution of two L² functions vanishes at infinity, by approximation with compactly supported L² functions and Cauchy–Schwarz.

Mere L¹ local output does not imply Rajchman. Take K the middle-thirds set on [0,1], v=min(1,D^α) and X with a strictly positive smooth density. Since K is Lebesgue-null, v(X)>0 almost surely and Tonelli gives a global probability density in L¹ for Y. But for ν the independent digit product, |ν̂(2π·3^k)|=∏_{j≥1}|cos(2π/3^j)|>0; every factor is nonzero and the sum of its tail deviations from 1 is finite. Thus some quantitative regularity at the relevant mass scale is necessary for this implication. No claim that H^β is the weakest sufficient function space is made.

**An elementary family of noninteger contraction ratios.** For every integer p≥2 set B=(p+√(p²+4))/2 and ρ=1/B. Take the two similarities x↦ρx and x↦ρx+1−ρ on [0,1]. For ξ_k=2π B^k/(1−ρ) the independent-digit product gives

    |ν̂(ξ_k)| = ∏_{h=−∞}^k |cos(π B^h)| → [∏_{h≥1}cos(πρ^h)]² > 0.

For h≥1, B^h+(−ρ)^h is an integer by the linear recurrence s_h=p s_(h−1)+s_(h−2), s_0=2,s_1=p. Thus each positive-index modulus equals cos(πρ^h); the h=0 factor has modulus one. All the remaining factors are positive since ρ<1/2, and their deviations from one are summable. This explicit non-Rajchman family is excluded at any local germ by the theorem, for arbitrary allowable coefficient trace and exterior variance zeros.
