# Detection of arbitrarily power-small oscillations on the intrinsic noise width

Use the exact identity D_e=p_{M+e} and g=p₁ as in r063-intrinsic-periodic.md. At a point x₀ suppose on the positive side u↓0
v(x₀+u)=u^q [c+η(u) H(ωu^{1−q/2})+o(η(u))],
where q>2, c>0, ω≠0, H is bounded continuous periodic, η>0 tends to zero. Write h=r^{q/2}, d=h/r, and assume η(r+ht)/η(r)→1 for each fixed t and sup_{r/2≤u≤3r/2}η(u)/η(r)≤C. If for some fixed positive integer k we have d^k=o(η(r)), then H is constant. A subsequence r→0 with this comparison and the same local ratio assumptions is sufficient. In particular every η(r)=r^β(log(1/r))^ν tending to zero, β≥0, ν∈ℝ, qualifies by choosing k large. There is no assumption on v on the opposite side of x₀ and no differentiability condition on H or on the o(η) remainder.

Let s be fixed. The scaled kth derivative is
h^k D_{h²}^{(k)}(x₀+r+hs)=∫ g(x₀+r+ht) p^{(k)}_{1+v(x₀+r+ht)/h²}(s−t)dt.
The exact target is O_k(h^k) since M>0. For T>0 and y≠0 one has |p_T^{(k)}(y)|≤C_k |y|^{−k−1} by maximizing the scaled Hermite polynomial times its Gaussian. Thus input outside |t|≤1/(2d) contributes O_k(d^k). On the remaining interval write B(t)=c(1+dt)^q for the baseline rescaled variance and Δ=v/h²−B. Uniformly there c₂≤1+B≤C₂ and |Δ|≤Cη(r). The first-order Taylor formula in the variance and ∂_T p^{(k)}_T=p^{(k+2)}_T/2 give an integrated O(η(r)^2) remainder, because all intermediate variances stay in a fixed positive compact interval and the Gaussian derivative envelopes are integrable in t. The linear term is ½∫g Δ p_{1+B}^{(k+2)}(s−t)dt.

The baseline term ∫_{|t|≤1/(2d)} g(x₀+r+ht)p_{1+c(1+dt)^q}^{(k)}(s−t)dt is O_k(d^k). To verify cancellation, introduce the real parameter b=dt and Taylor-expand g(x₀+r+rb)p_{1+c(1+b)^q}^{(k)}(s−t) at b=0 to degree k−1. On |b|≤1/2 its kth b derivative is bounded by C(1+|s−t|^L)e^{−C'|s−t|²}, uniformly for small r, so the integrated remainder is O(d^k). For each j<k the Taylor polynomial terms are d^j t^j times fixed linear combinations (bounded coefficients depending on r) of p_{1+c}^{(k+2l)}(s−t) with 0≤l≤j. After extending their exponentially small finite interval tails to the line, every integral vanishes by j<k and integration by parts. This justifies subtracting even a nonconstant power-law background at this precision without estimates on derivatives of v.

Divide by η(r) and take a subsequence for which the phase ω r^{1−q/2} converges on the period circle to θ. Since (1+dt)^q→1 and the derivative of the phase at r times h is ω(1−q/2), one has
Δ/η(r)→H(θ+ω(1−q/2)t)
on compact t sets pointwise. Dominated convergence in the linear term gives, for every s,
∫ p_{1+c}^{(k+2)}(s−t) H(θ+ω(1−q/2)t)dt=0.
The rescaled periodic function on the left has Fourier multipliers (iκ_n)^{k+2}exp[−(1+c)κ_n²/2], nonzero for every nonzero harmonic. Hence all nonconstant Fourier coefficients vanish. This does not constrain a constant H in this test. In fact the same reasoning applies to bounded measurable periodic H if the value asymptotic is meaningful: on compact t sets the phase maps converge C¹ to an affine map of nonzero slope, so approximation of H by continuous periodic functions in L¹ gives convergence in measure, with uniform Jacobian bounds after change of variables.
