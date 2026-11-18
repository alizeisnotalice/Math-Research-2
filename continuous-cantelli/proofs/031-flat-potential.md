# Flat contacts via the positive call potential

Let x₀ be an attained zero of the continuous bounded variance field v≥0. Write ρ=p₁. Exactness says X+√v(X)Z has law N(0,M) for independent standard normals X,Z. The following use only this law and bounds on v in a shrinking neighborhood. If M=1 then Ev=0 from the variance identity, so the nontrivial case has M>1.

## Positivity, analyticity and the exterior

For s≥0 define H_s(y)=E[(y+√s Z)_+]−(y)_+=½∫₀ˢ p_t(y)dt. The formula follows by the heat equation in distributions or by differentiating the explicit Gaussian call formula for s>0 and taking the limit at s=0; H₀=0. It is even, nonnegative, nondecreasing in s, with ∫H_s=s/2. Exactness gives
L(y)=∫ρ(x) H_{v(x)}(y−x)dx = ½∫₁ᴹ p_t(y)dt
for real y. Both signs of the call parameter yield this formula since the centered laws are Gaussian. L has an entire continuation. For small h>0 let S_h(y) be the same positive integral restricted to |x−x₀|<8h. It agrees on the real line near x₀ with L minus the exterior potential. The latter is holomorphic for |z−x₀|<3.1h, uniformly bounded there as h↓0: for each exterior mean x, Re[(z−x)²]≥[(8−3.1)²−3.1²]h²>0. Thus absolute complex kernels ≤(2πt)^−1/2; integrate in t up to v(x)≤M, and then against ρ. Zero-variance means contribute identically zero. Uniform domination and Morera give holomorphy. Consequently S_h(x₀+h·) on the disk of radius 3 has an analytic continuation bounded by a constant B≥1 independent of h. This construction does not analytically continue the integrand across its own mean.

If U(r)=sup_{|u|≤r} v(x₀+u)=o(r²), then for I=[5/8,7/8] and J=[1/2,1], all small h satisfy
∫_{x₀+hJ}S_h(y)dy ≥ c₀ ∫_{hI} v(x₀+u)du.
Indeed for u∈hI, all t≤v(x₀+u) are o(h²) uniformly, so a N(x₀+u,t) variable falls in x₀+hJ with probability ≥1/2; also ρ is bounded below. The same result holds with both intervals reflected. This lower bound is solely positivity and Tonelli. An upper bound on S_h on J thus gives a bound on the average of v over hI.

## Elementary propagation from short real segments

Let F be holomorphic near |z|≤3 with |F|≤B there. If |F|≤e≤B on [−α,α], 0<α≤1/4, then for |y|≤1
|F(y)|≤4B(e/B)^θ,  θ=log 2/log(12/α).
For N≥2 interpolate at N equally spaced points on [−α,α]. The sum of absolute Lagrange coefficients for |y|≤1 is at most (euler·(1+α)/α)^{N−1}≤(6/α)^N, using (N−1)!≥[(N−1)/euler]^{N−1}. Cauchy's contour formula for the interpolation remainder on |z|=3 bounds its modulus by (3B/2)[(1+α)/(3−α)]^N≤(3B/2)2^−N. Use N=floor(log(B/e)/log(12/α)); for N<2 the same stated bound follows from |F|≤B. Let e↓0 for the zero case. With initial segment [−1/4,−1/8], its length 1/8 gives coefficient bound [2·euler·(5/4)/(1/8)]^{N−1}≤64^N and the same remainder ≤(3B/2)2^−N, hence the same conclusion with θ₀=1/7.

## One-sided plateau and arbitrary opposite flatness

Suppose v(x₀+u)=0 on −r<u≤0 and U(r)=o(r²) as r↓0. For y=x₀+hs with −1/4≤s≤−1/8, all contributing local input means are at least h/8 away. Monotonicity of H and the Gaussian tail bound give
S_h(y)≤‖ρ‖∞ U(8h) exp[−h²/(128U(8h))].
Here ∫_{|w|≥b} H_U(w)dw=(1/2)∫₀ᵁP(|N(0,t)|≥b)dt≤U exp[−b²/(2U)]. Propagation and the local mass bound give constants C,c>0 for which
h^−1∫_{5h/8}^{7h/8} v(x₀+u)du ≤ C exp[−c h²/U(8h)].
This holds without a lower bound, regular variation or derivative assumption on the positive side. A direct corollary excludes positive-sided comparability to exp(−d/u^γ) or distinct exponents d within fixed multiplicative bounds on the action.

## Two-sided flat sources

For |y−x₀|≤αh split local means into |u|≤2αh and |u|>2αh. By monotonicity and ∫H_U=U/2 the first contribution is ≤‖ρ‖∞U(2αh)/2. The second is at distance ≥αh, so
e_h := sup_{|y−x₀|≤αh} S_h(y)
 ≤ ‖ρ‖∞[U(2αh)/2+U(8h)exp(−α²h²/(2U(8h)))].
Consequently the two-sided necessary condition is
h^−1∫_{5h/8}^{7h/8} v(x₀±u)du ≤ C (e_h/B)^{log 2/log(12/α)}
for small h. The constants B,C are independent of α∈(0,1/4].

In particular assume for all sufficiently small |u| that v(x₀+u)≤C₁ exp(−d/|u|^γ) with d,γ>0. If on either side v(x₀±u)≥c₁ exp(−D/u^γ) for some finite D, there is no exact graph. For any fixed α, e_h≤C exp[−d(2α)^−γ h^−γ] for small h, since its remaining term is double exponentially small on this scale. Choose α so d(2α)^−γ log 2/log(12/α)>D(8/5)^γ. The necessary integrated upper bound then contradicts the integrated lower bound from hI. In fact the upper-envelope hypothesis alone forces, on each side, the annular integral divided by h to be O_K(exp(−K/h^γ)) for every K>0. Pointwise lower bounds can be replaced throughout by a corresponding lower bound for that annular integral along any sequence h↓0.

This excludes canonical bilateral infinitely flat contacts, and yields a stronger necessary condition at plateau boundaries. It does not exclude arbitrarily thin spikes or infinitely flat profiles lying below every specified action.
