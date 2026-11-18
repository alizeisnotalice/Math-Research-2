# Exponentially accurate cancellation across two distinct continuous wells

For every prescribed κ>0 there is a bounded even continuous graph with a nonempty canonical core of height T>0, and two distinct right remote wells (reflected on the left), for which

    |E exp(tX−t²f(X)/2)−1|=O(exp[−Tt²/2+(b−κ)|t|]).

It can have Ef=1 and a nonzero fourth transform residual. Nevertheless its combined raw reciprocal-center density has an oscillation of nonvanishing relative amplitude on the √t scale. Round 144 obtained only leading-order cancellation by adding two cosines; here a positive first-passage kernel cancels their entire output in the constant-distance model.

Write K_vρ(t)=∫₀^∞p_{vm}(t−m)ρ(m)dm, where v is the reciprocal of the distance from the core edge. Transfer the canonical surviving measure to the coordinate u=b−y, writing η(du) there; near zero its density is Ju+o(u). The required edge output is A(t)=∫exp(−tu)η(du). Set

    ρ_v^*(m)=∫(1+vu)exp[−m(u+vu²/2)]η(du).

Then m²ρ_v^*(m)→J and K_vρ_v^*=A for t>0. The limit follows by rescaling u with m and using the linear edge bound. The identity follows by completing the square in the Gaussian integral: K_v[exp(−λm)](t)=(1+2vλ)^(−1/2) exp[−t(√(1+2vλ)−1)/v]. Here one can use ∫₀∞s^(−1/2)exp(−a/s−cs)ds=√(π/c)exp(−2√(ac)), a,c>0.

For v₁>v₂>0 define the substochastic first-passage kernel Π by the density

    (v₁−v₂)m exp[−(n−m)²/(2(v₁m−v₂n))]/√(2π(v₁m−v₂n)³),    0<n<v₁m/v₂.

Its exact identity on t>0 is K_{v₁}ρ=K_{v₂}(Π_#ρ). Indeed start Brownian motion at m and before time v₁m stop when it first meets the line n=(v₁m−s)/v₂. Every positive endpoint at time v₁m requires a prior hit. The hitting-time density is that of first passage to a positive constant by Brownian motion with drift 1/v₂ (reflection with exponential tilting); substituting s=v₁m−v₂n gives the displayed formula. At the hit the remaining transition is Gaussian with variance v₂n by the strong Markov property. This proves the density identity, first almost everywhere and then everywhere at t>0 by continuity, for the densities used below.

Fix k>0, h(m)=1_{m≥1}m⁻²cos(2k√m), and choose α_i>0 with sum 1. Take the two center densities

    α₁ρ^*_{v₁}+δh,         α₂ρ^*_{v₂}−δΠ_#h.

Their frozen output sums exactly to A. At large n,

    n²(Π_#h)(n)=exp[−(v₁−v₂)k²/2]cos(2k√n)+o(1).

For verification set m=n+√n w in the displayed density. On bounded w its rescaling tends to the centered Gaussian of variance v₁−v₂, and 2k√m=2k√n+kw+o(1). Where |m/n−1|≤c for small fixed c the rescaling is bounded by C exp(−c' w²); at m≥2n the integrand is exponentially small in m. On the remaining comparable region away from m/n=1 the exponent gives exponential smallness in n, including near the lower support endpoint by absorbing the prefactor (v₁m−v₂n)^(−3/2) in half the exponential. This proves the limit and boundedness of n²Π_#h at infinity. Both densities are therefore positive for m≥M with M large and δ>0 small. Their sum, times m², has the nonzero oscillatory term δ(1−exp[−(v₁−v₂)k²/2])cos(2k√m), in addition to J+o(1).

Truncate both to m≥M and realize them on actual continuous flanks at a_i=b+1/v_i by the terminal integral equation of round 144. The omitted bounded-center kernels cost O(exp(−c t²)); the densities before truncation are integrable on (0,M) since Π is substochastic. Write D=min_i 1/v_i. The actual distances differ from their frozen values by O(m^(−5/2)exp(−Dm/2)). On t/2≤m≤2t the mean value theorem on the Gaussian density gives an integrated error at most a polynomial times exp(−Dt/4). For m<t/2 and m>2t both frozen and actual kernels are bounded by a polynomial times exp(−Dt/4), since (t−m)²/(2m)≥t/4 and actual distances only increase. These bounds hold for large t with fixed M. Taking D>8(κ+1) makes the overall error O(exp(−κt)).

At the core edge and each unused flank use f−T=1/√log(1/u) for sufficiently small positive distance u. As in round 144 each costs at most exp(−c t^(4/3)+C t). The intervening compact pieces have a positive height gap, and the constant far tail is taken above 1+T. The mirrored distant wells on the opposite wing cost O(exp(−κt)) after edge normalization by their spatial distance; the mirrored core edge is again negligible by its high cusp. Complete the graph with the even compact mean adjustment and a mean-zero fourth-moment perturbation of round 144, choosing a sufficiently small initial core. All compact adjustments take place with a positive height gap. Multiplication by exp(−Tt²/2+b|t|) gives the claimed transform error. The estimate concerns both remote directions; the constructed graphs need not match the transform at finite tilts.
