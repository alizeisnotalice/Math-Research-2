# Quadratic growth of the local critical coupling on both wings

Let h=f−1, B=sup h, D=−inf h, both positive for a nonconstant bounded continuous exact normalized graph. Fix w>1 and v=(w−1)/w. At each center x there is a unique b=b_w(x)>0 for which P_w exp(−bh)(x)=1. This follows from the strict negativity for small positive b in round 101, positivity for large b from a fixed interval where h<0, and monotonicity of the divided difference J. Then for either x→+∞ or x→−∞,

    sup_{0<δ<1} min{δ²/[2w²(1+vB)], ((1+(1−δ)/(w−1))²−1)/(2wD)}
        ≤ liminf b_w(x)/x² ≤ limsup b_w(x)/x² ≤ 1/(2wD).

The lower bound is strictly positive. No shape, regularity at infinity or tail limits are imposed. The upper bound follows by testing a fixed interval where h<−D+ε; at the root its positive integral against the shifted Gaussian is at most one. Take logarithms, divide by x², then send ε to zero.

For the lower bound reflect to x>0, fix δ∈(0,1) and choose a positive c strictly below both displayed lower-bound quantities. Take b=c x². In the bridge identity of round 101, given y let m=(x+(w−1)y)/w, T∼N(m,v) and H=h(y). The integrand against p_w(x−y)dy is
K_b(H)−E K_{T²/2}(H), K_a(H)=(exp(−aH)−1)/a.
On y≥−(1−δ)x/(w−1), we have m≥δ x/w. For 0<u≤1,

    E exp(−uH T²/2)=(1+vuH)^(−1/2) exp(−uHm²/[2(1+vuH)]).

Since b<m²/[2(1+vB)], this expression is at most exp(−buH) if H>0 and at least exp(−buH) if H<0. Integrate after multiplying by −H to obtain a nonpositive bridge integrand throughout that region.

There is a strict negative reserve on any fixed bounded interval where H≥η>0. Here m=x/w+O(1) and E K_{T²/2}(H)=−(2+o(1))/m² uniformly. Indeed on |T−m|≤sqrt(m), T/m→1 uniformly and the exponential in K is negligible; on its complement |K|≤B and the normal tail is exponentially small in m. Meanwhile K_b(H)=−(1+o(1))/(c x²). Since c<1/(2w²), the reserve in absolute value after integration is at least c₀ x^(−2)exp(−x²/(2w)−C x) for all sufficiently large x.

On the remaining opposite tail the integrand is bounded above by B+(exp(bD)−1)/b, since K_a(H)≥−B for a≥0. The Gaussian weight of this tail is at most exp(−(1+(1−δ)/(w−1))²x²/(2w)). By the second strict constraint on c its entire possible positive contribution is exponentially smaller than the reserve. Therefore J(w,x,c x²)<0 for all sufficiently large x. This holds for every indicated δ and c, proving the bounds.

**Geometry at finite centers.** b_w(x) is analytic since J_b>0. Its global minimum is β(w) at the unique center in round 101. Every level above the minimum has precisely two simple spatial crossings. In fact J is negative in both tails and positive somewhere. If there were three distinct spatial zeros, the signed measure with density j_b(t) p_{w−1}(t) (j_b=J(1,·,b)) would annihilate three real exponentials with distinct exponents. A nonzero linear combination of them can be chosen to vanish at ±√(2b). Such a combination has at most two zeros counted with multiplicity: divide by its first exponential and apply Rolle to the derivative, which is a sum of two exponentials. Its signs therefore alternate exactly where j_b changes sign, contradicting the annihilation. A spatial zero with zero derivative already forces J negative at all other centers by round 101's affine-chord argument. This gives the simple crossings and strict monotonicity of b_w on each side of its minimum. At any finite extinction width w₀>1 with b=β(w₀), the two crossing centers for w↑w₀ are ξ(w₀)±√(w₀−w)+O(w₀−w): Taylor-expand at the contact with J_x=0, J_xx<0 and J_w=J_xx/2, then use analytic implicit functions in √(w₀−w) after dividing the expansion by w₀−w.
