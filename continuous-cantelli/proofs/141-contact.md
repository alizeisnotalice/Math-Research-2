# Joining the logarithmic contact estimate all the way to the origin

Assume f≥0 bounded continuous and positive on (0,3b], for some b>0. On the favorable side suppose, with ℓ=log(1/x),

    x²/(2f(x)) ≤ ℓ+β log ℓ+C+o(1)   as x↓0.

Use the signed resolvent formula of R132 with X standard normal and mark f(X). For each fixed n≥0 the free derivative dominates its subtracted marked derivative uniformly for 0≤y≤b and q large if β<1/2, and also if β=1/2 and C<log(2)/2. The same half-neighborhood therefore has a smaller fixed interval containing zero where the first derivative has its strict sign for all q>0 if f is exact and normalized. Only the indicated side needs the contact condition. The comparison holds in the strong sense that the marked/free ratio tends uniformly to zero for β<1/2, and has limsup at most exp(C)/√2 at β=1/2.

Put z=√(2q), L=log z and α=zy; the free derivative is exp(−α)P_n(α)/z^{2n+1}, in the notation of R132. For α≤L/2 the global marked bound C_n/z in that normalization gives uniformly vanishing relative error. From now on α≥L/2. The integrand ratio of the marked to the free derivative at input x is

    exp(α−|α−zx|−qf) [Σ_{j=0}^n binom(n,j)(z²f)^j P_{n−j}(|α−zx|)]/P_n(α).    (1)

For x≤0 this is bounded by C_n(1+z|x|/α)^n exp(−z|x|−(1−O_n(1/α))qf); its integral is O(1/z). For x>3b the geometric exponent is ≤−zx/3 and for x on any compact subinterval of (0,3b] the factor exp(−qf) beats the O(z) favorable tilt. Polynomial factors are harmless there by bounded f and Gaussian tails. The estimates are uniform in α.

Consider now u=zx positive with x in a small neighborhood of zero. The geometric exponent is u−2(u−α)_+, at most u. A crude upper bound for (1) is

    C_n(1+u/α)^n exp[u−2(u−α)_+−qf(1−O_n(1/α))].                (2)

This follows from the positive coefficients of P_j and (1+z²f/α)^n≤exp(2n qf/α). For 0≤u≤1 the integral is O(1/z). For 1≤u≤L the denominator ℓ+βlog ℓ+C+o(1) is ≤L+O(log L), so u−(1−O_n(1/L))qf≤(3/4+o(1))L uniformly. For u≥L the denominator is at most L+(β−1)log L+C+o(1) by monotonicity in ℓ near infinity, or the coarser uniform O(1) constant if the fixed neighborhood has not yet been shrunk. Complete the square in the exponent. Since α≥L/2, the (1−O_n(1/α)) factor changes the effective denominator by O_n(1) and the center is 2L+O_n(log L). Thus, outside the window |u−2L|≤√L log L, (2) has integral bounded by a polynomial in L times exp(−c(log L)²). The extra factor (1+u/α)^n is absorbed into the Gaussian tail moment, and the integrand is with respect to g(u/z)du/z. All these contributions vanish faster than any negative power of L, except for already exponentially small terms.

Within the window the denominator is uniformly bounded by

    D_L=L+(β−1)log L+C−log 2+o(1),

i.e. qf≥u²/(4D_L)=:B₀. Here the o(1) may be chosen as a single uniform positive error. Since B₀∼L, the functions B^j exp(−B) for 0≤j≤n are decreasing for B≥B₀. Consequently every f-dependent factor in (1) can be replaced in an upper estimate by its value at qf=B₀, including its exact polynomial multiplier. If α<2L−2√L log L, the window has u−α≥√L log L, and the reflection loss exp(−2(u−α)_+) makes the entire window negligible. Otherwise α≥2L−2√L log L, and

    (|α−u|+2B₀)/α = 1+o(1)

uniformly (when α grows much faster than L this is immediate; in general use 2B₀=u+O(√L log L) and (u−α)_+=O(√L log L)). Monicity of the P_j and α≥L/2 then show that the polynomial quotient in (1), at B₀, is at most 1+o(1): all lower powers contribute O_n(1/L) since B₀=O(L) and |α−u|=O(α). Keeping the favorable geometric exponent at most u gives the final Gaussian-integral upper bound

    (1+o(1)) g(0) z⁻¹ ∫ exp(u−u²/(4D_L))du
     = (exp(C)/√2+o(1)) L^{β−1/2}.

The integral on this final line is over the whole real line to give an upper bound. Combining with the discarded regions proves the stated uniform bound. For the fixed-interval all-q first-order consequence, combine large q with the central strict estimate and its spatial Lipschitz bound on the compact remaining q range from 141-center.md.
