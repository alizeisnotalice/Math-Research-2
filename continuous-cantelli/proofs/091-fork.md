# Uniqueness and symmetry of the early true exit, with endpoint bounds

Assume the established true-exit conclusion of rounds 83–84 on 0<t<T. Thus Brownian motion starts at zero, survives inside expanding bounded intervals I(t)=[−A(t),B(t)], and its paths stopped by time t have spatial law g(x)1_{I(t)}dx, with g the standard normal density. The endpoints are continuous on the early window under those results. Write q(t,x)dx for its surviving law. Then 0≤q(t,x)≤p_t(x), and q=0 almost everywhere outside I(t). All measures here use Lebesgue normalization.

## Uniqueness of the early barrier

More generally, two such early true-exit solutions with the same initial point mass and the same prescribed spatial density g≥0 are identical up to each common time. Assume g is locally bounded and strictly positive on the swept interval for uniqueness of the endpoints. Set

    H_i=q_i+g1_{I_i},     q_i=(H_i−g)_+.

The second equality follows directly from support and positivity, even if q_i happens to vanish inside I_i. The heat evolution satisfies ∂_t H_i=½∂_xx q_i with common initial measure δ₀. Indeed H_i is the law of B_(t∧τ_i), and Itô's formula for a smooth compact test gives H_i(t)−δ₀=½∂_xx∫₀ᵗq_i(s)ds in distributions. This is the monotone enthalpy equation for the expanding interval; no uniqueness theorem for a globally specified stopping time is required.

Let Q=q₁−q₂, K=H₁−H₂, W(t)=∫₀ᵗQ(s)ds. On every [0,T'] with T'<T, Q,K are in L²(dt dx): use ∫p_t²dx=O(t^(−1/2)), bounded union of the intervals and the bounded g there. We have W''=2K and W_t=Q. The map z↦(z−g(x))_+ is monotone and 1-Lipschitz, so KQ≥Q² pointwise. Spatial Fourier truncation to |ξ|≤N, temporal integration and W(0)=0 give

    (1/2)∥∂_x W_N(t)∥₂² = −2∫₀ᵗ∫ K_N Q_N dx ds.

Here the projection is orthogonal, all derivatives of the truncated functions are in L², and the right side tends to −2∫KQ by L² convergence. Thus ∫KQ≤0, forcing Q=0 a.e. and K=0 via the distribution identity. The weakly continuous stopped laws H_i are equal at every time; alternatively choose a full-measure sequence of matching times tending to each given one and use continuity of the endpoints. Hence the two intervals coincide pointwise in time, since g>0 distinguishes their interval indicators. No Brownian embedding continuation past T is part of this argument.

For even g and zero starting point the reflecting solution is another candidate. Therefore A(t)=B(t)=b(t), and q and H are symmetric for the full early interval. For the original graph the two inverse boundary functions agree at every swept radius, and f(x)=f(−x) on the swept interval. If the endpoints tend to a finite radius b as t↑T, the limiting surviving law at T is also symmetric and supported on [−b,b].

## Conditional survivor concentration

At each early time the ratio q(t,y)/p_t(y) has an even nonincreasing representative in |y|. To prove this without boundary regularity estimates, first test membership in closed symmetric intervals at finitely many positive times ≤t. Conditionally on B_u=y, the earlier B_s has law N((s/u)y,s(1−s/u)), s<u. Frozen-time multiplication by an interval indicator preserves even nonincrease. Convolution with an even nonincreasing Gaussian preserves the same property, as is seen from layer-cake reduction to symmetric intervals. Induction therefore proves the claim for the bridge probability at the finite time grid. On any [ε,t] impose margins b(s)−1/k, first check a countable dense grid including endpoints and then use continuity. Increasing k recovers strict survival on [ε,t] by compactness. Finally let ε↓0. All limiting probabilities preserve even nonincrease and give survival to t; hitting exactly at the deterministic time t is a null event. The same limits pass to T along the early interval.

It follows that conditionally on survival, |B_t| is stochastically dominated by |N(0,t)| conditional on being ≤b(t): the Radon–Nikodym derivative relative to the latter is nonincreasing in the radius, so its covariance with each increasing test is nonpositive. In particular with a=b(t)/sqrt(t) its conditional second moment is at most

    t [1−2a g(a)/(2Φ(a)−1)].

The future moment identity from the exact graph reads ∫y² q(t,y)dy=∫_{|x|>b(t)}[x²−(f(x)−t)]γ(dx). Therefore, now writing b=b(t),

    E_γ[f(X) | |X|>b] ≥ 1 + b g(b)/barΦ(b) + t·2a g(a)/(2Φ(a)−1).   (1)

Higher increasing radial moments of the surviving measure have the analogous truncated-Gaussian bounds. Non-strict inequalities suffice and require no smoothness of the boundary.

## A logarithmic correction to the lower remote wing

Suppose b=b(T−)<∞ as in the original reduction and T>0. At T the survivor satisfies the stronger endpoint density estimate

    q_T(y) ≤ p_T(y)−p_T(2b−y)
            = p_T(y)[1−exp(−2b(b−y)/T)]       (−b<y<b).

Indeed its path has always stayed below b, and the right side is the half-line killed heat kernel by reflection. This holds also by weak passage from t<T; deterministic time T cannot create a path hitting b with positive probability. Thus q_T(y)≤C(b−y), and

    ∫ exp(zy)q_T(y)dy ≤ C exp(bz)/z²     (z>0).                    (2)

The constants are finite for the fixed b,T. Reflection gives the left analogue. The full graph identity gives the matching future transform ∫_{f≥T} exp(zx−z²(f(x)−T)/2)dγ, with the prior exact boundary convention at T if a plateau occurs. Its positive integrand on every far-tail interval where f>T gives a lower bound independent of that convention. If the positive wing is eventually nondecreasing, choose x∈[z−1,z] beyond all such thresholds; using the shifted Gaussian under x−z gives

    c exp[(1+T−f(z))z²/2] ≤ C exp(bz)/z².

In particular, if its tail limit is H_*=1+T, then

    H_*−f(z) ≤ 2b/z−4log(z)/z²+O(1/z²).                          (3)

A corresponding reflected estimate holds for the negative wing if it is eventually nondecreasing in the outward radius. Without monotonicity (2) still bounds the exact exponentially weighted Gaussian window. These consequences use the genuine surviving-law region provided by earlier rounds and do not establish true stopping beyond that region.

## The first asymmetric level and its quantitative flatness threshold

In the single-well setting both outward functions F_+(x)=f(x) and F_−(x)=f(−x) are nondecreasing, x>0, with limits H_±>1. Set T=min(H_+,H_−)−1. The early result gives the same radius b>0 on both wings as t↑T, and F_+=F_− on [0,b]. The future transform of the limiting *preactivation* law Q_(T−) counts all f≥T, including level-T sources outside the already stopped interval:

    ∫e^{zy} q_T(y)dy = ∫_{f≥T} exp(zx−z²(f(x)−T)/2) dγ(x).

This follows by subtracting the true stopped part f<T from the full identity and taking the weak limit from t<T. Contributions at x=±b have zero γ mass. In particular for E_h={x≥b:T≤f(x)≤T+h}, inequality (2) yields γ(E_h)≤C exp(z²h/2)/z². Taking z²=2/h proves γ(E_h)=O(h), and reflection gives the same bound on the left. Thus there is no plateau at T and

    liminf_{u↓0} [F_±(b+u)−T]/u >0.                               (4)

Indeed the interval [b,b+u] lies in the appropriate E_h for h=F_±(b+u)−T, with Gaussian mass asymptotic to g(b)u. More explicitly limsup_{z→∞} z²e^(−bz)∫e^{zy}q_T(y)dy≤2b p_T(b)/T by the static-boundary kernel, so (4) is at least T g(b)/(e b p_T(b)).

Suppose now H_+>H_−=1+T, the reverse case following by reflection. Let c₀=liminf_{u↓0}[F_+(b+u)−T]/u∈(0,+∞] and Δ(u)=(F_−(b+u)−F_+(b+u))_+. Then

    liminf_{u↓0} u log(1/Δ(u)) ≤ b²/(2c₀),                       (5)

using +∞ when Δ=0 and 1/+∞=0. In particular, points with F_+<F_− occur arbitrarily near the right of b, the opposite inequality from the remote ordering. The equality germs cannot extend any positive distance past b.

To prove (5), symmetry of q_T makes the integral over b<x<∞ of

    sinh(zx) [exp(−z²(F_+(x)−T)/2)−exp(−z²(F_−(x)−T)/2)] g(x)

zero for all z>0. At large fixed R>b one has F_+≥1+T+δ and F_−≤1+T on x≥R for some δ>0. The integrand is negative there, and its integral on [z,z+1] is at most a fixed strictly negative constant for all large z: use exp(−z²(F_−−T)/2)≥exp(−z²/2), the exponentially small relative contribution of the F_+ term, and exp(zx−z²/2)g(x)=g(x−z). On each fixed [b+η,R] the positive part is O(exp(zR−c_η z²)), since neither wing has a plateau at T.

If (5) failed, choose k>0, 0<c<c₀ with 2ck>b² and, on sufficiently small u∈(0,η), Δ(u)≤exp(−k/u), F_+(b+u)−T≥cu. Using 1−exp(−v)≤v, the total positive contribution there is bounded by

    C z² exp(bz) ∫₀^η exp[−(cz²/2−z)u−k/u]du
       ≤ Cη z² exp[bz−2sqrt(k(cz²/2−z))] →0.

This contradicts the fixed negative mass, proving (5).

Consequently a nonconstant exact continuous single-well graph cannot have both reflected wings in a Gevrey class of order s<2 on a neighborhood of their first possible splitting radius b; it suffices that their difference belong to this class there. The equal-tail case was already excluded by the earlier connected-exit budget. In the unequal-tail case the difference and all its derivatives vanish on the left of b. If 1<s<2 and the local derivatives obey C A^n(n!)^s, Taylor's remainder from b gives C(Au)^n(n!)^(s−1). Taking n the integer part of (2Au)^(−1/(s−1)) gives a bound C' exp(−c' u^(−1/(s−1))), contradicting (5). For s=1 the same Taylor bound tends to zero for Au<1. This argument does not exclude order two at all possible constants.

At the order-two endpoint, if the flat wing difference h has a bound |h^(n)|≤C A^n(n!)² on some neighborhood of b, choosing n=floor(1/(Au)) in Taylor's remainder gives liminf_{u↓0} u log(1/Δ(u))≥1/A, by Stirling's estimate. Accordingly the local derivative type must obey A≥2c₀/b². In particular the Beurling (zero-type) Gevrey-2 germ cannot occur. If c₀=+∞, every finite order-two type is excluded.

There is a joint time-scale tradeoff sharpening (3) and (5). Write h(u)=F_+(b+u)−T and set

    κ=liminf_{u↓0} h(u) log(1/Δ(u)),   θ=limsup_{z→∞} z[H_−−F_−(z)].

Then sqrt(2κ)+θ/2≤b. Indeed θ∈[0,2b] by monotonicity and (3). For any 0<k<κ choose η small enough that Δ≤exp(−k/h) on (0,η). The integrated positive difference on this interval is at most Cη z² exp[(b+η−sqrt(2k))z], since z²h/2+k/h≥sqrt(2k)z. The entire middle positive part still decays quadratically in z in its logarithm. On the far interval [z−1,z], use monotonicity of F_− and the translated Gaussian mass to lower-bound the negative contribution in absolute value by c exp[z²(H_−−F_−(z))/2]. Comparison on a sequence approaching θ, followed by η↓0 and k↑κ, proves the stated tradeoff; when κ=0 it reduces to (3)'s leading estimate. In particular θ=2b forces κ=0, excluding every fixed exponential of the inverse local time as an upper bound for the required reversed fork.
