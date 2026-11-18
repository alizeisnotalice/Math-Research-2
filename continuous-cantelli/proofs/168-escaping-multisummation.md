# Arbitrary escaping Laurent–Puiseux branches using same-direction multisummation

Use the setting, real Gaussian density g, existing early cut y, and integral G(h) of [round 167](167-temporal-quasianalyticity.md). Finitely many escaping poles of any positive orders can be combined with any of the finite endpoints above. For a pole of order p>0 put k=1+2p and choose

    π/(2k)<θ<θ'<min(π, π/(4p)).

Such angles exist since 2p<k. Write R=|h|+u. On the Cauchy circle of radius c₀R about h the argument of h+u remains in the larger sector and its modulus stays comparable to R. After shrinking δ the integrand is bounded by

    C u⁻¹ᐟ² R^(−M) exp(−d/(u R^{2p}))
       ≤ C u⁻¹ᐟ² R^(−M) exp(−d/R^k)

for constants M,d>0. Indeed the real coefficients of the leading Puiseux pole and 2pθ'<π/2 give a strictly positive lower bound of order R^(−2p) on the real part of (y−F_j)². The density g(F_j) is bounded there as well; the remaining derivative is a fixed negative power at worst. Cauchy differentiation multiplies this estimate by n!(c₀R)^(−n). Maximizing R^(−M−n)exp(−d/R^k) over R>0, integrating u⁻¹ᐟ² and using Stirling gives

    |G^(n)(h)| ≤ C A^n n! Γ(1+n/k).

Limits of all derivatives at zero exist by the same estimates at each fixed order. Taylor's remainder along radial segments gives strong Gevrey asymptotics of order 1/k on a sector of opening greater than π/k; hence G is the k-sum of its Taylor jet in the positive direction. Finite endpoints have k=1 by the preceding argument. The remaining regular terms are convergent and can be included at any level. Group summands of equal k. By the same-direction multisummability characterization and independence of the sum from the decomposition, their analytic sum is uniquely determined by the sum of their jets. A zero jet thus forces a zero germ. This invocation of multisummability is needed here only to combine potentially disparate orders without imposing the common wide-sector estimate of the elementary proof. [Balser, sections 2–3, in particular equations (3.1)–(3.2), pp. 164–166](https://doi.org/10.7146/math.scand.a-12392)


## Countably many inverse branches

Within any one of finitely many orders k, the argument also applies to countably many branches provided the sector, its Cauchy-circle fraction and upper age bound are common in the group, and their integrand absolute values on all those circles have a summable majorant of the displayed form: replace C by C_j with Σ C_j finite and retain common M and d>0. For k=1 use the finite-endpoint majorant u^(-M−1/2) exp(−d/u), including the preceding factor C_j. The older part of the age splitting, whose endpoint in source time moves, also requires normal convergence of its holomorphic extensions on a common complex disk about h=0. These assumptions justify summing before the Cauchy and radial estimates; the entire group has the same certified summability order. No bound on the number of connected components of {f<1} is imposed.

For example a sufficient source description for the finite-endpoint group consists of analytic or sectorial branches with endpoint a_j→±∞, with |X_j(T−v)−a_j|≤D and |X'_j(T−v)|≤C|v|^(-M) throughout the common slightly larger sector and Σ exp(−a_j²/8)<∞. In a sufficiently small sector disk and away from the finitely many exceptional indices, Re(X_j²)≥a_j²/2 and Re[(y−X_j)²]≥1. Thus the Gaussian source weight supplies summable constants, while the short-age kernel provides the common exp(−d/u). For the moving-endpoint positive-age segment, parameterize the contour by affine interpolation of its endpoints as in round 167. Both age and v stay uniformly near positive real numbers bounded away from zero. The same separated endpoint estimates and Gaussian summation give normal convergence there. The finitely many remaining indices require the usual separation and sector conditions at the chosen cut. One only applies this criterion when the branch family exhausts the relevant source band up to null sets.

If every possible extinction level meets these hypotheses, the contradiction with round 166 proves constancy. This covers mixtures of distinct positive escaping powers through a common temporal test even when the behavior of f below one admits no branch expansion.
