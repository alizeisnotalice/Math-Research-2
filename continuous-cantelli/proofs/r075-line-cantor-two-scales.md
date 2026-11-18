# Two-scale exclusion at every equal-ratio strongly separated Cantor germ

Use the real-line model and notation of r075-line-rajchman.md, including reflections, arbitrary extra variance zeros outside the local germ, and the uniform positive continuous coefficient trace in v=D^α(a+o(1)), 0<α<2. Write d=log(m)/log(1/ρ)<1, β=2(1−d)/α. A locally smooth output at even one point x₀∈K is impossible, for every such self-similar set, irrespective of the Fourier decay of its natural measure ν. It suffices that near this point the output have density F satisfying

    F(y)=P(y−x₀)+o(|y−x₀|^s),     s=(1−d)(2/α−1)=d+β−1>0,        (1)

for a polynomial P. A bound O(|y−x₀|^s') with s'>s works as well. The argument derives |ν̂(t)|≤C/|t| for all |t|≥2, which would make the singular probability ν absolutely continuous.

**Local geometry.** Work inside a fixed neighborhood U of x₀ where the density g is continuous and strictly positive, v may be compared above and below with D^α, and both traces exist. Choose the depth-j cylinder K_j containing x₀, write K_j=b_j+e_j R K, with R=ρ^j and e_j=±1, and let j→∞. Then |b_j−x₀|=O(R). If Δ>0 is the first-level set separation, the distance to every other depth-j cylinder is ≥(Δ/ρ)R. Take φ∈C_c^∞(R) equal to one in a fixed tube about K and supported in the tube of radius <Δ/(3ρ). In support of φ((x−b_j)/(e_jR)), D(x)=R dist((x−b_j)/(e_jR),K). No disjointness of the cylinder hulls is required.

The measure ν is d-Ahlfors regular: at any scale l the depth-k cylinders with ρ^k≍l carry mass m^(−k)≍l^d, are contained in intervals of length O(l), and are separated by ≳l. This gives both ball bounds for centers in K, using a suitably smaller cylinder for the lower bound. Thus, for 0<u≲1, one can cover K∩B(b,r+u) with at most C((r+u)/u)^d intervals of length O(u), by taking a maximal u-separated set and the two ball bounds. Consequently

    |K_u ∩ B(b,r)| ≤ C(r+u)^d u^(1−d).

The same inequality for u>1 and bounded r follows by enlarging C (or use 2r). Positive layer-cake integration and dyadic spatial annuli now give, uniformly for 0<σ≤R≪1 and M>d,

    ∫_U 〈(x−b_j)/R〉^(−M) exp[−c(D(x)/σ)^α] g(x)dx
                  ≤ C_M R^d σ^(1−d).                              (2)

Indeed the ball integral for radius l≥R is ≤C(l+σ)^d σ^(1−d), and the geometric shell sum converges. Restrict U to a fixed bounded interval where g is bounded. Outside U the same spatial weight integrated against the finite full source is O_M(R^M), regardless of its variance.

**A band-limited test.** Fix t≥2 first, set N=t/R, σ=N^(−2/α), z=(x−b_j)/(e_jR). Choose χ∈C_c^∞((−1/2,1/2)) identically one on [−1/4,1/4], and put φ̂_t(u)=φ̂(u)χ(u/t). Test the output against exp(it z)φ_t(z) in its own y-coordinate. For each fixed t this modulated Schwartz function has all its polynomial moments equal to zero because its Fourier support avoids the origin. Formula (1), |b_j−x₀|=O(R), and dominated convergence give

    ∫ exp(it(y−b_j)/(e_jR)) φ_t((y−b_j)/(e_jR)) dLaw(Y)(y)
                       = o_t(R^(1+s)).                              (3)

Indeed the global integral against P vanishes, local Taylor errors integrate to o_t(R^(1+s)), and the Schwartz tails outside any fixed neighborhood cost O_{t,M}(R^M) for the finite output and the polynomial for arbitrary M. Smoothness gives the same statement by Taylor expansion.

Given the source x, set A=N²v(x)/2, and write out the conditional expectation of this test by Fourier inversion. The difference from exp(it z)φ(z)exp(−A) satisfies the pointwise bound

    |error| ≤ C_M t^(−1) 〈z〉^(−M) exp(−cA),                       (4)

with C_M independent of t≥2, j, x and A≥0. In fact its Fourier amplitude in z after division by exp(it z) is φ̂(u) times

    χ(u/t)[exp(−A(1+u/t)²)−exp(−A)] + [χ(u/t)−1]exp(−A).

On |u|≤t/2 the undifferentiated bracketed exponential difference is bounded by C (|u|/t)exp(−cA), by the mean value theorem. Each positive-order u derivative is bounded by C_l t^(−l) exp(−cA); polynomial factors in A are absorbed into the exponential with c<1/4. Outside |u|≤t/4 the Schwartz tails of all derivatives of φ̂ give arbitrary inverse powers of t for the last term. The resulting Fourier amplitude and its derivatives through any fixed order have L¹ norms O_M(t^(−1)e^(−cA)). Integration by parts proves (4), also at deterministic sources. The reflection sign drops out of the Gaussian multiplier.

Since α<2, σ/R→0 for fixed t as R→0. Combining (2)–(4) shows

    |∫ exp(it z)φ(z)exp(−N²v(x)/2)g(x)dx|
         ≤ (C/t) R^d σ^(1−d) + o_t(R^(d+β)).                     (5)

Here in (2) use local comparability of v and select M>d+β to dispose of exterior sources. The integral on the left has source support near K_j. The order in (3) equals d+β. Although R^d σ^(1−d)=t^(−β)R^(d+β), t is fixed in this limit. The constant C multiplying 1/t in (5) is genuinely independent of this fixed t.

**The model cell.** Let a₀=a(x₀), g₀=g(x₀). On the support of φ(z), g(x)=g₀+o(1) uniformly, and where D(x)>0 one has v(x)=D(x)^α(a₀+o(1)) uniformly as R→0. The associated exponential mass error is o_t(R^d σ^(1−d)), by comparison with the weaker Gaussian damping as in (2). On changing variables the frozen integral on the left of (5) equals

    g₀ R ∫ exp(it z)φ(z)exp[−(a₀/2)(dist(z,K)/h)^α]dz,
                 h=σ/R →0.

For completeness, the fixed-reference packet masses normalized by their total in a fixed bounded tube T about K converge weakly to ν. Their total is ≍h^(1−d) by the tube bounds. For any fixed cylinder depth k take disjoint congruent small tubes about all its cylinders, each contained in T. Their masses are exactly equal by similarity and nearest-set separation, and the complement in T is exponentially negligible. Shrinking the fixed tube radii as k increases (e.g. proportional to ρ^k) proves the weak convergence by uniform continuity of tests. Choose T containing supp φ and a small neighborhood of K. It follows that the frozen integral divided by g₀R times this total tends to ν̂(t), since φ is one on K. The denominator is bounded below by c R h^(1−d)=c R^d σ^(1−d) with c independent of t. Divide (5), let R→0 along the nested cylinders, and obtain |ν̂(t)|≤C'/t for every t≥2. There is no restriction on the orientation signs because both source and output use the cylinder's oriented coordinate. Negative frequencies follow by conjugation.

Finally ν̂∈L²(R) would imply that ν is represented by an L² Lebesgue density, by Plancherel and uniqueness in distributions. The compact support K has zero Lebesgue measure by d<1. This contradicts ν(K)=1.

## Unequal contractions, including arbitrary arithmetic relations

The full two-scale conclusion also holds when the finitely many similarity ratios ρ_i∈(0,1) are unequal. Set d by ∑ρ_i^d=1, assume d<1 and use the natural cylinder probabilities p_w=ρ_w^d, where ρ_w is the product of its contraction ratios. For a scale r, stop each symbolic branch when ρ_w first drops below r. These cylinders have ratios in [r min_i ρ_i,r]. They cover K, have diameter O(r), carry masses comparable to r^d, and are pairwise separated by ≥Δr: two distinct stopping words have common prefix of contraction >r. The Ahlfors and tube bounds used above follow. For the occupied cylinder containing x₀ choose its stopping partition at r=ρ_w (or use the first words with contraction ≤r). Its distance from every other stopping cylinder is ≥Δr, so the identical local cutoff construction applies with R=ρ_w. If using a strictly-less stopping rule, take r slightly larger than ρ_w but smaller than its parent's contraction.

It remains only to justify the fixed-model tube probability limit used in the last step. Let

    Z(h)=∫_R exp[−(dist(x,K)/h)^α]dx,  F(s)=exp((1−d)s) Z(exp(−s)),
    a_i=−log ρ_i,  p_i=ρ_i^d.

Positive separation gives, as h↓0,

    Z(h)=∑_i ρ_i Z(h/ρ_i)+O(exp(−c h^(−α))).

To verify this, partition into disjoint fixed-radius tubes of the first-generation pieces. Inside each tube the nearest set is that same piece, and dropping or adding its exterior changes the integral only by the stated exponentially small amount. On the bounded complement every distance has a positive lower bound. On the two unbounded outer components integrate the one-dimensional tail exp(−(u/h)^α) from a fixed positive endpoint, absorbing any polynomial factor by decreasing c. Affine scaling of each whole-piece integral gives the formula. By the tube bounds, F is bounded above and bounded away from zero for large s. Its derivative is bounded there as well: differentiating the positive kernel introduces a multiple of (D/h)^α, absorbed in a weaker Gaussian-power kernel of the same order Z(h). Therefore translations of F along sequences escaping to +∞ have locally uniformly convergent subsequences, and every limit H is bounded continuous on R and solves

    H(s)=∑_i p_i H(s−a_i).

For completeness this identity forces each a_i to be a period of H, without any rational-independence or lattice hypothesis. Iterate it n times, writing S_n for a sum of n independent steps with probabilities p_i. For i≠j the total variation distance between the law of S_n and its translate by a_i−a_j tends to zero. Condition on which steps are of types i and j and on all remaining step values. The number k in those two types tends to infinity in probability (binomial with success parameter p_i+p_j). Conditional on k, their sum is k a_j plus (a_i−a_j) times a binomial variable with parameters k,p_i/(p_i+p_j). If a_i=a_j the shift is zero. Otherwise the total variation of that binomial law and its unit translate equals its largest atom by unimodality, and tends to zero, for example by Stirling's formula. Bounded convergence proves the assertion for the unconditional mixture. The equation H(s)=E H(s−S_n) thus gives H(s+a_i−a_j)=H(s). All the values H(s−a_i) are consequently equal, and the one-step equation gives H(s−a_i)=H(s). Compactness of translates and their positive lower bound now imply F(s−a_i)/F(s)→1 for each i, and likewise for every fixed finite sum a_w.

For a fixed word w the heat mass in its separated cylinder tube equals ρ_w Z(h/ρ_w) up to exponentially small error. Its ratio to Z(h) tends to ρ_w^d by the last identity. Fixed-depth separated tubes cover all but negligible mass, and their diameters tend uniformly to zero as the depth grows; the normalized heat probabilities therefore converge to the natural measure with cylinder masses ρ_w^d. This step also applies at the constant coefficient a₀/2 by rescaling h. All remaining steps, including |ν̂(t)|≤C/t and the null-support contradiction, follow unchanged. Thus no strongly separated real self-similar Cantor germ of any finite system with d<1 and subquadratic distance-power variance can have the one-point regular output specified in (1).

## One-point versions of the hypotheses

The proof only needs the positive finite limit v(x)/D(x)^α→a₀ as x→x₀ through points outside K, with v=0 on K locally. This already gives both local comparability and the uniform freezing used on the deep cylinder packets. The input need only have a local density continuous and strictly positive at x₀.

One can replace (1) by the total-variation mean condition at the single point

    |Law(Y)−P(y−x₀)dy|(B(x₀,l))=o(l^(d+β))      as l↓0.

For each fixed t the test in (3) annihilates the polynomial globally. Its contribution against the difference measure on the shell |y−x₀|≍2^h R is at most C_{t,M}2^(−Mh) times the total variation in the corresponding ball, because b_j−x₀=O(R). Choose M>d+β and apply dominated geometric summation in a small fixed neighborhood, using the assumed little-o bound; the tails outside cost arbitrary powers of R. This yields precisely (3). Equivalently the proof produces an explicit type of obstruction: since ν is singular, t|ν̂(t)| is unbounded, so for some fixed t the leading term g₀ν̂(t) times the normalized cell mass strictly dominates the uniform C/t error in (4). Along every sufficiently deep chosen cylinder this fixed modulated packet has absolute output pairing at least a positive multiple of R^(d+β).

The exclusion is preserved under a local C¹ diffeomorphism L with a nonzero derivative at x₀. For the image germ L(K) use the physical packet origin L(b_j) and oriented scale L'(x₀)e_j R. Within the bounded scaled packet, L(b_j+e_jR z) differs from this affine approximation by o(R), uniformly. The inverse image of its support remains in the separated model tube after an arbitrarily small enlargement. Upon changing variables, dist(L(x),L(K))=(|L'(x₀)|+o(1)) dist(x,K) on this tube: the candidate nearest preimage lies locally since the local map and inverse are Lipschitz and x→x₀. The input Jacobian similarly tends to |L'(x₀)|. For fixed t, the phase and cutoff errors from linearization are o(1) there, while the Ahlfors upper bound (2) is preserved by local bi-Lipschitz equivalence. The frozen-model limit and contradiction are therefore unchanged.
