# Independent review of the controlled quartic jet argument

Reviewed round035-general-quartic-controlled-jet.md and its two prerequisite coefficient-identity notes. I find the signed exterior lower bound and the stated strict-Gevrey and critical-zero-type exclusions valid. The near-diagonal and exterior contributions both have a positive lower bound on the selected even orders. Their normalizations agree: each half-axis exterior has prefactor L_k/2 after normalizing its signed Hermite measure.

## Direct checks of the sign and scale

Integrating Rodrigues by parts k times against z^(2k)exp(az) leaves a positive base density X^k exp(−X²/2) and the polynomial P_k(aX)=∑ binom(k,j)k!/(k+j)! (aX)^j; the edge terms vanish at zero. The Laguerre recurrence, Jacobi matrix and Gershgorin bound in the draft give all reciprocal roots O(1/k), and the exact first two sums are

    ∑r_i^−1=k/(k+1),
    ∑r_i^−2=k(2k+1)/[(k+1)²(k+2)].

A strengthening is T_k(a)→1 for fixed complex a. Under the normalized base density Y=X−√k tends to a centered Gaussian of variance 1/2: on bounded Y the quadratic expansion of k log X−X²/2 gives −Y²+o(1), and the global Gaussian majorant justifies normalization and passage to the limit. On bounded Y the log product is aX−a²+o(1); the cubic remainder is O(k^(3/2)∑r_i^−3)=O(k^(−1/2)). Thus the full integrand tends to exp(2aY−a²), whose limiting expectation is 1. For bounded complex a and X∈[√k/2,2√k], the modulus is ≤C exp(C|Y|). Outside this central strip the product is bounded by exp(CX), and the shifted Gaussian majorant gives exponentially small tails despite the exp(O(√k)) centering. This justifies the asserted local boundedness, positive real limit, and bounded derivatives by Cauchy's formula.

In the exterior variable z=1/(x√w_σ(x)), the signed weight is exp(kE_σ(1/z))G_σ(1/z), with G including the positive Jacobian reciprocal. Expanding at z=2√k+Z gives the precise linear tilt −a_σ Z/4 and shift a_σ√k/2+b_σ/4. The O(k^(−1)) weighted Taylor remainder costs O(k^(−1/2)) because |σ_k| is dominated by C√k times a Gaussian centered at 2√k. The quadratic correction costs O(k^(−1/2)) against the signed measure using bounded T_k derivatives, not its absolute mass. On the intermediate tails the normalized weight is ≤C exp(C|Z|). The far tails are O(exp(−ck)) after choosing the fixed spatial cutoff sufficiently small. Uniform sign therefore survives separately on each half-line. In fact the exterior contribution from side σ is

    (L_k/2) exp(a_σ√k/2+b_σ/4) [g(0)+o(1)].

I also checked the diagonal reindexing, beta coefficient formula, inverse coordinate residue change and even-part averaging. Standard pointwise Gevrey conventions shift the exponent on Taylor coefficients by one relative to derivatives.

## Additional critical positive-type range

The argument excludes every critical Taylor type
τ=limsup_j(|[u^j]w|/√(j!))^(1/j)<2√(e c), c=w(0). No optimality or equality-case claim is made.

Indeed for h=(4ck)^(−1/2), the product for β_i/h^(2i) is bounded by a universal constant uniformly in 0≤i<k/2. A factor exceeds 1 only for j(j+1/2)<k/2; there are O(√k) such indices and each factor is at most 1+O(1/k). The other factors are ≤1. The approximation to 1 for i≤T√k remains as in the draft.

Choose τ<A<2√(ec) and a fixed ε>0 small enough that (1+3ε)A<2√(ec). Then |[u^j]b|≤C_A A^j√(j!) for b=w/c. Log-convexity of j! gives for 2≤j≤k

    √(j!) ≤ √2 K_k^(j−2),     K_k=(k!/2)^[1/(2(k−2))] ∼ √(k/e).

Thus (1+3ε)h A K_k stays strictly below 1. The absolute Taylor sum after the linear term on this disc is O(h²) by a geometric series, and the derivative is uniformly bounded there by the analogous differentiated series. The large-i tail (i>T√k) of both diagonal sums is bounded by C(1+ε)^(−2T√k) exp((1+ε)λ√k+O(1)) using the uniform beta bound and a majorant on the (1+ε)h circle. Choose T large enough. The inverse map is uniformly defined on |z|≤(1+2ε)h using the (1+3ε)h circle; the analytic inverse-series tail from i≥k/2 costs exp(O(√k))(1+2ε)^(−k). Thus the positive real-axis evaluation and contradiction are unchanged.

## Independent computation

r035_independent_signed_audit.py checks four exact Rodrigues identities, 128 exact reciprocal sums, and the beta maxima in exact rational arithmetic up to k=1024. A separate quadrature of the transformed probability shows T_1024(−1)=1.0059953 and T_1024(1)=0.9955004. These floating computations illustrate the limiting proof; they do not establish the sign for arbitrary k.
