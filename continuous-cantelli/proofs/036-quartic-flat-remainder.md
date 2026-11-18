# A quartic variance germ plus an arbitrary flat remainder is impossible

Let a continuous bounded exact variance graph produce a Gaussian output from a standard Gaussian source, as in round 31. At a point x₀ write x=x₀+u and w=y−x₀. Assume some d>0 satisfies, for every integer J≥1,

    v(x₀+u)=d u⁴+O_J(|u|^J)             (u→0 from either side).

Then the exact identity is impossible. Smoothness of the remainder, and an exponential rate of flatness, are unnecessary. The argument works for this quartic Taylor jet. A general quartic contact with a nonzero fifth or higher Taylor coefficient is outside the conclusion.

Put g(u)=p₁(x₀+u), and take a fixed small r so v(u) and d u⁴ are comparable for 0<|u|<r (writing v(u) for the translated graph). In computations below N tends to infinity along even integers. Define

    C_N=binom(2N,N)/4^N,       B_N=(4√d)^N C_N E|Z|^N,   Z~N(0,1).

In particular (B_N)^(1/N)→∞. The proof consists of a lower bound for the quartic reference coefficients and an absolute upper bound with exactly the same exponential base for the error.

## Quartic reference coefficients

Use the whole-line quartic model with source g; removing |u|≥r makes an analytic difference near w=0. For any source graph whose variance is zero at some exterior points, its zero-variance exterior law is supported away from zero. The positive-variance complex kernel is uniformly bounded on |w|<r/4 for |u|≥r, since Re((w−u)²)≥u²/2 and sup_{s>0}p_s(u/√2)≤C/|u|. Thus under exactness the perturbed local output admits an analytic representative near zero.

The full quartic reference conditional on b=√d Z is the image of g(u)du under u↦u+b u². On 1+4bw>0 its density is

    K_b(w)=[g(h_b+√(1+4bw)/(2b))+g(h_b−√(1+4bw)/(2b))]/√(1+4bw),
    h_b=−1/(2b);

on 1+4bw<0 it is zero. For a fixed cutoff A≥1 the |b|≤A contribution is analytic on |w|<1/(8A). In this disc both roots have uniformly bounded imaginary part: the near root is 2w/(1+√(1+4bw)) and the far root is −1/b minus the near root. Hence both Gaussian factors and the inverse discriminant are uniformly bounded.

For |b|≥A the coefficient of w^N, divided by 2(4b)^N(−1)^N C_N, is

    g(h_b)+∑_{k≥1} [g^(2k)(h_b)/(2k)!](4b²)^(−k)
                          binom(k−1/2,N)/binom(−1/2,N).

In fact the Taylor series about the root midpoint leaves only even terms in the square root. All the changes of series converge locally normally. The binomial ratio has modulus ≤3^k: expand (1+z)^(k−1/2) by multiplying (1+z)^k and (1+z)^(−1/2), and use C_{N−j}/C_N≤2^j. The Gaussian Cauchy bound on the unit circle centered at h_b gives |g^(2k)(h_b)|/(2k)!≤e^(1/2)/√(2π) uniformly. Enlarge A until the sum has modulus <g(0)/4 and g(h_b)>3g(0)/4. For even N all large-b conditional coefficients are therefore ≥g(0)4^N C_N |b|^N.

These are also the coefficients of the integrated real asymptotic expansion, to every finite order. Indeed truncate first at |b|≤|w|^(−1/2). For A≤|b|≤|w|^(−1/2) the function of q=4bw is uniformly bounded and holomorphic on |q|≤3/4, so Taylor remainder through order k is ≤C_k|bw|^(k+1), integrable against the Gaussian in b with the asserted order. The discarded Taylor polynomial tails are superpolynomially small in w. The actual discarded density is at most C exp[−1/(4d|w|)] times the integral of exp(−b²/(4d))(1+4bw)^(−1/2) over 1+4bw>0. The last integral is O(1+1/|w|) by separating discriminants below and above 1/2. This proves existence and identification of all asymptotic coefficients. Since ∫_{|b|<A}|b|^N dP≤A^N and the small-b analytic coefficients are O(C^N), the full quartic coefficients satisfy

    [w^N]_{asym} density_quartic ≥ c B_N                  (N even large).  (1)

## Flat perturbations and the full coefficient formula

On |u|<r set δ(u)=v(u)−d u⁴. For real nonzero w the difference of the two local densities is the integral of g(u)[p_{v(u)}(w−u)−p_{d u⁴}(w−u)]. Its order-N Taylor coefficient as a real asymptotic expansion is exactly

    D_N=(1/N!)∫_{|u|<r}g(u)[p_{v(u)}^(N)(−u)−p_{d u⁴}^(N)(−u)]du.      (2)

To prove this, for each fixed Taylor degree k split at |u|=8|w|. On the outer part the two kernels are analytic on |z|≤|u|/4 and bounded by C|u|^(−2)exp(−c/u²), so Cauchy's estimate on the Taylor remainder gives an integrable O_k(|w|^(k+1)) bound after integration. On the inner part differentiate the variance: ∂_s p_s=p_s''/2, whose absolute value is ≤C s^(−3/2) for real arguments. The integral of the actual difference is bounded by C∫_{|u|≤8|w|}|δ(u)| |u|^(−6)du, smaller than every fixed power of |w| by flatness. The missing integrals of all fixed Taylor terms are superpolynomially small by the same off-diagonal Cauchy bound. This proves (2), without a uniform-in-k Taylor remainder assumption.

## Hermite majorant at the required exponential base

For the probabilists' Hermite polynomial He_k one has

    |He_k(z)| exp(−z²/4) ≤ C(k+1)^(1/4) √(k!)           for all real z.  (3)

Indeed h_k=He_k exp(−z²/4) has squared L² norm √(2π) k! by orthogonality, while h'_k=(k h_{k−1}−h_{k+1})/2 has squared norm (2k+1)/4 times that of h_k. Integrating the derivative of h_k² and applying Cauchy–Schwarz gives (3). Consequently

    |p_s^(k)(u)| ≤ C(k+1)^(1/4)√(k!) s^(−(k+1)/2)exp(−u²/(4s)).

Use the mean value inequality in (2) and heat equation at k=N+2. On N^(−1/4)≤|u|<r, every intermediate s≥c/N, so this part of |D_N| is ≤C^N, using √(N!)≥(N/e)^(N/2) and bounded δ.

On |u|≤N^(−1/4), use |δ(u)|≤C|u|¹², so every intermediate s=d u⁴(1+O(N^(−2))) and |δ(u)|≤C|u|^(4+8). Taking the Nth power of this variance comparison loses only a bounded factor. For a fixed C' and β_N=1−C'/N²>0, one obtains

    |D_N|_{small} ≤ C [N^(5/4)/√(N!)] d^(−(N+3)/2)
                    ∫₀^∞ u^(8−2N−2) exp[−β_N/(4d u²)]du
                  ≤ C_d (4√d)^N N^(5/4) Γ(N−7/2)/√(N!).

Here the elementary substitution gives one half times (4d/β_N)^(N−7/2)Γ(N−7/2), and β_N^(−N) stays bounded. The exact Gaussian moment is E|Z|^N=2^(N/2)Γ((N+1)/2)/√π. Stirling's formula gives

    N^(5/4) Γ(N−7/2)/[√(N!) C_N E|Z|^N] = O(N^(−5/2)).                 (4)

It is essential to keep the factor 4√d in the error bound. A generic exponential constant in its place would not justify (4). Thus |D_N|=o(B_N). The absolute majorant is sufficient; no stationary-phase or oscillation assumption is used.

By (1) and (2), adding the local perturbation and analytic exterior terms leaves coefficients ≥c' B_N on every sufficiently large even order. Any analytic representative would have these same asymptotic coefficients by uniqueness along a sequence in the full-measure agreement set, whereas its Taylor coefficients grow at most geometrically. This proves the theorem.
