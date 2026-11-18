# Every convergent quartic Taylor jet is excluded, even after a flat perturbation

Let g(u)=p₁(x₀+u). Suppose a bounded continuous variance graph yields the exact Gaussian output. If d>0 and H is real analytic near zero, H(0)=1, the graph cannot satisfy

    v(x₀+u) = d u⁴ H(u)² + O_J(|u|^J)     for every J≥1 on both sides.

Thus a C∞ exact graph with a quartic contact cannot have a convergent Taylor series at that point (even if it does not coincide with that series). Indeed any convergent real quartic Taylor series with positive leading coefficient d admits this analytic square root factorization after dividing by d u⁴. This includes H(u)=(1−a u)⁻¹ for arbitrary real a, and arbitrary allowed higher coefficients if the first correction is of order six or higher. The conclusion imposes no regularity on the flat remainder beyond its value bounds. Divergent Taylor series remain outside the theorem.

Use N even, C_N=binom(2N,N)/4^N, B_N=(4√d)^N C_N E|Z|^N for standard normal Z. Write a=H'(0), A=|a|/√d. We prove a lower bound c B_N exp(A√N/2) for the asymptotic Nth density coefficient of the local reference, and an o(B_N exp(A√N/2)) estimate for the flat difference. Exact Gaussian output would force the local density to have an analytic representative: exterior mean positions |u|≥r give a complex Gaussian kernel bounded by C/|u| on |w|<r/4; any zero-variance exterior law is supported away from this disc. Thus coefficients growing faster than geometrically are impossible.

## Conditional density of the analytic reference

Take a fixed small r>0 so H is holomorphic on a neighborhood of [−2r,2r], H>0 and (uH(u))'>1/2 on this real interval. Restrict the reference to |u|<r and condition on b~N(0,d); the output is w=u+b u²H(u). Introduce e=1/b and F_e(u)=u²H(u)+eu. For |e| small there is a unique critical point h(e) near zero. Analytic implicit differentiation and expansion at zero give

    h(e)=−e/2+O(e²),
    F_e(h(e))=−e²/4−a e³/8+O(e⁴),
    κ(b)=−e/F_e(h(e))=4b−2a+O(1/b),     κ'(b)=4+O(1/b²).

Indeed writing h=−e/2+O(e²) in u²+eu eliminates its cubic correction, leaving a h³. Factor F_e(u)−F_e(h)=(u−h)²D(e,u) with D analytic, D(0,0)=1. The map ζ=(u−h)√D(e,u) has a jointly holomorphic inverse u=U_e(ζ) on a fixed small bidisc in e,ζ, real and increasing for real small arguments. This follows from the analytic inverse theorem at (e,ζ)=(0,0); take the square root positive there. For large real |b| and small real |w| any root in |u|<r must lie in this inverse coordinate patch. Indeed |u|²H(u)≤|e|(|w|+r) at a root, H≥1/2, uniformly forcing |u|→0. Conversely the two coordinate roots belong to |u|<r whenever their real square root exists and |e| is small. Set δ_b=−F_e(h(e))>0, ζ²=δ_b+ew=δ_b(1+κ(b)w).

For 1+κ(b)w>0 the conditional density is

    K_b(w)=G_b(w)(1+κ(b)w)^(−1/2),
    G_b(w)=|e|/(2√δ_b) · [g(U_e(ζ))U_e'(ζ)+g(U_e(−ζ))U_e'(−ζ)].

It vanishes for 1+κ(b)w<0. The bracket is even and holomorphic in ζ, so it is holomorphic in ζ². In particular G_b is holomorphic and bounded uniformly in a fixed physical disc |w|≤r₀, for all |b|≥B with B sufficiently large, and G_b(0)→2g(0). For even N divide the coefficient of K_b by C_N κ(b)^N. The leading term is G_b(0); by Cauchy's bound and C_{N−j}/C_N≤2^j the remaining sum is bounded by C∑_{j≥1}(2/(r₀|κ(b)|))^j. Enlarging B makes every large-b coefficient ≥g(0) C_N |κ(b)|^N for all even N.

These conditional coefficients integrate to the asymptotic coefficients to each finite order, even with the moving fold singularity. On B≤|b|≤|w|^(−1/2), expand in q=κ(b)w, using uniform holomorphy and boundedness on |q|≤3/4 to bound the degree-k remainder by C_k|bw|^(k+1). The Gaussian b-moments are finite, and extending this Taylor polynomial to the full b-tail costs a superpolynomial error. On |b|>|w|^(−1/2), the actual density is bounded in the b-integral by C exp(−1/(4d|w|)) times the integral of exp(−b²/(4d))(1+κ(b)w)^(−1/2) over its positive-discriminant region. This last integral is O(1+1/|w|): on discriminants ≥1/2 use the Gaussian bound, while below 1/2 change variables to the discriminant on the two b-rays, using κ'(b)≥2. Thus it contributes no finite-order term.

## The spatial and conditional cutoffs for bounded b

Ensure B above is large enough that |r H(±r)|>2/B. After fixing B choose 0<r'<r/4 very small such that |(b u²H(u))'|≤1/2 for |u|≤2r', |b|≤B, and |r'H(±r')|<1/(2B). For bounded b and |u|≤r', u↦w is monotone, and its endpoint images lie beyond ±r'/2. For small |w| the inverse at u=0 and its Jacobian are uniformly analytic in w for b∈[−B,B] by the analytic implicit theorem and compactness. This part of the bounded-b integral is therefore analytic.

On the remaining compact real annuli r'≤|u|≤r, first integrate b rather than fixing a branch. The conditional cutoff gives exactly

    ∫ g(u) ρ_d((w−u)/(u²H(u)))/(u²H(u))
          · 1_{|(w−u)/(u²H(u))|≤B} du.

The integrand without the indicator extends jointly holomorphically to a neighborhood of (w=0,annuli). At w=0 equality in the indicator occurs once strictly inside each of the two annuli: uH(u)=±1/B, with nonzero derivative. Equivalently, the crossing derivative in u of −1/(uH(u)) is nonzero there. The four fixed endpoints ±r,±r' are separated from equality by the previous choices. For small real w the two simple crossing positions are analytic functions of w and no other crossings appear, by compactness. The integral is a sum of integrals of a jointly analytic function with fixed or analytic endpoints. Each is analytic in w (split off a small neighborhood of each varying endpoint and integrate the convergent power series there). This proves the bounded-b assertion with both cutoffs in place.

## Size of the reference coefficients

Over |b|≥B the even-order contributions all have positive sign. Select b with sign opposite a (either sign if a=0) and |b| in [b_*,b_*+1], b_*=√(dN). Uniformly there,

    |κ(b)|^N ≥ (4|b|)^N exp(A√N/2−C).

This follows by taking logarithms of |κ|/(4|b|)=1+|a|/(2|b|)+O(|b|^(−2)). The chosen interval carries a fixed positive fraction, independent of N, of the full Nth absolute Gaussian moment. Indeed on the positive half-axis the log-density N log s−s²/(2d) is maximized at b_*, with second derivative ≤−1/d everywhere, and ≥−2/d on [b_*,b_*+1]. A Gaussian upper bound on the total and a lower bound on this unit interval prove the claim. The bounded-b analytic piece contributes at worst −C^N, negligible here. Thus the reference Nth asymptotic coefficient is ≥c B_N exp(A√N/2).

## The flat perturbation has a smaller coefficient

Write v₀=d u⁴H(u)² and δ=v−v₀. The argument for the quartic jet applies verbatim to show that the difference of the two local output densities has the real asymptotic coefficient

    D_N=(1/N!)∫_{|u|<r} g(u)[p_v^(N)(−u)−p_(v₀)^(N)(−u)]du.

Indeed on |u|≤8|w| the actual density difference integrates to ≤C∫|δ||u|^(−6), smaller than every power; on the complement Taylor remainders are controlled on the discs |z|≤|u|/4 by C|u|^(−2)exp(−c/u²). For the probabilists' Hermite functions, L² orthogonality and the first derivative recurrence give |He_k(t)|exp(−t²/4)≤C(k+1)^(1/4)√(k!). Consequently

    |p_s^(k)(u)| ≤ C(k+1)^(1/4)√(k!) s^(−(k+1)/2)exp(−u²/(4s)).

Differentiate in the variance between v and v₀ (k=N+2), multiply by |δ|/(2N!). On N^(−1/4)≤|u|≤r every intermediate variance is ≥c/N, so this integral is O(C^N). On |u|≤N^(−1/4), use |δ|≤C|u|¹². Each intermediate variance is v₀(1+O(N⁻²)); its high power ratio is bounded and the Gaussian factor can be bounded with β=1−C'/N². Put

    z=1/(√d u),      y=z/H(u),      z=y+a/√d+O(1/|y|).

The last relation and |dz/dy|≤C follow by the convergent Laurent expansion and inversion at infinity on either real ray. Then u²/v₀=y² and 1/√v₀=√d |zy|. The transformation |δ|du≤C_d |z|^(−14)|dz| gives

    |D_N|_{small} ≤ C_d d^(N/2) N^(5/4)/√(N!) ·
          ∫ |z|^(N−11)|y|^(N+3) exp(−β y²/4) dy.

The integral runs over large positive and negative y only. For h=|y| sufficiently large there is a fixed K≥0 with |zy|≤(h+A/2)²+K. Set t=√[(h+A/2)²+K]. Then t≥h, dh/dt≤2, and h²=t²−At+O(1) uniformly. Writing the power as |zy|^(N−11) h^14 bounds it by t^(2N−8). Extend the resulting positive integral to t∈(0,∞), giving

    ∫ |z|^(N−11)|y|^(N+3)e^(−β y²/4)dy
       ≤ C ∫₀^∞ t^(2N−8) exp(−βt²/4) exp(βAt/4)dt
       ≤ C exp(A√N/2) (4/β)^(N−7/2) Γ(N−7/2).

For the last inequality, the log base weight (2N−8)log t−βt²/4 is maximized at t₀=√[(4N−16)/β]. Its second derivative is ≤−β/2 everywhere and ≥−β on [t₀,t₀+1]. Hence the expectation of exp(β A(t−t₀)/4) under its normalized density is bounded uniformly by a Gaussian moment divided by the mass estimate on that unit interval. Also β A t₀/4≤A√N/2. The integral of the base weight is the usual gamma integral, and β^(−N) is bounded.

Thus, by Stirling's formula as in the pure quartic case,

    |D_N| ≤ C^N + C_d B_N exp(A√N/2) N^(−5/2).

In particular the flat difference cannot cancel the reference coefficients. Their Nth roots diverge, whereas any analytic representative must have at most geometric Taylor coefficients. This contradiction proves the theorem.

For the explicit H(u)=(1−a u)⁻¹, the root discriminant is exactly 1+(4b−2a)w+a² w², independently confirming the linear term −2a in κ. For H'(0)=0 the subexponential factor is identically one, so the sixth and higher analytic coefficients do not alter the required base.
