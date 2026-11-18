# Strict first time sign at the origin for all Laplace parameters

Let (X,F) have X standard normal, F≥0, EF=1 and E exp(tX−t²F/2)=1 for every real t. A bounded normalized exact graph satisfies these hypotheses. Let R be its signed heat resolvent and S=−∂q R. Then S(q,0)>0 for every q>0. The same conclusion at q=0 follows by continuity for bounded F. The estimates below are uniform in the joint law.

Set z=√(2q), a=z|X| and B=qF. Write g for the standard Gaussian density. The ordinary absolutely convergent signed heat-resolvent formula gives

    z³ S(q,0)=1−E ψ(a,B),     ψ(a,B)=exp(−a−B)(1+a+2B).             (1)

For a≥0 its supremum over B≥0 is W(a)=2exp(−(1+a)/2) for a≤1 and W(a)=(1+a)exp(−a) for a≥1. The two pieces join and decrease. Thus z≥9/16 is covered by a single integral E W((9/16)|X|)<1, certified below.

For 0<z≤9/16, put θ=q, h=1−z² and λ=h exp(−θ). For x≥0 put

    V(x)=(1−exp(−zx))/z + exp(−zx) min(x/h,1/z),
    v(x)=zV(x)=1−exp(−zx)(1−zx/h)_+ .

For a<h, the tangent to B↦ψ(a,B) at θ majorizes ψ for all B≥0. Indeed B↦∂Bψ=exp(−a−B)(1−a−2B) decreases up to where it becomes zero, is positive at θ, and is nonpositive after its zero. Its tangent slope is λ(1−v(|X|)). For a≥h use the maximum W(a), with tangent slope zero. Consequently

    E ψ(a,B) ≤ E ψ(a,θ) + θλ E[(1−v(|X|))(F−1)] + ε,
    ε=E[1_{a≥h}(W(a)−ψ(a,θ))].

Since EF=1, the middle term is θλ E[v(|X|)(1−F)]. For h≤a≤1 integrate −∂Bψ from the maximizer (1−a)/2 to θ to get W−ψ(a,θ)≤exp(−a)θ². For a≥1 integrate from 0 to θ to obtain W−ψ(a,θ)≤exp(−a)[θ(a−1)+θ²]. Dropping exp(−a) gives

    ε/z³ ≤ (z/4)P(|X|≥h/z) + (1/2)E[(|X|−1/z)_+]
         ≤ z² g(h/z)/(2h) + z² g(1/z)/(1+z²).                       (2)

The last line uses the Gaussian Mills bounds g(t)t/(1+t²)≤barΦ(t)≤g(t)/t, which follow by integration by parts. For the constant mark F=1 the signed heat flow is exactly p_s up to s=1 and zero afterwards. Thus (1) gives (1−Eψ(a,θ))/z³=g(0)∫₀¹√s exp(−qs)ds≥(2/3)g(0)exp(−q).

Write A=∫₀∞Vg, M=∫₀∞xVg, B₂=∫₀∞V²g. Apply the entropy inequality with respect to the original marked law reweighted to the positive wing by V(X)/A. Its relative entropy is ≤log(B₂/A²) by Jensen, and the exact real transforms give the bound tM/A−t²E_reweighted F/2≤log(B₂/A²) for every real t. Optimizing in t, multiplying by A, and repeating on the reflected wing yields

    E[F V(|X|)] ≥ U := M²/[A log(B₂/A²)].

Denominators are positive since B₂≥2A² by Cauchy–Schwarz on the half axis. Combining these facts gives

    S(q,0) ≥ (2/3)g(0)exp(−q) + (λ/2)(U−2A)
              − z² g(h/z)/(2h) − z² g(1/z)/(1+z²).                 (3)

All terms have finite limits at zero, with V(x)=2x. This form permits a single interval certificate down to z=0.

## Reproducible bounds

The standard-library script certificate.py produced certificate.json under Python 3.13.18. Run it with python3 -B certificate.py. All arithmetic is integer fixed point with scale N=2⁶⁴ and outward division. For the Gaussian factor it uses 0.39894<g(0)<0.39896, following from the indicated Archimedean rational bounds on π. The elementary exponential and logarithm calls use the correctly nearest-rounded Decimal exp/ln at 96 digits, enlarge by one decimal ulp on each relevant side, and convert the exact decimal endpoint to an outward integer via its rational representation. Every integer argument divided by N is exactly represented in the decimal precision used here.

The rectangles have Δx=Δz=1/512 on 0≤x≤12 and 0≤z≤9/16. In each rectangle (1−exp(−zx))/z increases in x and decreases in z, with limit x at z=0. The other summand of V is bounded by separate endpoint products. The Gaussian g decreases in x. These bounds give lower and upper sums for A, a lower sum for M, and an upper sum for B₂. The omitted tail uses V≤3x and the Gaussian tail moments past 12. The loss terms in (3) are increasing in z; all bounds are rounded in the direction adverse to (3). The large-z check uses upper sums for the product of the two decreasing functions g(x) and W((9/16)x), and bounds W≤2 beyond the last panel.

Actual integer inequalities across 1,769,472 rectangles were

    S(q,0) ≥ 13438013023005694/N                      (0<z≤9/16),
    z³ S(q,0) ≥ 1−18134464667853051311/N              (z≥9/16).

Both right-hand sides are strictly positive. The JSON gives the minimum across each block of 32 consecutive z intervals. No random sampling enters the verification.

## Spatial extensions

On |y|≤1 the spatial Lipschitz constant of S(q,y) is at most 5, uniformly in q>0. Differentiating the signed formula leaves the term

    [−y exp(−z|y|)+E((y−X) exp(−z|y−X|−qF))]/z

and another term of magnitude at most EF. Cancel the constant parts using EX=0. The inequality 1−exp(−z|s|−qF)≤z|s|+z√(F/2), followed by Cauchy–Schwarz, bounds the magnitude of the first term by y²+(1+y²)+√((1+y²)/2). Differentiation in y can equally be read as a weak derivative, so the bound implies Lipschitz continuity. Direct differentiation without cancellation also gives the global spatial Lipschitz bound 4/(e z²), using sup_{v≥0}v exp(−v)=1/e and sup_F F exp(−qF)=1/(eq). Thus (3) and the large-z estimate imply a universal all-q region around the origin; for z≤9/16 take |y|<(13438013023005694/N)/5, and for z≥9/16 take |y|<(1−18134464667853051311/N)e/(4z). R132 enlarges the outer boundary to order log(z)/z at large z.

If f is continuous and f(0)>0 there is a fixed full neighborhood of 0 on which S(q,y)>0 at every q≥0. In fact choose r,c>0 with f≥c on |x|≤3r. For |y|≤r the past term in any fixed signed resolvent derivative is negligible relative to the free derivative as z→∞: on |x|≤3r use exp(−q c) up to polynomial factors, and on |x|>3r one has |x−y|−|y|≥r. Marked powers are controlled by the exponential in F. The low- and compact-q positivity neighborhoods follow from the preceding uniform estimates and continuity.
