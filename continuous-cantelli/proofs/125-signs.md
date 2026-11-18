# Ordinary time derivatives of the single-cut resolvent

Let f≥0 be a bounded continuous exact normalized graph, H=sup f. Write k_q(u)=sinh(√(2q)u)/√(2q), k₀(u)=u, and for y≥0 let R_f(q,y)=2∫₀^∞k_q(u)exp(−qf(y+u))g(y+u)du, with g the standard Gaussian density. On the other half-line use the reflected graph. All q derivatives extend to q=0. The estimates below concern ordinary powers of time in the signed trace; all derivatives in the formulas are with respect to q.

## All cuts for a positive range of q

Assume h=f−1 is nonnegative outside an interval [a,b] and nonpositive inside, with a<b. Exactness gives ∫h dγ=∫xh dγ=0. For every real y the convex function (x−y)₊ minus its affine interpolant at a,b is nonnegative outside [a,b] and nonpositive inside. Consequently ∫(x−y)₊h(x)dγ≥0, and the reflected inequality follows by the same argument. For every y≥0,

    −R_f′(q,y)/R_f(0,y) ≥ exp(−qH)−(2/3) exp(q).                    (1)

Thus strict first-order positivity holds simultaneously at every spatial cut on both wings for 0≤q<log(3/2)/(H+1), without any monotonicity of f within its sign regions. At q=0 the stronger comparison −R_f′(0,y)≥−R₁′(0,y)>0 holds.

Indeed the positive term after differentiating is twice ∫e^(−qf) f k_q dγ≥2 exp(−qH)∫u dγ by the convex test above. The absolute negative term is at most 2∫k_q′ dγ. The constant graph has Q_t=p_t for 0<t<1, hence

    2∫k_q(u)g(y+u)du=∫₀¹ exp(q(1−t))p_t(y)dt.

Differentiating gives ∫k_q′≤exp(q)∫k₀′, and k₀′=u³/3. Finally (∫u³g(y+u)du)/(∫u g(y+u)du)≤2 by the decreasing tilt exp(−yu) from the Rayleigh law at y=0. This proves (1).

For all q>0 and all cuts, the same two-crossing condition implies the zero-order upper comparison R_f≤R₁. For z=√(2q) the signed measure σ=(exp(−q)−exp(−qf))dγ has exterior positive and interior negative signs and annihilates exp(±zx) by exactness. Multiply k_q((x−y)₊) by exp(zx), and view it as a function of v=exp(2zx). It is the hinge (e^(−zy)v−e^(zy))₊/(2z). Subtract its affine interpolant at v=exp(2za),exp(2zb), then multiply back by exp(−zx) and integrate against σ. The subtracted homogeneous solution integrates to zero and the product of signs is nonnegative, as claimed.

## Monotone tails and all q at first order

If f is nondecreasing on [y,∞), the decreasing tilt exp(−qf(y+u)) of the base density proportional to k_q(u)g(y+u) reduces the mean of the increasing function d_q(u)=k′_q(u)/k_q(u). Its base mean is

    τ(q,y)=1−[∫₀¹ t exp(−qt)p_t(y)dt]/[∫₀¹ exp(−qt)p_t(y)dt] < 1.

Thus −R_f′/R_f≥f(y)−τ(q,y). In particular every monotone tail cut with f(y)≥1 has strict first-order sign at every q. Monotonicity of d_q in u follows from d_q=(zu coth(zu)−1)/z² for q>0, and directly from u²/3 at zero. The bound also applies to cuts with f(y)<1 at the values of q for which f(y)>τ(q,y).

## Every fixed finite number of signs at all q on a sufficiently remote wing

Suppose merely that f(y) tends to a finite F>1 as y→+∞. For each fixed integer n there is a finite Y_n such that

    (−1)^j ∂_q^j R_f(q,y)>0          (0≤j≤n, q≥0, y≥Y_n).           (2)

No ordering, derivative, opposite-wing or crossing hypothesis is needed. The direct kernel proof also works for nonexact profiles as an integral estimate. The cutoff can depend on the order; (2) does not establish complete monotonicity at any fixed cut.

Choose δ>0 small enough that f≥1+3δ on a tail. Write d=k′_q/k_q. The product for sinh gives k_q(u)=u∏_{l≥1}(1+2qu²/(π²l²)); it and its derivatives converge locally uniformly. Let c_l=(2u²/(π²l²))/(1+2qu²/(π²l²)). Thus d=Σc_l≤min(u²/3,u/z), z=√(2q), where the second bound follows from coth(v)−1/v≤1. The logarithmic derivative ratios k_q^(j)/k_q sum the products over ordered distinct j-tuples of the c_l. Removing repeated indices from d^j therefore gives

    0≤d^j−k_q^(j)/k_q≤binom(j,2) q^{−1}d^{j−1}   (j≥2),

since max c_l≤1/q. For 0<u≤(1+δ)z, it follows by Leibniz that (−1)^j∂_q^j(k_q exp(−qf))/(k_q exp(−qf)) differs from (f−d)^j by O_{H,n,δ}(1/q). Since f−d≥2δ, all ratios for j≤n are bounded below by one fixed positive number for q sufficiently large. Everywhere their absolute values are bounded by (H+u/z)^j for q>0.

To bound the excluded region under the normalized base density proportional to k_q(u)g(y+u), note that y≥0 adds a decreasing tilt exp(−yu) relative to y=0. The increasing function 1_{u>(1+δ)z}(1+H+u/z)^n therefore has expectation at most its y=0 value. At y=0 the normalized density on u>0 is [g(u−z)−g(u+z)]/[2Φ(z)−1]. For z≥1, a standard Gaussian tail estimate (or Cauchy–Schwarz with the bounded 2n-th moment and Φ̄(δz)) bounds its excluded expectation by C exp(−cδ²q) with an absolute c>0. The further tilt exp(−qf) can increase this expectation by at most exp(q osc_{[y,∞)} f). Increase the spatial cutoff so that the oscillation is at most cδ²/2. The excluded contribution then tends uniformly to zero, proving all signs for sufficiently large q simultaneously.

For 0≤q≤Q fixed, the factor exp(−qf) is bounded above and below. Under the normalized base density the moments of u^{2j} are O_{Q,j}(y^{−2j}) uniformly in q as y→∞: use k_q(u)≤u exp(√(2Q)u) in the numerator, k_q(u)≥u on 0<u<1/y in the denominator, and cancel the common g(y). Since k_q^(j)/k_q≤(u²/3)^j, every nonleading term in the signed Leibniz derivative has vanishing relative mean uniformly, whereas the leading f^j is bounded away from zero. This finishes (2). The same proof permits a remote ray with fixed lower bound 1+3δ and oscillation bounded by a sufficiently small absolute multiple of δ².

## Independent check of the unequal monotone wings closure

R083 §6 gives both limiting heights above one and the unique zero at zero for a hypothetical nonconstant single-valley graph. Let L be the lesser height, T=L−1. R084 applies for every t<T and at the limiting budget, as both limiting heights exceed T by at least one. R088 symmetry on the early band gives f(r)=f(−r) for 0<r<b, where b is the finite positive limiting radius, so ∫xf dγ=0 identifies the two r-weighted exterior means. On the lesser wing f≤L; on the other wing its nondecrease makes Cov(r,f(r))≥0 for r distributed as a standard normal conditioned to exceed b. Both unweighted exterior means are ≤L, contradicting the strict R084 gap 2b[g(b)−bΦ̄(b)]>0. This includes plateaus at any positive height, since the argument cancels the symmetric open interior and the endpoint itself has zero Gaussian mass. It confirms that the general single-valley class reduces to f≡1; the first two sections above also treat genuinely nonmonotone two-crossing profiles.

## Why monotone exterior values ≥1 alone cannot give the second sign

Here is a nonexact monotone test profile for the right-hand kernel, so the obstruction concerns only arguments using its monotone exterior input. Take y=1/10, f_B(x)=1 until x=B and 2 afterwards, for B large. At z=√(2q)=B−1 one has R_f″(q,y)<0. The same holds for continuous monotone smoothings at a sufficiently small transition width, and one can add a sufficiently small positive constant to make the lower plateau strictly above 1.

For the piece on [y,B] the exact formula is

    z R_prefix = e^(−zy)[Φ(B−z)−Φ(y−z)]−e^(zy)[Φ(B+z)−Φ(y+z)].

The remaining [B,∞) piece has the analogous differences of upper tails multiplied by e^(−q). On |z−(B−1)|≤1/2 the full R and its first two z derivatives equal those of e^(−zy)Φ(B−z)/z up to O(exp(−c z²)); this follows by the Gaussian tail bound and the explicit first two derivatives of Φ, including for the omitted term exp(zy)Φ̄(y+z). Since ∂²_q=(∂²_z/z²)−(∂_z/z³), at z=B−1 this yields

    z³ e^(zy) R_f″(q,y) → y² Φ(1)+(2y−1)g(1) <0.

For fixed B,q the integrated ordinary derivative varies continuously under bounded pointwise-almost-everywhere smoothing of f_B by Gaussian domination. Thus the violation persists. Exactness or further structure is essential to extend the all-q outer-cut argument to higher orders at arbitrary fixed cuts. In particular this test is compatible with the far-cut theorem, whose cutoff depends on the entire tail profile and the derivative order.

## A wider all-cut first-order range and finite-order stability

Retain only the two sign crossings and the exactness hypotheses of the first section. At every q>0 all ordinary first derivatives have the desired strict sign whenever

    H(exp(q)−1)/(exp(Hq)−1) > τ(q,0).                               (3)

Here τ(q,0)=1−∫₀¹ t^(1/2)exp(−qt)dt/∫₀¹ t^(−1/2)exp(−qt)dt. In fact, for μ_q(y)=∫₀^∞ k_q(u)g(y+u)du,

    −R_f′(q,y)/(2 exp(−q)μ_q(y)) ≥ H(exp(q)−1)/(exp(Hq)−1)−τ(q,y).

To see this, let K_j(x)=∂_q^j k_q((x−y)₊). In the distributional space variable (D²−2q)K₀=δ_y and (D²−2q)K_j=2j K_{j−1}≥0 for j≥1. For z=√(2q), v=exp(2zx), writing K_j=e^(−zx)U_j(v) gives (D²−z²)K_j=4z²exp(−zx)v² U_j″(v) for smooth portions and preserves the derivative jump sign. Hence U_j is convex and the earlier generalized secant argument gives ∫e^(−qf)K_j dγ≤e^(−q)∫K_j dγ for every j. The elementary concavity of −s log(s)/q between s=exp(−qH) and 1 gives f exp(−qf)≥H(1−exp(−qf))/(exp(qH)−1). Integrate against K₀, and use the upper comparison against K₁ for the negative term. Finally τ(q,y)≤τ(q,0) because its base-mean integrand d_q increases in u. The left side of (3) decreases from 1 to 0 for H>1, while τ(q,0) increases from 2/3 to 1; thus the interval has a unique upper endpoint. For each finite Q this criterion includes [0,Q] if H is sufficiently close to one from above.

For any fixed n≥1 and Q≥0, an upper bound on H−1 alone also gives positivity of every order through n on [0,Q] at every cut, even when f reaches zero. A sufficient bound is

    0≤H−1≤e^(−2Q)/[4(2n+1)4^n(n+Q+1)].                            (4)

For q≥0 each K_l has a nonnegative series in the convex hinges (x−y)₊^(2m+1). Consequently ∫h K_l dγ≥0 by the ordinary affine secant at the crossings, and ∫|h|K_l dγ≤2(H−1)∫K_l dγ. The bound (4) implies f≤2. For m≤n, q≤Q, the functions v^m exp(−qv) on 0≤v≤2 have Lipschitz constant at most 2^n(n+Q+1). The constant-graph transform gives ∫K_l≤∫K₀ for each l by differentiating its representation with exp(q(1−t)). Hence by Leibniz and the preceding weighted L¹ bounds the error in the jth signed derivative from the constant graph, j≤n, is at most 4(H−1)4^n(n+Q+1)μ_q(y). On the other hand the jth signed derivative of R₁ is at least 2 exp(−2Q)μ_q(y)/(2n+1): tilt the time law with density exp(−qt)p_t(0) on [0,1] by the increasing factor exp(−y²/(2t)) to see its t^j mean is minimized at y=0, where it is ≥exp(−Q)/(2j+1). Applying (4) retains at least one half of this lower bound. Reflection proves the same result throughout the line. This is uniform over the class with the given upper bound on H; it uses no lower bound on f.
