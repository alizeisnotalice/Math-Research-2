# Critical parabolic logarithmic contacts

Use signed coordinates at an attained maximum, let g(u)=p₁(x₀+u), and put K_s(y)=∫₀ˢp_t(y)dt. The positive potential S(y)=∫g(u)K_{v(u)}(y−u)du is real analytic since distributionally S'' equals twice the difference of the Gaussian output and source densities. Put L_r=log(1/r).

## Every positive log exponent is excluded under a bilateral upper envelope

Assume v(±u)≤C u²/L_u^k for all sufficiently small u>0 and some k>0, with no spatial regularity assumed beyond continuity. For every A>0,

    ∫_{|u|≤r} v(u)du = o(r³ L_r^(−A)).                          (1)

Thus a positive coefficient c u²/L_u^α on either side, with any finite α>0, is impossible if the opposite side admits any positive-log-power upper envelope (including a zero plateau). The two exponents need not agree.

On |y|≤|u|/4, for j=0,1,2,3 one has

    |∂_y^j K_{v(u)}(y−u)| ≤ C_j |u|^(1−j) exp(−b L_|u|^k).

For j=0 use the heat integral and √v=O(|u|); for j=1 the derivative away from its cusp has absolute value twice a Gaussian upper tail. For j≥2 use K_v''=2p_v off the input mean (and zero convention at v=0) and absorb Hermite powers into one half of the Gaussian exponent. The coefficients a_j=∫g ∂^j K_{v(u)}(−u)/j!, j≤2, are integrable even for k≤1. Beyond a fixed neighborhood use the usual off-diagonal derivative bounds. Split at 8r for |y|≤r. The exterior Taylor remainder through degree two is

    O(r³)+O(r³∫_{8r}^{r₀} u^(−2) exp(−b L_u^k)du)
      = O(r^(5/2)+r²exp(−b' L_r^k)),

by splitting at √r. Restoring the coefficients from |u|≤8r costs o(r² L_r^(−A)) for every A: the j=2 error is controlled by r²∫_{L_r−log8}^∞exp(−b t^k)dt and j=0,1 have an additional positive power in their Gaussian-weighted integrals. Hence the exterior potential equals a₀+a₁y+a₂y²+o(r² L_r^(−A)) uniformly. By monotonicity of K_s and ∫K_s=s, the inner potential is ≤C U(8r)=o(r²). Comparing coefficients against the analytic S at y=0, ±r forces its first three Taylor coefficients to equal the displayed a_j. Its O(r³) remainder and the exterior estimate therefore bound the inner potential by o(r² L_r^(−A)). Integrate it over |y|≤r: for each |u|≤r the Gaussian window has probability at least Φ(1)−1/2 for all 0<t≤v(u)≤r². This proves (1).

## What the unprescribed opposite side must supply

Suppose only v(u)∼c u²/L_u^k on u↓0 (c,k>0) and leave negative u unrestricted. Restrict first to a fixed small symmetric core; omitted means have analytic densities nearby and affect derivative estimates only by geometric Taylor bounds. Define, on the opposite half,

    D(r)=∫_{v(−x)>0,0<x<r₀} g(−x)p_{v(−x)}(x)
             [1−exp(−(2xr+r²)/(2v(−x)))]dx.

Exactness necessarily gives

    D(r) ∼ g(0)c L_r^(−k).                                         (2)

The off-diagonal integral J=∫g(−x)p_{v(−x)}(x)dx is finite by Fatou on the smoothed output. At positive observations the negative-source potential F_− has F_−''(r)=2(J−D(r)) by differentiation off support. Monotone convergence and the uniform first-derivative bound extend this expansion to r=0 from the right, giving

    F_−(r)=F_−(0)+F_−'(0)r+Jr²−2∫₀ʳ (r−s)D(s)ds.

The positive-source potential F_+(r) equals its own off-diagonal jet through order two plus g(0)c r²/L_r^k(1+o(1)): apply the previous exterior estimates on that half; on 0<u<8r the part near r freezes to g(0)v(r) since the regular variation v(u)∼c u²/L_u^k is locally uniform in u/r near 1. More explicitly take |u−r|≤r L_r^(−k/4), where v(u)/v(r)→1 uniformly, dominate by K_(1±o(1))v(r) above and below, and extend the Gaussian kernels to the whole line with relative error o(1). The rest of 0<u<8r has Gaussian distance at least r L_r^(−k/4) and total variance ≤C r²/L_r^k, giving a negligible contribution. Since both corrections are o(r²), matching the three polynomial coefficients against the analytic total, then subtracting its O(r³) remainder, gives ∫₀ʳ(r−s)D(s)ds∼g(0)c r² L_r^(−k)/2. The monotone density argument applied twice (first to the convex primitive, then to its derivative, using difference quotients at (1±ε)r and finally ε↓0) proves (2). A different fixed remote cutoff changes D by O(r). This is a necessary condition for an arbitrary opposite half, not its general impossibility.

## Every regularly varying opposite side is excluded

A separate odd-difference estimate permits an opposite half with no prescribed regular-variation rate beyond index two. For a positive-side source variance w(x)=x²ℓ(x), ℓ positive slowly varying at zero with ℓ→0, let F be its fixed local potential. Then

    F(r)−F(−r)−2rF'(0) = g(0) w(r)(1+o(1)).                       (3)

Indeed ∂K_s is bounded by an absolute constant, yielding F'(0) under the integral. Set q=1/ℓ. Potter's bound for q gives, for fixed 0<δ<1 and small x,r,

    q(x) ≥ C q(r) min[(x/r)^δ,(r/x)^δ].

On x≥8r, the odd Taylor remainder is at most C r³∫x^(−2)exp(−bq(x))dx plus a far-core O(r³). Substitution t=r/x and Potter's bound show this costs O(r² q(r)^(−1/δ))+O(r³)=o(w(r)). Within x≤8r, the linear subtraction costs at most C r∫exp(−bq(x))dx=o(w(r)) by the other side of Potter's bound. The negative-output term K_{w(x)}(r+x) there is negligible: sup_{x≤8r} w(x)≤C w(r) by the positive index two, and the tail distance is ≥r. Finally the positive-output term from |x−r|≤rℓ(r)^(1/4) is g(0)w(r)(1+o(1)) by uniform convergence of the variance ratios and heat-integral monotonicity. Its complement in x≤8r is negligible by a Gaussian tail bound, as above. These estimates establish (3).

For two index-two, ratio-to-quadratic-zero sides, their contributions to the total odd difference have opposite signs and satisfy (3). Analyticity gives O(r³) after its linear part is subtracted, so w_+(r)/w_−(r)→1. In particular if the positive side is the prescribed critical log power, the negative side must have the identical scale/coefficient, contradicting (1). If the putative opposite has regularly varying index >2, it already meets a bilateral log upper bound. For index <2 it has v/u²→∞, contradicting round 20. For index two with ratio not tending to zero, observe that the prescribed flat side forces phase θ=1 by round 21; intermediate ratios bounded between two positive constants have vanishing logarithmic mass on translated fixed annuli, while ratios tending to infinity along a subsequence contradict θ=1 on such an annulus. Uniform convergence of slowly varying functions on each fixed annulus then forces the opposite ratio to tend to zero. Thus (3) applies. The regular-variation argument requires a positive eventually nonzero opposite variance; an identically zero opposite side is covered by (1).

## Every monotone opposite side is excluded

Suppose the positive half has the prescribed contact c r² L_r^(−k) and the opposite variance w(x)=v(−x), x>0, is nondecreasing near zero. Then no exact law is possible, even without regular variation of w. The same proof works under the following weaker pointwise thickness hypothesis. For some fixed 0<a<b and d>0, for every small x there is an interval J_x⊂[ax,bx] of length ≥dx on which w(t)≥d w(x). For monotone w use [x,2x]. This also includes mirror-image monotonicity in the signed coordinate. If w vanishes identically near 0 use (1).

First q(x)=w(x)/x²→0. Otherwise, whenever q(x)≥ε, the thickness interval has q(t)≥dε/b² throughout. The positive side forces the round-21 phase θ=1. On t/x∈[a,b] the middle-ratio set {dε/b²≤q(t)≤1} has vanishing rescaled measure by the finite logarithmic measure estimate, and {q(t)>1} has vanishing rescaled measure by phase 1. This contradicts the length bound.

If B(x)=q(x)L_x^k is bounded, apply (1). Otherwise take r_n↓0 with B(r_n)→∞. With δ=1/4, on 0<t≤r₀ choose x_n attaining the maximum of B(t)exp(−δ|log(t/r_n)|). A maximum exists by continuity and B(t)=O(L_t^k) at zero. Then B(x_n)≥B(r_n)→∞, hence x_n→0. The triangle inequality implies

    B(t) ≤ B(x_n) exp(δ |log(t/x_n)|).

For x=x_n, Q=q(x)→0 and the elementary slow variation of L_t^(−k) yields

    q(t) ≤ C Q max[(t/x)^η,(x/t)^η],    0<t≤r₀,      η=1/2.       (4)

Indeed choose r₀ so small that |d log(L_t^(-k))/d log t|≤δ. No lower bound on q away from x is used.

Consider the odd potential of this half in its own orientation, at y the midpoint of J_x, thus ax≤y≤bx. Subtract its linear term 2yF'(0). The portion from t≥C₀x, C₀>4b+2 fixed, has absolute value at most

    Cx³∫_{C₀x}^{r₀} t^(−2) exp(−c/q(t))dt + O(x³)
        = O(x² Q^(1/η))+O(x³)

by (4) and s=x/t. In t≤C₀x, the linear subtraction is bounded by Cx∫exp(−c/q(t))dt=O(x² Q^(1/η)), by (4) and the substitution s=t/x on (0,x); the fixed annulus [x,C₀x] is exponentially small. In this same inner region sup w(t)≤Cx² Q by (4) and η<2. The negative-observation term K_{w(t)}(y+t) there has Gaussian distance at least ax, so its integral is o(x²Q). The positive-observation term is nonnegative everywhere. On J_x itself, K_{w(t)}≥K_{d w(x)} and g is bounded below by a positive constant. Since the midpoint lies ≥dx/2 from both ends and √w(x)=o(x), this term is ≥c' w(x) for all sufficiently large n. Finally x³=o(x²Q) along our sequence by B(x)→∞. These estimates show

    F(y)−F(−y)−2yF'(0) ≥ c'' x²Q.

The known side has an odd remainder O(x²L_x^(−k)) by (3), uniformly for y∈[ax,bx], and the analytic total has odd remainder O(x³). The orientation reversal changes the unknown side's odd difference by a sign. Consequently Q≤C L_x^(−k), contradicting B(x)→∞. This proves the assertion. The thickness condition is assumed for every small x, so a sparse set of tall spikes with narrowing relative widths remains outside this result.

## A sharp double-logarithmic necessary threshold for an arbitrary opposite side

Let Q(x)=x²/v(−x), with Q=+∞ on the zero set. Under the prescribed single-side contact, put

    d_* = 2 log( sqrt((k+1)/π)/(k c) ).

Then

    liminf_{x↓0} [Q(x)−2(k+1)log L_x−log log L_x] ≤ d_*.        (5)

In fact, suppose Q(x)≥H(L_x)=2(k+1)log L_x+log log L_x+d near zero. Since z^(1/2)e^(−z/2) is decreasing for z>1 and H tends to infinity,

    g(−x)p_{v(−x)}(x)dx ≤ [g(0)+o(1)] sqrt((k+1)/π) exp(−d/2)
                             L_x^(−k−1) |dL_x|.

Fix a∈(0,1). For r^a≤x≤r₀ the bracket in D(r) is at most 2Qr/x, so that this portion is O(r∫_{r^a}^{r₀}dx/x²)=O(r^(1−a)), uniformly over all Q, because Q^(3/2)e^(−Q/2) is bounded. For x<r^a bound the bracket by one and integrate the last display. Thus limsup L_r^k D(r) is at most g(0)sqrt((k+1)/π)exp(−d/2) a^(−k)/k. Let a↑1 and apply (2). Any d strictly larger than d_* is impossible as an eventual lower bound; (5) follows.

For the explicit opposite variance Q(x)=H(L_x), the upper estimate is attained at this order. Given b>1, on x≤r^b the bracket tends to one uniformly, since its exponent is at least H(L_{r^b})(r/r^b)²/2. Integration gives the corresponding lower bound with b^(−k) in place of a^(−k); let b↓1. In particular (5) is sharp for the scalar condition (2) alone. This regular-variation example still fails the exact density law by (3).

## An arbitrary opposite side must attain the universal constant 2

There is a stronger necessary threshold using the analytic potential on both signs. If, on both halves near zero, Q(u)=u²/v(u)≥A log L_|u| for a fixed A>2, then (1) holds for every log exponent. In particular, with one prescribed side c r²/L_r^k,

    liminf_{x↓0} [x²/v(−x)] / log log(1/x) ≤ 2.                    (6)

The constant 2 in (6) comes from summability of off-diagonal second derivatives. Its optimality under the exact law is open.

Choose R large and 1<p<A(1−1/R)²/2. For |u|≥R|y|, j=0,1,2,3 and small |u|,

    |∂^j K_{v(u)}(y−u)| ≤ C v(u) |u|^(−j−1) L_|u|^(−p).          (7)

To obtain (7), use the previous Gaussian derivative formulas and write √v=|u|/√Q. After factoring v|u|^(−j−1), the remaining polynomial factor in √Q is absorbed using the strict inequality p/A<(1−1/R)²/2. For j=0 use K_s(z)≤C√s exp(−z²/(2s)); for j=1 use its Gaussian tail. For j=2,3 the factors are O(Q^(3/2)) and O(Q^(5/2)) respectively.

Let M(r)=∫_|u|≤r v(u)du. Initially M(r)=O(r³). Suppose M(r)≤C r³ L_r^(−β) for β≥0; we prove this with β+p−1 in place of β. In particular a_j=∫g ∂^jK_{v(u)}(−u)/j! exist for j≤2, since the only borderline sum, at j=2, is dominated by ∑_{n≥0}(L_{r₀}+n log2)^(−β−p). For the exterior |u|>Rr, |y|≤r, Taylor remainder through degree two is bounded by

    C r³ [1+∫_{Rr<|u|≤r₀} v(u)|u|^(−4)L_|u|^(−p)du]
        = O(r² L_r^(−β−p)).

Indeed shells [ρ,2ρ] between Rr and √r bound the integral by a constant times the sum of ρ^(−1) L_r^(−β−p), and those beyond √r cost O(1/√r) using the initial area bound. Restoring the near coefficients costs O(r² L_r^(−β−p)) for j=0,1 by geometric summation of r^j ρ^(2−j)L_ρ^(−β−p) over ρ=2^(−n)Rr. For j=2 it costs O(r² L_r^(1−β−p)) by the borderline sum. Consequently the exterior equals a₀+a₁y+a₂y² with uniform error O(r² L_r^(1−β−p)).

The full analytic potential has these same three Taylor coefficients: the inner |u|≤Rr potential is O(sup_{|u|≤Rr}v(u))=o(r²), by kernel monotonicity and ∫K_s=s. Compare at y=0,±r as before. Subtraction from the full potential, and positivity of the inner kernel, now bound its integral over |y|≤r by O(r³ L_r^(1−β−p)). Its lower bound c M(r) follows by the Gaussian-window estimate used for (1). This proves the claimed induction. As p−1>0, taking arbitrarily many finite iterations proves (1), and yields (6) by contradiction. No eventual smoothness or monotonicity of the opposite half is assumed.

## Weak matching and the refined arbitrary-side threshold

Suppose the positive side has the prescribed contact c r² L_r^(−k) and the opposite w(x)=v(−x) obeys only w(x)=o(x²). Then for every continuous compactly supported φ on (0,∞),

    r^(−3)L_r^k ∫ φ(x/r)w(x)dx → c ∫ φ(t)t²dt.                  (8)

Write q=w/x² and F(y)=∫g(−x)K_{w(x)}(y−x)dx for the opposite source in its own orientation, and O(y)=F(y)−F(−y)−2yF'(0). By (3) on the prescribed side and analyticity of the total potential,

    O(y) ~ g(0)c y²L_y^(−k)    (y↓0).

Fixed remote means cost O(y³) after odd linear subtraction. The following estimates prove (8) even for rapidly oscillating q.

Set ρ_j=2^(−j) and m_j=ρ_j^(−3)∫_[ρ_j,2ρ_j] w, for j≥j₀. These m_j are bounded. For y/r bounded above and away from zero, split at C₀r with C₀>4 sup(y/r)+2 fixed. Write O(y) as the positive forward potential ∫_(x≤C₀r)g(−x)K_w(y−x)dx plus an error. On x≤C₀r, both the subtracted linear part and negative-observation kernel are bounded in integral by

    C∫ (r w/x² + w/x) exp(−b/q(x)) dx.

For the first estimate use the Gaussian tail expression for K'_w(−x); for the second use |−y−x|≥x and K_w(z)≤C√w exp(−z²/(2w)). Polynomial powers of 1/q are absorbed by decreasing b>0. On x>C₀r the odd Taylor remainder is bounded by C r³∫ wx^(−4)exp(−b/q(x))dx+O(r³), using the third derivative. Thus for r=ρ_j (or r within a fixed factor), the sum of errors has size at most

    C r² ∑_{n≥j₀} 2^(−|n−j|) ε_n m_n + O(r³),                   (9)

where ε_n=sup_[ρ_n,2ρ_n] exp(−b/q)→0, with slightly smaller b if necessary. For n≥j this follows by integrating each inner weight over its shell, which gives factors (ρ_n/r) and (ρ_n/r)², respectively. For n≤j the exterior gives r/ρ_n; finitely many shells at the splitting radius only change C.

First m_j=O(j^(−k)). If not, choose j_l→∞ with m_{j_l}j_l^k→∞, fix a small δ>0 with 2δ<log 2, and choose n_l maximizing m_n n^k exp(−δ|n−j_l|). Boundedness of m_n guarantees attainment; m_{n_l}n_l^k→∞ and n_l→∞. At such a contact shell j=n_l,

    m_n / m_j ≤ (j/n)^k exp(δ|n−j|) ≤ C exp(2δ|n−j|).

Formula (9) is therefore o(r²m_j) by dominated geometric summation; r³=o(r²m_j) there as well. Now integrate O(y) over y∈[r/2,3r]. Its positive forward portion captures at least a fixed multiple of ∫_[r,2r]w=r³m_j because w=o(x²); the predicted odd difference integrates to O(r³j^(−k)). This is a contradiction. All m_j are thus bounded as claimed.

With that bound established, (9) is o(r² L_r^(−k)) at every small r, uniformly for y/r in any fixed compact subinterval of (0,∞); use the same dominated geometric estimate with m_n≤Cn^(−k). Integrate against φ(y/r)dy. The positive forward portion equals ∫ g(−x)φ(x/r)w(x)dx+o(r³L_r^(−k)). Indeed its total source area in x≤C₀r is O(r³L_r^(−k)) by dyadic summation, and sup_{x≤C₀r} √w(x)/r→0. For the normalized kernels K_w/w, which are uniform averages of Gaussian laws with variance ≤w, uniform continuity of φ bounds the error by o(1) times this total area. Extend φ as zero to the line and choose C₀ past its support. The known limit of O(y), uniform on compact scales, yields (8).

In particular, for a completely arbitrary opposite side the following stronger necessary condition holds (with Q=+∞ if w=0):

    liminf_{x↓0} [Q(x)−2 log L_x−3 log log L_x]
              ≤ 2 log(2/(k sqrt(π))).                              (10)

Indeed, if Q≥2 log L+3 log log L+d eventually, then w=o(x²), and for large L

    Q^(3/2) exp(−Q/2)/sqrt(2π) ≤ [2/sqrt(π)+o(1)] exp(−d/2)/L.

On writing s=log(1/x), (8) and the nonatomic limiting density give
∫_T^(T+1) q(exp(−s))ds ~ c T^(−k); apply (8) to t^(−3) times the indicator of [exp(−1),1] using continuous upper and lower approximations. Hence ∫_T^∞ q(exp(−s))ds/s ~ c T^(−k)/k by summation over unit intervals. The baseline Gaussian density equals g(−x)q(x)Q(x)^(3/2) exp(−Q(x)/2)|ds|/sqrt(2π). For x<r^a this bounds D(r) by [g(0)c+o(1)] 2 exp(−d/2)(aL_r)^(−k)/(k sqrt(π)); the remainder is O(r^(1−a)) as above. Now a↑1 in (2) gives d≤2 log(2/(k sqrt(π))), proving (10). The proof gives a necessary threshold, and does not assert its sharpness for the full exact law.

The necessary condition (2) also determines the unshifted opposite return mass, even for a completely arbitrary continuous opposite side. Set J(r)=∫_(0<x<r,w(x)>0) g(−x)p_w(x)dx. For fixed 0<a<1<b the estimates above imply

    J(r^b)−O(r^(b−1)) ≤ D(r) ≤ J(r^a)+O(r^(1−a)).

For the left hand inequality, the mass omitted from the bracket below r^b is at most ∫₀^(r^b)g(−x)p_w(x+r)dx≤C r^(b−1), because sup_s p_s(z)≤C/|z| for z≠0. For the upper bound use the large-x estimate following (5). Replace r by r^(1/b) and r^(1/a), use (2), then let a,b→1 to obtain

    J(r) ~ g(0)c L_r^(−k).                                         (11)
