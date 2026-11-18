# Local mass in the recentered Esscher law and analytic rigidity with sublinear derivative growth

Let 0≤f≤M be bounded and real analytic on R and satisfy
F(t)=∫exp(tx−t² f(x)/2)γ(dx)=1 for every real t. Hence Eγ f=1 by differentiating twice. Suppose |f'(x)|≤C(1+|x|)^m for some m≥0. There exist K,c>0 such that for all sufficiently large t>0 and

    K log(t)/t ≤ r ≤ 1,

the probability ν_t(dx)=exp(tx−t² f(x)/2)γ(dx) satisfies

    ν_t(|X−t|≤r) ≥ c r.                                           (1)

Reflection gives the negative-tail version. In consequence, if |f'|=o(|x|/log|x|) at both infinities, f is identically 1. This rules out analytic functions with globally bounded derivative, or more generally with polynomial derivative growth of any degree p<1. No prior tail limits and no uniform holomorphic strip are assumed.

**Deficits and small source positions.** Put g=f−1 and let γ_t=N(t,1). Then ν_t=e^(−t²g/2)γ_t. If h=−g(x)>0, on [x−s,x+s] with s=min(1,h/(2C'(1+|x|)^m)) we have g≤−h/2. Positivity at t=x gives 1≥c₀s exp(x²h/4). For |x| large, either h≤|x|^(−2) or −log(s)≤(m+2)log|x|+O(1). Thus

    (−g(x))_+ = O(log|x|/x²).                                    (2)

Analytic f≥0 can have no zero outside the origin: any such zero x₀ has f(x₀+u)≤C u² and, with t of the sign of x₀, a window of size 1/|t| would make F(t) diverge. If f(0)=0, it must be quadratic. Otherwise analyticity and nonnegativity give f(x)=O(x⁴), or f≡0, and on [2 log t/t,3 log t/t] the integrand in F(t) is ≥c t² for large positive t, a contradiction. At a quadratic zero f(x)≥c x² near 0, hence the local ν_t mass is bounded by C∫exp(tx−ct²x²/2)dx=O(1/t). On a fixed remaining compact set f is bounded below. For x≤−1 use exp(tx)≤e^(−t); for R≤x≤t/8 choose fixed R large so (2) ensures f≥1/2 there and use exp(tx−t² f/2)≤exp(−t²/8). Together,

    ν_t(X<t/8)=O(1/t).                                             (3)

**Truncation and Gaussian smoothing.** Write W=X−tf(X). The entire identity at t+iu gives

    E_ν exp(iuW + u² f/2) = 1.

Set η=σ²=80 log(t)/t², and retain f≤1+η. On this event let
S=W+√(1+η−f+σ²)Z with an independent standard Gaussian. Its (subprobability) law has characteristic function

    exp(−σ²u²/2) E_ν[1_{g≤η} exp(iuW−(η−g)u²/2)].

The same expression without the indicator equals exp(−(1+η+σ²)u²/2). For |u|≤t/2, the absolute difference is at most exp(−3t²η/8)=t^(−30), because ν_t=e^(−t²g/2)γ_t and g>η in the omitted expectation. For |u|>t/2 both characteristic functions are bounded by exp(−σ²u²/2). Fourier inversion gives a uniform density discrepancy O(t^(−9)/log t)+O(t^(−29)); here ∫_a^∞ exp(−σ²u²/2)du≤exp(−σ²a²/2)/(aσ²). In particular the density of S is bounded below by an absolute constant on [−1,1] for large t.

On the retained event with X≥t/8, (2) gives −g≤C'log t/t²; hence |tg|≤C''log t/t and the variance added to W is ≤C''log t/t². Since W=X−t−tg, there is a fixed A such that the retained joint mass of |S−(X−t)|>A log t/t is O(1/t): (3) bounds the low-source part, and a Gaussian tail bound ≤2t^(−2) handles the other part. For r≥2A log t/t the probability ν_t(|X−t|≤r) is at least the retained mass of |S|≤r/2 minus O(1/t). The density lower bound gives (1).

**Attainment and rigidity.** If f≥1+δ throughout [t−r,t+r] then ν_t(|X−t|≤r)≤exp(−δt²/2). By (1) such an interval cannot have r≥K log t/t at large |t|. If along a sequence |t|→∞ we had f(t)>1+2δ under the sublinear derivative assumption, the mean-value theorem would give f≥1+δ throughout those intervals. Therefore limsup f≤1 at both tails; (2) gives liminf f≥1. A nonconstant f with Eγ f=1 takes a value >1, so it would attain its supremum in a compact set, contrary to the analytic-maximum theorem (round 23). Hence f=1. Contrapositively, any nonconstant exact real-analytic graph, with no derivative growth assumption, must satisfy limsup_{|x|→∞}|f'(x)|log|x|/|x|>0 (possibly infinite).

The argument up to (1) only uses regularity in the deficit width bound (2), plus absence of zeros away from zero and the local quadratic lower bound if zero is attained. Analyticity is sufficient for the latter assertions under the exact identities.

The localization also gives an actual near-unit point in each such window. Within the retained event one already has f≤1+η; hence the lower bound for P(|S|≤r/2) minus the exceptional probability applies to ν_t(|X−t|≤r, f≤1+η). For r≥K log(t)/t it is strictly positive, with η=80 log(t)/t² as above.
