# Time and space sign inequalities for every signed slice trace

Use the assumptions and notation of r101-effective-tail.md. Put z=√(2q), q>0. For 0<a<b, d=b−a, let h_d(t)=d exp(−d²/(2t))/√(2πt³) for t>0 and zero otherwise. Then as time measures on (0,H),

    (h_d * μ_a)(dt) − μ_b(dt)
       = ∫_{x>a, f(x)<t} [p_{t−f(x)}(x−b)−p_{t−f(x)}(x+b−2a)] dγ(x) dt ≥ 0.   (1)

To prove this use the R083 positive one-sided resolvent for both traces. Multiplication by the Laplace transform exp(−zd) of h_d gives, for all x>a,

    2 exp(−zd) sinh(z(x−a)) − 2·1_{x>b}sinh(z(x−b))
       = exp(−z|x−b|) − exp(−z(x+b−2a)).

Divide by z, multiply by exp(−qf(x)) and invert on t>0. The right side is nonnegative since x,b>a. Time-local inversion is legitimate by local integrability and equality of Laplace transforms on q>0, taking μ_y zero after H. In particular the Laplace transform of the difference on the full half-line is completely monotone. There can be a nonzero right side also after H, where h_d*μ_a need not vanish. If f(x₀)<t for any x₀>a, continuity and the strict reflected-kernel inequality make the right-side density strictly positive at t (with its extended integral value).

One implication is a necessary prior positive trace between the source origin and any activated outer source. If μ_a restricted to (0,T) is nonpositive, with T≤H, then f(x)≥T for every x>a, and both μ_a and all μ_b for b>a vanish throughout (0,T). Indeed (1) first gives μ_b≤0 before T. For almost every t<T, Q_t restricted to the entire open ray (a,∞) is nonpositive, by Fubini and countably many compactly supported tests. The identity

    ∫ exp(λx) Q_t(dx)=∫_{f>t} exp(λx−λ²(f(x)−t)/2) dγ(x) ≥ 0

then forces the restriction to vanish. If such a slice charged a strictly negative mass in [a+ε,M], that term has magnitude at least const·exp(λ(a+ε)), whereas the full positive part lies on x≤a and contributes at most 2exp(λa) for λ>0. Deterministic zero-age slices can only occur on the at most countable collection of times t with γ{f=t}>0, so do not enter this almost-everywhere conclusion. Thus μ_b=0 for almost every b>a on (0,T). Equation (1) with any one such b now expresses h_d*μ_a, nonpositive on this interval, as a nonnegative measure there. The strict positivity of h_d on (0,∞) forces μ_a=0 there; the strict reflected kernel forces f≥T throughout the ray. Then (1) also gives μ_b=0 for each b>a. The reflected statement holds to the left. This uses individual-slice Laplace positivity and imposes no restriction on the geometry or number of source islands.

There is a smooth family of global sign inequalities at all time orders. The functions

    U_{m,q}(y)=∫ (H−t)^m exp(−qt) dμ_y(t),    m=0,1,…, q≥0,

are strictly positive and strictly decreasing on y>0, with strictly positive second derivatives there. On the left they are strictly positive, increasing and convex. Indeed exp(Hq)L_y(q) is the integral over x>y of

    2 exp(q(H−f(x))) ∑_{k≥0} (2q)^k (x−y)^{2k+1}/(2k+1)!.

All q derivatives have nonnegative coefficients, and the term k=m proves strict positivity at q≥0. For the first spatial derivative replace the sinh kernel by −2cosh(z(x−y)); strict negativity follows by the same series, and the second spatial derivative is

    ∂²_y(exp(Hq)L_y(q)) = 2g(y)exp(q(H−f(y))) + 2q exp(Hq)L_y(q).

Taking m derivatives in q gives strict positivity on y>0, including q=0 (for m=0 the source g(y)>0; for m≥1 the differentiated factor q contributes). Termwise differentiations follow by Gaussian domination on compact q-intervals. Finally multiply by exp(−Hq) to obtain the asserted signs. At q=0 one has the exact moment formula

    U_{m,0}(y)=2 ∑_{k=0}^m 2^k m!/[(m−k)!(2k+1)!]
                         ∫_{x>y}(x−y)^{2k+1}(H−f(x))^{m−k} dγ(x).

These inequalities control integrated traces in all spatial directions. They do not imply the positivity of Q_t at a specified intermediate time.
