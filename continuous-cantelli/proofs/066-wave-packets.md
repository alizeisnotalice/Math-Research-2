# Intrinsic oscillations at stretched-exponentially small amplitudes

Let q>2, p=q/2−1, c>0, ω≠0, g the standard Gaussian density, and suppose D_e=∫g(x)p_{v(x)+e}(y−x)dx=p_{M+e}(y) with M>0, v≥0. Suppose on one side alone, as u↓0,

v(x₀+u)=u^q [c+η(u)H(ωu^{−p})+o(η(u))],    η(u)=exp(−a u^{−pγ}),   a>0, 0<γ<1.

Then every bounded continuous periodic H is constant. The proof also works with η(u)=u^β(log(1/u))^ν exp(−a u^{−pγ}) for arbitrary fixed real β,ν. Fix a nonconstant coefficient of H and choose the sign n such that κ=−nκ₀ℓ>0, where κ₀ is the Fourier unit of H, ℓ=−pω, and Ĥ(n)≠0 (H is real). Take r↓0 with ω r^{−p}→θ modulo its period. Set h=r^{1+p}, d=r^p, a₀=x₀+r, choose γ<b<1, L²=d^{−b}, and λ=h(logη)'(r)=O(d^{1−γ}). Use the complex frequency z=κ+iλ.

Test D_{h²}(a₀+hs) against p_{L²}(s)exp(izs)ds. Absolute integrability holds since each smoothed density is bounded for fixed h and the Gaussian window dominates the real tilt. The resulting kernel at input a₀+ht with rescaled variance T=1+v/h² is

K(t,T)=p_A(t) exp(izt L²/A) exp(−z²L²T/(2A)),    A=L²+T.

This follows from completing the square first at real z, then complex continuation. The exact target is h^{−1}K(−a₀/h,1+M/h²). For real T>0 and arbitrary t, absorbing the term |λt|L²/A into the Gaussian by Young's inequality gives

|K|≤C L^{−1} exp[−t²/(4A)−(κ²−λ²)L²T/(2A)+Cλ²L²].

In particular if T≥L² then |K|≤C L^{−1}exp(−c₁L²), since λ=o(1). If T<L² and |t|≥1/(2d) the contribution is at most CL^{−1}exp(−c₂/(d²L²)), because dL²→0. Since ∫g(a₀+ht)dt=h^{−1}, both the exact target and the entire input tail |t|≥1/(2d) are o(η(r)), as L²≫d^{−γ}.

On I=[−1/(2d),1/(2d)] put T₀(t)=1+c(1+dt)^q. The integral of the baseline g(a₀+ht)K(t,T₀(t)) over I is O(poly(h^{−1})exp(−c₃L²)). For details move the segment to t=x+iεL², with 0<ε a fixed small constant. In the intervening rectangle Re(1+dt)≥1/2 and |Im dt|≤ε dL²=o(1), so the principal powers T₀ and A are holomorphic, |T₀|≤C, and A=L²+O(1). The factor g is uniformly bounded since hL²=o(1). For t=x+iy there, direct use of 1/A=L^{−2}+O(L^{−4}) and izt L²/A=izt+O((|x|+|y|)/L²) gives the modulus bound
C L^{−1} exp[−x²/(8L²)+Cy²/L²−κy+Cλ²L²+C].
Here |λx| and the linear |x|/L² error are absorbed into the quadratic Gaussian. On the upper edge this is integrable with e^{−c₃L²} decay by fixing ε sufficiently small. On each vertical edge, x²/(8L²)≫L². Cauchy's theorem gives the baseline estimate, which is o(η(r)).

Write Δ=T−T₀. With t∈I, both variances stay in one fixed compact interval, and Δ=(1+dt)^q η(r+ht)[H(ω(r+ht)^{−p})+ε_r(t)] with sup_I|ε_r|→0. Taylor expansion in T to first order introduces the integral of O(Δ² sup_T |∂_T² K|). For bounded T, direct differentiation of K gives, for j≤2,
|∂_T^j K|≤(C/L)(1+(|t|/L²)^4) exp[−c₄t²/L²−λt+C|t|/L²+C].
Indeed the logarithmic first derivative is
−1/(2A)+t²/(2A²)−iztL²/A²−z²L⁴/(2A²),
and its T derivative obeys the analogous polynomial bounds. On I, Taylor's theorem in r+ht gives
log(η(r+ht)/η(r))=λt+O(d^{2−γ} t²),
uniformly, including the stated fixed power/log prefactors. In the quadratic remainder a Gaussian envelope therefore absorbs the O(d^{2−γ}t²) term since d^{2−γ}L²→0; completing the square for the remaining λt bounds the remainder by C η(r)² exp[Cλ²L²+C]. This is o(η(r)), since λ²L²=O(d^{2−2γ−b})=o(d^{−γ}) and −logη(r)∼a d^{−γ}. Polynomials from differentiation can be absorbed into the Gaussian with a smaller fixed c₄. After division by η(r), the remainder ε_r(t) in the linear term contributes o(1) by the same bound with the λt factors cancelling.

In the linear term let t=Lw. Uniformly in bounded w, T₀→1+c, g→g(x₀), (1+dt)^q→1, and
ω(r+ht)^{−p}=ωr^{−p}+ℓLw+O(d L²w²).
Moreover, uniformly on such w, the ratio η(r+ht)/η(r) times e^{−λt} is 1+o(1). The factor e^{iztL²/A} after removal of e^{iκt−λt} is 1+o(1). The above logarithmic derivative tends to −κ²/2, and the other factors in LK tend, after removing the modulation, to p₁(w)exp[−κ²(1+c)/2]. The scaled integrands (after exact η-tilt cancellation) are dominated by C(1+|w|^4)exp(−c₅w²). Uniform continuity and these bounds reduce the normalized linear term to
−(κ²/2)g(x₀)exp[−κ²(1+c)/2] ∫ p₁(w)e^{iκLw} H(θ+ℓLw)dw +o(1).
Trigonometric polynomial approximation of H and the Gaussian Fourier transform make the integral converge to exp(inκ₀θ)Ĥ(n). By construction this is nonzero, whereas baseline, target, complement, and nonlinear remainder all contribute o(η(r)). This proves the assertion.

## Critical exponential rate

For η(u)=u^β(log(1/u))^ν exp(−a u^{−p}), the same theorem holds if the real H has a nonzero harmonic at an intrinsic frequency κ=|nκ₀pω|≥1024(1+p)a. In particular this excludes every nontrigonometric-polynomial continuous H for any a>0, and excludes any fixed nonconstant H when a is sufficiently small. Choose n with κ=−nκ₀ℓ>0, let λ=ap, z=κ+iλ, B=32a/κ², L²=B/d. Then λ/κ≤1/1024, κB≤1/32, λB≤1/32, aB≤1/32768 and λ²B≤a/32768. Those inequalities are intentionally loose. The power/log factor has uniformly bounded ratio on I and its ratio tends to one for t=Lw with fixed w. We give estimates at fixed a,κ,p,B and d↓0 so constants independent of d may depend on these fixed parameters.

The kernel and Gaussian variances are as above. Complete the absolute tilt via |λt|L²/A≤t²/(4A)+λ²L⁴/A. For real T≥L² this gives |K|≤CL^{−1}exp(−κ²L²/8); for T<L² and |t|≥1/(2d) it gives |K|≤CL^{−1}exp[(−1/(32B)+λ² B)/d]. Both exponential rates strictly exceed a; factors polynomial in h^{−1} or the logarithm of r are harmless. This controls exterior and target.

For the baseline translate the interval by iκL²/2. On the rectangle one has d Im t≤κB/2≤1/64, |Re t|≤1/(2d) and |t|/L²≤1/(2B)+κ/2, uniformly bounded. Thus T₀=O(1), the source is uniformly bounded, and the logarithm of the modulus of L K(t,T₀(t)) has the upper bound
−x²/(2L²)+y²/(2L²)−κy−λx+O(1).
Indeed, approximation of 1/A by 1/L² costs O(|t|²/L⁴)=O(1), the phase correction is O(|t|/L²)=O(1), and the remaining bounded-variance exponential is O(1) in its logarithm. On the top, the bound is ≤−x²/(4L²)−(3κ²/8−λ²)L²+O(1). On the vertical sides it is ≤[−1/(8B)+|λ|/2]/d+O(1), because y²/(2L²)−κy≤0. The strict gaps over a in both rates again absorb all polynomial factors.

The exponential factor in η has concave logarithm as a function of t throughout I, so η(r+ht)/η(r)≤C exp(λt) there. For bounded real T, |t|/L²≤1/(2B) on I bounds all logarithmic derivatives and phase corrections used in the quadratic remainder, and hence
|∂_T^j K|≤CL^{−1}exp[−t²/(3L²)−λt]  (j≤2).
The integrated quadratic error is ≤Cη(r)²exp(3λ²L²/4)=o(η(r)); the linear error from the uniform o(η) term is o(η(r)) using exact tilt cancellation. These estimates avoid depending on a Taylor estimate with a large second derivative.

At t=Lw the original phase is θ+ℓLw+ω p(p+1)B w²/2+o(1) modulo the period, uniformly on bounded w. The amplitude ratio after cancelling exp(λt) tends to exp[−a p(p+1)B w²/2]. This follows by a degree-two expansion of (1+dt)^{−p}, since d²L³→0. The Gaussian envelope just proved allows passage by compact truncation in w. The prefactor from ∂_T K converges after removing the modulation to −z²exp[−z²(1+c)/2]/2 times the standard Gaussian. Trigonometric polynomial approximation and the Riemann–Lebesgue lemma for each nonresonant harmonic give the exact normalized linear-response limit

−(z²/2) g(x₀) exp[−z²(1+c)/2] Ĥ(n) exp(inκ₀θ)
× [1+B p(p+1)(a−i nκ₀ω)]^{−1/2}.

The square root uses the branch of positive real part. Every factor is nonzero. The critical rate follows by contradiction. The finite threshold 1024(1+p) is sufficient and has not been optimized.

**Quasiperiodic phases.** All of the weak-amplitude conclusions have the following extension. Replace H(ωu^{−p}) by F(Ωu^{−p}), where F is a real continuous function on the m-dimensional torus of coordinate periods 2π, and k·Ω≠0 for every nonzero k∈ℤ^m. For 0<γ<1 the exact identity forces F constant. At γ=1 every nonzero coefficient F̂(k) with |p k·Ω|≥1024(1+p)a is excluded. To verify it, take a subsequence of Ωr^{−p} converging on the compact torus. Fix one chosen wave vector k and use κ=p k·Ω>0, switching k to −k when needed. A continuous torus function is uniformly approximated by trigonometric polynomials. Distinct integer vectors produce distinct intrinsic frequencies −p k·Ω by the assumption, so only the chosen k survives the final large-window average. At the critical rate replace nκ₀ω by k·Ω in the complex square-root factor and replace nκ₀θ by k·θ in the phase. In dimension greater than one an infinite set of integer vectors need not have unbounded projected frequencies, so the critical corollary is correctly expressed through |k·Ω|. The derivative test for algebraic amplitudes extends by the same reasoning, using distinct frequencies under Gaussian convolution.
