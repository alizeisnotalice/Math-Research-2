# Classical slope in the true early interval and continuation of the higher wing

Use r088-early-stefan-symmetry.md. In the continuous single-well case let H_+>H_−>1, T=H_−−1, I(t)=[−b(t),b(t)] be the symmetric true early interval and b=b(T−). Let q_t be the survivor density, taking q_T to be the preactivation limit. Then the early boundary is C¹ at every positive early time and from the left at T, with strictly positive finite speed. If j_t denotes the inward boundary slope of q_t, then

    j_t=lim_{u↓0} q_t(b(t)−u)/u >0,
    b'(t)= j_t/[2g(b(t))],                   0<t≤T.

The derivative at T means the left derivative. On the higher, positive wing the inverse slope matches across the transition:

    lim_{u↓0} [F_+(b+u)−T]/u = 2g(b)/j_T.                        (1)

On the lower wing (1) also holds if limsup_{x→∞}x[H_−−F_−(x)]<2b. No derivatives of the original variance graph are imposed. To avoid ambiguity, (1) is existence of the right tangent at the single radius b; continuity of derivatives on an open exterior neighborhood is not asserted.

**Lipschitz boundary speed away from zero.** For a compact positive time band in (0,T], reflection at the fixed upper barrier b(t) gives q_t(y)≤C(b(t)−y) throughout its support, uniformly in t, by the half-line kernel bound. The mass newly stopped on the right between t and t+h, t+h<T, is ∫_{b(t)}^{b(t+h)}g(x)dx. Any such surviving path must cross the static level b(t) during the next h time units, so this mass is at most

    ∫_{y<b(t)} q_t(y) 2 barΦ((b(t)−y)/sqrt(h))dy ≤ C'h.

The Gaussian density g has a positive lower bound on the bounded swept region. Thus b is locally Lipschitz up to T, after passing to the limit. The stopped mass on either wing has time measure dμ(s)=g(b(s))db(s), with bounded Lebesgue density away from zero.

**Existence, continuity and positivity of the heat slope.** The true-exit density has the free-minus-past representation

    q_t(y)=p_t(y)−∫₀ᵗ [p_{t−s}(y−b(s))+p_{t−s}(y+b(s))]dμ(s).

This has a continuous representative in y at each t>0: near s=t dominate by C(t−s)^(−1/2) times the bounded local time density, and early s use separation in time. This representative is zero at and above b(t). Compute q_t(b(t)−u)/u by subtracting q_t(b(t)+u)=0. For r=t−s, d=b(t)−b(s), the kernel [p_r(d−u)−p_r(d+u)]/u equals a constant times r^(−1/2)exp(−(d²+u²)/(2r))sinh(du/r)/u. Since d≤C r for late s, it is O(r^(−1/2)) uniformly for 0<u<1. On the opposite wing the argument stays separated from zero for late s, hence derivatives are uniformly bounded. Dominated convergence gives

    j_t=−2p'_t(b(t))+2∫₀ᵗ [p'_{t−s}(b(t)−b(s))+p'_{t−s}(b(t)+b(s))]dμ(s).

This quantity is continuous through every positive t up to T from below. Indeed the near-diagonal integral on s∈[t−η,t] is O(sqrt(η)) uniformly on any such compact band, since d=O(t−s) and dμ/ds is bounded. On the complement use dominated convergence and continuity of b. Positivity can also be established without an appeal to parabolic boundary regularity. Fix t>0, choose a short final interval [s₀,t] in the Lipschitz band, with bound K for the speed. The line d(s)=b(t)−K(t−s) lies below the true right boundary. The strip (d(s)−a,d(s)) lies within the true interval if s₀ is close enough to t and a>0 is small. There is positive surviving mass in this strip at s₀: condition at any earlier positive time with surviving mass, and use positive Brownian probability for a tube around a continuous piecewise linear path to a compact interior part of the strip, staying strictly inside the expanding interval. Such a path exists by positivity and continuity of the early radii. Thereafter retain only paths confined to the strip. Subtracting d(s) produces Brownian motion with constant drift in a fixed interval of width a; its killed density is the driftless one times a strictly positive bounded exponential factor on the finite strip. For duration l=t−s₀ the Dirichlet interval heat-kernel sine series is bounded below by a positive multiple of sin(πx/a)sin(πy/a) in translated coordinates 0<x,y<a: choose a sufficiently small so ∑_{n≥2} n²exp(−(n²−1)π²l/(2a²))<1, and use |sin(nz)|≤n|sin z| in the tail of the series. Thus the retained density at y=b(t)−u is ≥c_t u near u=0. This proves j_t>0.

**Recovering the inverse speed.** For 0<t≤T put b_t=b(t), and write h(u)=F_+(b_t+u)−t. The linear onset derived in the main note gives h(u)≥c_t u for small u; its proof applies identically at each t since the surviving law vanishes at least linearly at its upper edge. The density slope implies

    z² e^(−b_t z)∫e^{zy}q_t(y)dy → j_t            (z→+∞).

For t<T every remote wing has a strict tail margin over 1+t; at T the positive wing still does since H_+>1+T. In the positive-frequency future transform the negative wing costs O(exp(−b_t z)), so only the local piece x=b_t+u matters after the displayed scaling. On the positive wing a compact interval separated to the right has a strictly positive time gap; beyond a fixed large radius the strict tail margin and completion of the square give O(exp(−δz²/2)). The local expression is z²∫₀^η exp(zu−z²h(u)/2)g(b_t+u)du. Replacing exp(zu) by one costs O(1/z) by h(u)≥c_tu and exp(zu)−1≤zu exp(zu). Replacing g(b_t+u) by g(b_t) costs O(z^(−2)). Consequently with λ=z²/2,

    λ∫₀^η exp(−λh(u))du → j_t/[2g(b_t)] = C_t>0.

This gives |{u∈(0,η):h(u)≤a}|/a→C_t as a↓0. One short Tauberian argument rescales the positive pushforward of du under h by λ in its argument and multiplies its mass by λ. Its Laplace transforms at all v>0 tend to C_t/v. Exponential testing gives bounded local masses and controls each exponential tail via a test at half the parameter. Every vague subsequential limit therefore has the same Laplace transforms as C_t times Lebesgue measure, hence is that measure. Testing [0,1] proves the claimed sublevel asymptotic. For t<T this sublevel length is b(t+a)−b(t), so the right derivative b'_+(t)=C_t is continuous. Lipschitz absolute continuity gives b∈C¹ and the left derivative at T, by integration of the continuous speed. At t=T the same sublevel asymptotic and monotonicity yield (1), allowing plateaus at higher levels.

For the last low-wing claim, suppose θ=limsup x(H_−−F_−(x))<2b. In the reflected positive-frequency calculation write x=z+s. Choose ε>0 small so that (θ+ε)/(2(1−ε))<b. On |s|≤εz, for large z the Gaussian completion bounds the scaled far integrand by exp(−δz)g(s), after the factor exp(−bz) is included. On |s|>εz beyond a fixed large x choose H_−−F_−(x)<ε²/4 and use the Gaussian exp(−s²/2) to absorb this contribution exponentially in z². The remaining compact region is handled by the time gap as before. Thus the same Laplace inversion proves (1) for the reflected wing.
