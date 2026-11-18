# Positive-time analyticity for the symmetric Gaussian expanding exit boundary

Use the C¹ symmetric canonical radius R(t) and the positive slope from 090-critical-slope.md. The canonical future transform late-tail estimate follows from stopping at t+d, as explained in r088-unsigned-canonical-comparison.md. These arguments apply at every finite positive time for the symmetric Gaussian embedding. We prove directly from its global symmetric Volterra equation that R, and therefore its inverse F, are real analytic at every strictly positive radius and time.

Write k(r,a)=∂_a p_r(a)=−a(2π)^(−1/2)r^(−3/2)exp(−a²/(2r)), using the principal branch for Re r>0. The equation for u=R' is

    g(R(t))u(t)=−p'_t(R(t))
       +∫₀ᵗ [k(t−s,R(t)−R(s))+k(t−s,R(t)+R(s))] g(R(s))dR(s).    (1)

Consider initial times t₀ in a small fixed compact neighborhood of a desired positive time T. The real history solves (1), R is bounded and bounded away from zero near T, and the historical derivative is bounded on a somewhat larger compact positive band. The kernel on remote historical times is uniformly bounded with all first derivatives. We extend the unknown future u to a sector S_δ={0<|w|<δ, |arg w|<π/6}, with w=t−t₀. For a bounded holomorphic candidate v(w), define

    Y(w)=R(t₀)+w∫₀¹v(θw)dθ.

This function is holomorphic on the sector, extends continuously to 0 with value R(t₀), and satisfies |Y(w)−Y(θw)|≤K(1−θ)|w| in the ball ∥v∥∞≤K. The same bound applies to differences of such increments with K replaced by ∥v₁−v₂∥∞.

Define the Picard map by the right side of (1) divided by g(Y(w)). Keep the history on [0,t₀] real. Continue the new integral on the segment s=t₀+θw, 0≤θ≤1, replacing R(s) by Y(θw), dR(s) by v(θw)w dθ. The analytic Gaussian g is bounded above and away from zero at all these complex Y for small fixed δ. Here are the uniform kernel estimates giving convergence and contraction.

* In the old co-oriented layer let r=w+t₀−s and a=Y(w)−R(s). For s in the recent real history, |a|≤C_K(|w|+t₀−s)≤C'_K|r|. Consequently |k(r,a)|≤C_K|r|^(−1/2) and |∂_a k(r,a)|≤C_K|r|^(−3/2). The bounded real history density and |Y₁(w)−Y₂(w)|≤|w|∥v₁−v₂∥ make its variation at most C_K sqrt(δ)∥v₁−v₂∥, using |r|≳|w|+t₀−s. The fixed older history remains a finite measure and contributes O_K(δ) to the variation.
* In the new co-oriented layer r=(1−θ)w and a=Y(w)−Y(θw), so the same kernel estimates hold. The difference of a between candidates is ≤|r|∥v₁−v₂∥; variation of the source strength g(Y(θw))v(θw) is O_K(∥v₁−v₂∥). Integrating |w| |r|^(−1/2) over θ gives a total O_K(sqrt(δ)) Lipschitz bound.
* In the opposite-oriented layer near s=t₀, the real part of a=Y(w)+R(s) or Y(w)+Y(θw) is bounded below by a fixed positive number and its imaginary part is O_K(δ). The argument of r is at most π/6, so Re(a²/r)≥c/|r| for δ small. Thus k and its a derivative are bounded even as r→0; its variations in the old and new layers cost O_K(δ)∥v₁−v₂∥. Free and denominator variations obey the same bound. All integral maps are holomorphic by the preceding locally uniform integrable domination (or Morera).

To check the ball condition without any circular bound, take first the test candidate v=0. Its new layer vanishes. In the old layer Y=R(t₀), the constants just used are uniform in t₀ on the chosen band and independent of K, so the image of zero has norm ≤C₀. Choose K>2C₀ plus a uniform bound for the true real derivative. Then fix δ small such that the Lipschitz constant C_K sqrt(δ)<1/2. The map preserves the ball and contracts in its sup norm, and therefore has a unique bounded holomorphic fixed point. The same estimates give contraction when restricted to bounded continuous functions on the real segment, so the real restriction coincides with the true u. The choices of the ball and sector radius were uniform for starting times t₀ near T. Taking t₀<T with T−t₀<δ places T in the interior of the complex sector. Thus R is analytic near T, and R'(T)>0 makes its local inverse F analytic as well.

**Consequence for just the lower wing.** Write h₀(u)=F(b+u)−T for the canonical time increment and θ=limsup_{x→∞} x(H_−−F_−(x)). The unsigned comparison gives

    liminf_{u↓0} h₀(u) log(1/|F_+(b+u)−F(b+u)|) ≥ 2b².

Indeed evaluate the estimate on the radii at time T+h₀(u), then use the locally Lipschitz inverse of the actual high radius (whose derivative is positive by continuity near T). The individual displacement is O(h₀^(5/2)exp(−2b²/h₀)). The reversed fork theorem of 091-fork.md gives the complementary obstruction

    liminf_{u↓0} h₀(u) log(1/(F_−(b+u)−F(b+u))_+)
           ≤ (b−θ/2)²/2.

To see the transfer, h₀(u)∼F_+(b+u)−T and on a sequence realizing the earlier liminf for (F_−−F_+)_+ that difference exceeds exp(−k/h₀) for any k>(b−θ/2)²/2, chosen below 2b². The high-to-canonical displacement is negligible on that sequence. In particular the actual lower wing cannot be of Gevrey order s<2 (or zero-type at s=2) across b, with *no assumption on the high wing at b*. Its difference from analytic F has zero derivatives from the common early interval, and Taylor optimization would contradict the second displayed bound. Quantitatively a local order-two type A for the lower wing must satisfy A≥2F'(b)/(b−θ/2)²; interpret this as ruling out every finite A when θ=2b. The analytic canonical germ can be included in any larger positive order-two type without changing the infimal type.
