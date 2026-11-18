# Equality of Brownian first-exit laws at the local extinction graph

Let f:ℝ→[0,H] be continuous and exact, with X∼γ=N(0,1) and signed heat remainders Q_t as in round 169. Let E be the continuous local last activity time, so f≤E≤H and the spacetime measure Q_t(dx)dt vanishes in the open set {(t,x):t>E(x)}. For t≥H the remainder is zero. Write P_v for the Brownian semigroup including the identity at v=0.

For a Brownian motion launched from (s,x) with s≤E(x), let

    τ_{s,x}=inf{t≥s:t≥E(x+W_{t−s})},    Y_{s,x}=x+W_{τ_{s,x}−s}.

In the second experiment sample X and then launch an independent Brownian motion from (f(X),X). The two experiments give exactly the same joint law:

    Law(τ_{0,0},Y_{0,0}) = ∫ Law(τ_{f(x),x},Y_{f(x),x}) γ(dx).       (1)

Every stop happens by H and obeys τ=E(Y). In particular starting on the graph E means a zero-duration stop. Equality (1) also applies with any continuous bounded majorant e≥f for which Q vanishes on the strict exterior {t>e(x)}, replacing E by e throughout the proof. The quantitative contact conclusions below use E; the same bounds hold for the corresponding contact set of any such e.

## Discrete clipping and passage to first contact

For each deterministic time t, Q_t vanishes as a measure on {E<t}. In fact each smooth spatial test compactly supported there has a neighborhood of its support lying in a column of zero spacetime measure just after t. Its pairing with Q_t is right continuous by bounded-test continuity of the heat semigroup and dominated convergence in the graph integral, including the marks at exactly t. This gives zero at t itself. This step does not presume any trace vanishing on {E=t}.

Use finitely many nested dyadic partitions of [0,H+1], then let their meshes tend to zero. At every grid time v, among the already launched paths still alive just before the decision, stop precisely those at y with E(y)<v. A path launched at v is eligible. Induct on the grid times. Before the first decision the difference of the two free position measures is Q_v. At later decisions the difference before clipping is again Q_v: from each free measure subtract the freely propagated spacetime law of paths stopped at earlier grid times. These propagated submeasures coincide by the induction hypothesis and the Brownian Markov property. Unlaunched sources enter neither free measure in the second experiment. The restriction of Q_v to {E<v} vanishes, so the two position laws clipped at v coincide, including their mass. This proves equality of the entire discrete stopped laws. The terminal grid time H+1 clips all remaining paths.

Couple each sequence of grids on the same launched Brownian path. The first grid hit decreases to the infimum of times in the relatively open strict exterior after launch; openness and density of the nested grids prove this even if the infimum itself is outside that open set. This infimum equals τ almost surely. At τ, continuity and the launch assumption give τ=E(Y), including a launch already on the graph. By strong Markov, the Brownian increments immediately after τ start afresh from zero and almost surely return to zero at positive times arbitrarily close to zero (by immediate oscillation on both sides). At each such time τ+d the position is again Y, so τ+d>E(Y). The argument also covers a tangency at τ=H since the grids extend past H. Brownian path continuity gives convergence of the stopped positions as well as times. Bounded continuous testing of the discrete laws proves (1).

## E is the least barrier for matching the two exit laws

For any nonnegative bounded continuous function e, launch both paths as above, stopping a source above e at once. Their complete first-exit laws for e coincide if and only if e≥E pointwise. For sufficiency use the clipping proof since the strict exterior of e is contained in that of E. For necessity, positive γ-mass with f(x)>e(x) would exit strictly above the graph of e immediately, whereas the path from (0,0) always exits on the graph of e. Continuity then gives f≤e everywhere. Subtracting the common freely evolved stopped measures, as in the next section, shows each Q_t vanishes as a measure on {e≤t}. Thus S lies in the closed region {t≤e(x)} and E≤e. When e≥f, equality of the terminal *spatial* laws is already equivalent to equality of the spacetime laws because both times equal e at the terminal position.

For a quantitative test write Δ_e for the difference of the two stopped spacetime laws. Regardless of equality, the surviving position laws at t are supported on {e>t}. Therefore, with total variation denoting variation mass,

    |Q_t|({e≤t}) ≤ |Δ_e|([0,t]×ℝ) ≤ |Δ_e|([0,∞)×ℝ).

In fact the restriction on the left is the restriction of the freely propagated signed stop law from times ≤t; Gaussian propagation contracts variation. In particular a mismatch has an observable lower bound at any time for which the tested exterior contains signed heat mass. Among all matching barriers, E has the least exit-time distribution in stochastic order and the least exit-position law in convex order: on a common Brownian path from zero the exit under any larger barrier occurs no earlier, and bounded-time optional sampling gives convex order.

## Consequences for individual slices

Let π be the common terminal spatial law and D={x:f(x)=E(x)}. Each source located in D is immediately stopped, so the common exit measure dominates the image of 1_Dγ under x↦(E(x),x). Also γ≼cxπ by conditional Jensen for Brownian motion with bounded duration. From the start at zero one has ∫y dπ=0 and ∫y² dπ=Eτ_{0,0}. If f is nonconstant, E(x)>1 at every finite x by round 169. Thus τ_{0,0}>1 almost surely and the same exit law is obtained by launching a fresh Brownian motion at absolute time 1 from γ.

More strongly, for each t, Q_t is the difference of the two position measures of paths genuinely alive at t under this first-exit rule, counting only sources launched by t. The part of each free position measure coming from paths stopped by t is the common exit law restricted to times ≤t, freely evolved to t. The strong Markov property gives the decomposition including duration zero. In particular

    Q_t restricted to {E≤t} = 0                                 (every t),
    Q_t^−(ℝ) ≤ γ{f≤t} − P(τ_{0,0}≤t),
    ∫₀ᴴ Q_t^−(ℝ)dt ≤ Eτ_{0,0}−1.

The first assertion includes the whole equality level {E=t} even if it has positive spatial measure; it asserts zero measure there, not pointwise continuity of a density. For the second, the negative part is bounded by the surviving past-source measure, whose mass is the right side by (1). Integrate in t and use ∫f dγ=1 for the third.

## Quantitative contact restrictions

For any Borel spatial set K⊂D, its graph measure is dominated by the Brownian exit law. In particular, for each s>0,

    γ restricted to D∩{E=s} ≤ p_s(x)dx restricted to D∩{E=s}.

Indeed a path stopped on that level has terminal position B_s. Thus a positive-measure contact plateau of height s>1 cannot lie where g(x)>p_s(x), equivalently |x|<√(s log(s)/(s−1)).

There is a finite-scale bound useful for any Borel spacetime contact selection. Suppose its spatial positions are in I=[a,b] and its launch heights are in [m,m+d], m>0. Its γ-mass is at most

    ∫_a^b p_m(y)dy
      +∫_{y<a} 2Φ̄((a−y)/√d)p_m(y)dy
      +∫_{y>b} 2Φ̄((y−b)/√d)p_m(y)dy.                              (2)

At d=0 the last two terms are zero. This is the probability, starting at zero, of visiting I at some time between m and m+d: condition on B_m and apply the Brownian reflection principle at the nearer endpoint. Only a subset of these visits can be first exits. Formula (2) sums over arbitrary countable rectangle covers of a spacetime selection; one may take the infimum over such covers. A coarse summand is (|I|+2√(2d/π))/√(2πm). For any compact spatial K in place of I the exact visit upper bound is ∫p_m(y)P_y(T_K≤d)dy, where T_K counts an immediate hit. The bounds apply equally to a Borel subset of the full exit measure, with the left side given by the μ-mass of that selection for μ=Law(τ_{0,0},Y_{0,0}).

For the local form take any Borel K⊂D, a point x with s=E(x)>0, and its lower Lebesgue density θ=liminf_{r↓0}|K∩[x−r,x+r]|/(2r). If θ>0, put d_r=osc_{K∩[x−r,x+r]}E. Then

    liminf_{r↓0} d_r/r² ≥ (π/2) ((θg(x)/p_s(x)−1)_+)².             (3)

Apply (2) with m_r equal to the infimum on the selected set and d=d_r. By continuity m_r→s. On any sequence with d_r/r² bounded, change variables by √d_r in the two exterior integrals. Since 2∫₀^∞ Φ̄(u)du=√(2/π), they sum to 2p_s(x)√(2d_r/π)+o(r); the central term is 2r p_s(x)+o(r). The left side is at least 2r(θg(x)−o(1)). Rearranging gives (3), and unbounded ratios cannot lower the liminf. For a density point of D take θ=1. The use of contact heights alone in d_r can strengthen the comparison with the oscillation over the entire interval. These restrictions bound how much deterministic graph mass Brownian first contacts can carry; they do not assert that the signed remainder is nonnegative everywhere inside the domain.
