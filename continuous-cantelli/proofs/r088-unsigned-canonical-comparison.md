# Propagation from the lower to the higher outer wing across distance 2b

Assume the single-well setting in 090-critical-slope.md, H_+>H_−=1+T, b=b(T−)>0. Write R_1(t) for the actual higher/right inverse and R_2(t) for the symmetric expanding Brownian exit boundary with standard Gaussian terminal law. The latter exists by Cox and Peskir, Theorem 1; their centered-target UI result makes the stopped Brownian motion a uniformly integrable martingale. Gaussian full support precludes spatial jumps and Gaussian absence of atoms precludes fixed-radius plateaus (a fixed boundary segment would receive positive exiting mass). The finite positive early radius cannot diverge immediately after T; more generally a finite-time explosion of the symmetric boundary would force all paths to have stopped by that time, whereas confining a positive surviving mass in a fixed interior interval gives positive survival beyond it. Its early truncation stops spatial density g on its swept interval. By the enthalpy uniqueness the true and canonical curves agree through T. [Cox and Peskir, Embedding laws in diffusions by functions of time](https://personalpages.manchester.ac.uk/staff/goran.peskir/embedding.pdf)

For all sufficiently small h>0, with a finite constant uniform in h,

    |R_1'(T+h)−R_2'(T+h)| ≤ C h^(−1/2) exp(−2b²/h),
    |R_1(T+h)−R_2(T+h)| ≤ C h^(3/2) exp(−2b²/h).                   (1)

No ordering or sign restriction on the difference is assumed.

**The unilateral flux identity in the signed region.** Choose a compact positive time band containing T in its interior and lying below H_+−1. The one-sided resolvent in round 83 gives supp Q_t⊂(−∞,R_1(t)] for the full signed remainder throughout this band: for an observation to the right of R_1(t) its resolvent uses only the right tail, where the tail margin is strict. The weak preactivation limit Q_(t−) is similarly supported on (−∞,R_1(t−)]. All past means in {f<t} are at most c=R_1(t−). Subtract its reflection about c (zero on the left) in the free-minus-past formula. On y<c, each age-positive past kernel difference p_r(y−x)−p_r(2c−y−x) is nonnegative for x≤c. Any age-zero source on the left gives the same measure inequality. Thus

    Q_(t−)(dy) ≤ [p_t(y)−p_t(2c−y)]dy on y<c.

There is no atom at c since the source is absolutely continuous. The matching positive future representation ∫_{f≥t} exp(zx−z²(f(x)−t)/2)dγ(x) is therefore ≤C e^(cz)/z² for z≥1, uniformly on the compact time band (use the bounded derivative of p_t on the whole line). Taking z²=2/h gives γ{x≥c:t≤f(x)≤t+h}≤Ch. Monotonicity of the high wing excludes a plateau at t and makes R_1 locally Lipschitz across the whole band. For the canonical curve the true stopping argument yields the same Lipschitz bound. All statements about pointwise densities below are local to the positive edge. Lower-branch plateau atoms in time do not affect an off-diagonal spatial kernel at the right edge.

As proved in the slope note, subtract the right exterior value at R_i(t)+u from the interior value at R_i(t)−u. The same-wing kernel has a difference quotient dominated by C(t−s)^(−1/2) near s=t by Lipschitzness, and on the opposite wing the arguments are separated. For the latter wing use its finite stopped time-space measure μ_i^−(ds,dx) with x≤0; no time density there is needed. Its integrand p'_r(d) extends continuously by zero at r≤0 for d separated from zero, uniformly on compact times. This gives a continuous interior slope j_i(t), whether or not Q_t is a positive measure on its whole support. The signed negative far part of Q_t contributes negligibly to the right exponential transform: off any fixed right-edge neighborhood its total variation is bounded and its support is a fixed distance to the left. The future Tauberian argument in the slope note consequently applies to the actual higher wing on the entire band, with its strict remote margin. Even if j_1(t)=0, the scaled positive Laplace tail gives the zero right derivative (bound its sublevel measure by e times the transform). For the canonical curve there is no need for a lower variance bound on the remote tail: at any slightly later time t+η the future transform is bounded by exp[−z²η/2+R_2(t+η)z], because the survivors at t+η are supported in that interval. The optional exponential identity used here follows from UI of B_(τ∧s) and integrability of exp(pz B_τ) for p>1, by conditional Jensen. Thus both radii are C¹ and, writing u_i=R_i', obey

    g(R_i(t))u_i(t) = −p'_t(R_i(t))
       +∫₀ᵗ p'_(t−s)(R_i(t)−R_i(s)) g(R_i(s))u_i(s)ds
       +∫_{s≤t} p'_(t−s)(R_i(t)−x) μ_i^−(ds,dx).                   (2)

The first integral before the local Lipschitz band means integration against g(R_i(s))dR_i(s), a common finite measure. The slope factor on both terms is j_i/2. In particular j_1(T)=j_2(T)>0, and both stay positive for a small further time.

**Comparison in (2).** For small h>0 set E_h=sup_{T≤s≤T+h}|u_1(s)−u_2(s)|. The two radii differ at time t∈[T,T+h] by at most (t−T) E_h. The denominator factor g and the free term give O(h E_h) in the difference equation. The old sources s<T coincide. Early times separated from T cost O(h E_h) when moving the observation point, since the second-derivative kernel is uniformly bounded there. The old right sources in a compact band ending at T have bounded time density, and |p''_(t−s)|≤C(t−s)^(−3/2); integrating up to T and multiplying by (t−T)E_h costs O(sqrt(h) E_h). The old left sources are separated spatially.

For the new right self-interaction s∈[T,t], writing d_i=R_i(t)−R_i(s) we have |d_1−d_2|≤(t−s)E_h. The change of kernel thus costs C(t−s)^(−1/2)E_h. The change of source strength is O(E_h); the untouched first derivative kernel obeys |p'_(t−s)(d_i)|≤C(t−s)^(−1/2), since |d_i|≤C(t−s). Both integrated losses are O(sqrt(h) E_h).

All new left sources have x≤−b and R_i(t)≥b, so their distance from the observation is at least 2b. Uniformly for 0<r≤h sufficiently small and d≥2b,

    |p'_r(d)| ≤ C h^(−3/2)exp(−2b²/h).

Each new left time-space measure on [T,T+h] has mass O(h): for the canonical true stop use its local Lipschitz bound; for the actual graph use the *symmetric preactivation bound at T*, γ{x≤−b:T≤f(x)≤T+h}≤Ch. Neither has a plateau at T. The new opposite-wing terms are therefore bounded by C h^(−1/2)exp(−2b²/h) in absolute value. The positive continuous factor g is bounded below on the whole fixed comparison band. Taking the supremum leaves E_h≤C sqrt(h) E_h + C h^(−1/2) exp(−2b²/h); absorb the first term. Integrating over h (or substituting v=2b²/h in the endpoint tail) proves (1).

**The low-wing direction.** In contrast, the lower-wing F_− must exceed the canonical inverse F_2 at points in every spatial interval (b,b+η). Otherwise compare both negative-frequency future transforms at T. The common surviving law has the same transform. The actual outward-collar contribution on [b,b+η] is at least the canonical one there if F_−≤F_2. For the canonical tail beyond this collar take d>0 so small that R_2(T+d)<b+η; the transform of every stopping event later than T+d is bounded by exp[−z²d/2+R_2(T+d)z]. The canonical inward-frequency opposite wing costs O(exp(−bz)). Meanwhile the actual outward interval [z−1,z] costs at least a positive constant under the tilted Gaussian integral, since F_−≤1+T by monotonicity. This is a contradiction. These results give no true-stopping assertion for the actual low wing past T.

A further polynomial factor follows by retaining the cumulative form of the new opposite-wing mass bound. For 0≤a≤h each such new measure satisfies M([T,T+a])≤C a. For r sufficiently small k(r)=r^(−3/2)exp(−2b²/r) is increasing and k(0)=0. Stieltjes integration by parts with the decreasing function s↦k(h−s) yields

    ∫_[0,h] k(h−s)dM(T+s) ≤ C∫₀ʰ k(r)dr
                               = O(h^(1/2)exp(−2b²/h)).

This bounds the new source term uniformly also for intermediate t≤T+h. The same absorption and integration sharpen (1) to

    |R_1'(T+h)−R_2'(T+h)| ≤ C h^(1/2)exp(−2b²/h),
    |R_1(T+h)−R_2(T+h)| ≤ C h^(5/2)exp(−2b²/h).
