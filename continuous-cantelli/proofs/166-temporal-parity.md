# Infinitely many temporal sign changes at almost every low cut, with a one-point exceptional set

Assume a bounded, continuous, exact normalized graph f. Let C be a component of {f<1}. By R161 and the common-extinction proof in round 164 there is σ(C)∈(1,H] such that q(σ,·)=0 on C, while every nonempty open subinterval of C has a nonzero slice at times arbitrarily close to σ from the left. On (1,H]×C the density q is smooth, solves the homogeneous heat equation, and is spatially real analytic at every fixed time.

**Theorem.** At most one x in C has q(t,x) of eventual weakly fixed temporal sign on a left neighborhood of σ. At every other x both signs occur arbitrarily close to σ, so by continuity the trace has infinitely many temporal sign alternations with times tending to σ. The possible exceptional trace may be identically zero near σ.

For a<b in C put L=b−a, k_n=nπ/L, λ_n=k_n²/2 and

    M_n(t)=∫_a^b q(t,x)sin(k_n(x−a))dx.

For each 1<s<σ integration by parts twice, followed by integration of the scalar ODE with M_n(σ)=0, gives

    ∫_s^σ exp(λ_n(t−s))[q(t,a)−(−1)^n q(t,b)]dt = −2 M_n(s)/k_n.      (1)

In particular the absolute difference of the two signed endpoint transforms, with plus for odd n and minus for even n, is bounded by 2∥q(s,·)∥_{L¹(a,b)}/k_n. Spatial differentiability at s gives a sharper exact expansion if desired:

    RHS of (1) = −2[q(s,a)−(−1)^n q(s,b)]/k_n² + O_{s,a,b}(k_n⁻⁴).

Indeed integration of M_n(s) once by parts yields its boundary term and k_n⁻¹∫q_x cos(k_n(x−a)); integrating once more replaces the last integral by −k_n⁻¹∫q_{xx}sin(k_n(x−a)). The latter integral is O(k_n⁻¹) by one additional integration by parts. Bounded derivatives on [a,b] suffice.

If both endpoint traces had eventual weakly fixed signs, choose a common s. Choose odd n if their signs agree, even n if they differ; an identically zero trace may be assigned either sign. On this infinite sequence the bracket in (1) has a fixed sign and absolute value |q(t,a)|+|q(t,b)|. Since its integral against a weight at least one is O(1/k_n), both traces vanish throughout (s,σ). Use (1) again for every n and with any starting time r∈(s,σ). Its left side now vanishes, so all Dirichlet Fourier coefficients M_n(r) vanish. Completeness and spatial continuity imply q(r,·)≡0 on (a,b), contradicting its nonzero slices arbitrarily close to σ.

Without any fixed signs, (1) also gives quantitative cancellation for every n. For B_n(t)=q(t,a)−(−1)^n q(t,b), put P_n=∫e^{λ_n(t−s)}(B_n)_+dt and N_n=∫e^{λ_n(t−s)}(B_n)_−dt on [s,σ]. Then |P_n−N_n|≤2∥q(s,·)∥₁/k_n and in fact O(k_n⁻²). On a parity for which the trace combination is nonzero on a closed subinterval of (s,σ), the larger of P_n,N_n grows at least exponentially in λ_n. Thus P_n/N_n→1 exponentially fast up to a polynomial factor along that parity. This estimate is for each of the two endpoint combinations separately, and does not assert that the two raw endpoint transforms have nonvanishing magnitudes.
