# Two separated Lipschitz trajectories determine a signed source

Let X₁,X₂ be real Lipschitz functions on [a,b] with X₂−X₁≥δ>0. Let ν₁,ν₂ be finite signed Borel measures on [a,b]. If for all real z one has

    Σ_j ∫ exp(z X_j(s)−z²s/2) dν_j(s) = 0,

then both measures vanish. Time atoms, including at a and b, are allowed. The assumptions can equivalently be restricted to a common time window containing the compact supports.

Translate a to zero, extend both trajectories constantly to the left and right, and keep both measures zero off [0,b−a]. Fix a larger finite time B>b−a. For the standard heat kernel p define

    U(t,y)=Σ_j ∫_{s<t} p_{t−s}(y−X_j(s)) dν_j(s).

For any y outside the full spatial range of the paths, the function t↦U(t,y) is continuous and has zero Laplace transform. Indeed, for q=z²/2, z>0, the time Laplace transform of p_t(d) equals exp(−z|d|)/z. The sign of y−X_j(s) is constant on either exterior ray, so the source hypothesis gives zero there. All integrals converge absolutely against the total variation when weighted by exp(−qt). Laplace uniqueness gives U(t,y)=0 there for every t≥0.

For fixed t the function U(t,·) is real analytic at points different from both X_j(t). For small ages, the real part of the squared complex distance to the source is uniformly positive in a small complex disk around such a point. Larger ages give bounded holomorphic kernels there. This continues the two exterior zero rays to (−∞,X₁(t)) and (X₂(t),∞). For almost every t in [0,B] one has Σ_j∫_{s<t}(t−s)^(−1/2)d|ν_j|(s)<∞ by Fubini. At these t, dominated convergence yields spatial continuity up to both paths and hence U(t,X_i(t))=0 for i=1,2.

Write A_ij(t,s)=exp(−(X_i(t)−X_j(s))²/[2(t−s)]) for t>s. After multiplying out the common factor (2π)^(−1/2), each trace is an Abel integral with amplitude A_ij. Put e_ij=A_ij−1_{i=j}. For 0≤s<t≤B, e_ij(t,s) extends to 0 at t=s and is Lipschitz in t, uniformly in s. On the diagonal branch, d=X_i(t)−X_i(s) satisfies |d|≤Lip(X_i)(t−s). At a point of differentiability the derivative of the exponent is −d X_i′(t)/(t−s)+d²/[2(t−s)²], which is bounded. On the other branch |d|≥δ/2 for sufficiently small ages by Lipschitz continuity, and the Gaussian suppresses all inverse powers of age. On the remaining compact range the same derivative is bounded. Absolute continuity on intervals away from the diagonal and the continuous extension give the claimed uniform Lipschitz estimate.

Integrate the zero trace for t<T against (T−t)^(−1/2). Fubini holds for every 0<T≤B since the double Abel kernel integrates to π against dt for each s<T. This gives

    πν_i([0,T)) + Σ_j ∫_{s<T} E_ij(T,s)dν_j(s)=0,
    E_ij(T,s)=∫₀¹ e_ij(s+θ(T−s),s) dθ/√(θ(1−θ)).

Extend E to zero for T≤s. The kernels are then uniformly Lipschitz as functions of T. It follows that ν_i([0,T)) is Lipschitz, with value tending to zero at T=0, so that also ν_i({0})=0. Each measure has a Lebesgue density v_i on [0,B), and almost everywhere |v_i(T)|≤C Σ_j|ν_j|([0,T]). This last estimate follows either by differentiating the displayed integral in T almost everywhere using Fubini, or by taking its difference quotient and using the uniform Lipschitz bound and E=0 for T≤s. With V(T)=Σ_j|ν_j|([0,T]), one has V(0)=0 and V′≤2CV almost everywhere. Grönwall yields V=0 throughout [0,B). Padding the time window beyond the original endpoint includes any terminal atoms.
