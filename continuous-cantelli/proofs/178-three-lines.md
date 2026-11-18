# The nullspace threshold for two and three affine source segments

Write Cμ(z)=∫exp(zx−z²t/2)dμ(t,x). All source measures here are finite, signed, and compactly supported in real spacetime. Their transforms are entire. Three distinct affine worldlines x=a_i+v_i t, present on a common nondegenerate time interval, support a nonzero real smooth-in-time null source on any smaller open subinterval on which the lines are pairwise disjoint. Two distinct affine worldlines support no nonzero signed spacetime null measure, whether parallel or intersecting.

## Two lines

Write H_i(w)=∫e^{wt}dh_i(t) for the temporal transforms, absorbing coordinate signs in h_i. For distinct parallel lines with common slope v, set w=vz−z²/2. Evaluating the null identity at z and 2v−z leaves a two-by-two system with nonzero determinant on an open real interval, because the offsets differ. Hence both entire H_i vanish.

For nonparallel lines let (t₀,x₀) be their intersection and change the sign of one time measure, so that the two individual spacetime transforms coincide. Divide both by exp(zx₀−z²t₀/2). The common entire G(z) can now be expressed as exp(−t₀w_i)H_i(w_i) for either i, w_i=v_i z−z²/2. Hence G is invariant under both real reflections z↦2v_i−z and has nonzero real period 2(v₁−v₂). It is bounded on the real axis. For either representation w_i ranges over every sufficiently negative real value. The compact signed Laplace support lemma of R083 implies h_i supported on [t₀,∞) (apply it after translating the left endpoint of any initial containing compact interval). Thus the integral for G tends to h_i({t₀}) as z→+∞ by dominated convergence. A periodic continuous function admitting such a limit is constant. Uniqueness of the time Laplace transforms forces both h_i to consist of the same atom at t₀, which cancels at their common spatial point. This proves injectivity for spacetime measures and also covers disjoint or truncated segments.

## A moving interior point in a fixed interval

First take outer lines x=0,L and middle line X(t)=A+Bt on a compact time interval [s₀,s₁], with a positive distance to both walls. Choose any nonzero real compactly supported Gevrey bump ψ in (s₀,s₁) of order 1<r<2; thus ∥ψ^(m)∥∞≤C R^m(m!)^r. Set D=∂_t and define

    h=[cosh(2LB)−cosh(2L√(B²+2D))]ψ,

where cosh(c√p)=Σ_{n≥0}c^{2n}p^n/(2n)!. For each fixed j, applying D^j to the series gives convergence in the uniform norm since the nth differentiated term is at most C_j A^n ((n+j)!)^r/(2n)!. This follows directly from the binomial expansion of (B²+2D)^n and the bump bound. The support of h is contained in supp ψ and h is real and smooth. Integration by parts gives

    H(w)=∫e^{wt}h(t)dt
        =P_B(w)∫e^{wt}ψ(t)dt,
    P_B(w)=cosh(2LB)−cosh(2L√(B²−2w)).

Neither entire factor is identically zero, so h≠0. At k_n=nπ/L, λ_n=k_n²/2 and w=λ_n±iBk_n, one has B²−2w=(B∓ik_n)², whence P_B(w)=0 for every n≥1, also when B=0.

Let u solve the Dirichlet heat equation (∂_t−½∂²_x)u=h(t)δ_{X(t)} in (0,L), initially zero at s₀. It is explicitly the time convolution of h with the Dirichlet heat kernel and the moving source position. At s₁ its nth sine coefficient, up to a fixed normalization, is

    exp(−λ_n s₁) ∫ exp(λ_n t) sin(k_n(A+Bt))h(t)dt=0

by the two zeros of H. Dirichlet completeness gives u(s₁,·)=0. The same observation at any time after the final point of supp h (and before s₁) gives zero there as well; likewise u=0 before supp h. Define u to be zero outside the closed spatial interval and outside [s₀,s₁]. There are no time boundary measures since the initial and terminal data both vanish. The interior solution has zero Dirichlet value, hence no spatial derivatives of point masses occur at the walls. Its full heat source is exactly

    μ=h(t)dt δ_{X(t)} − ½u_x(t,0+)dt δ₀ + ½u_x(t,L−)dt δ_L.          (1)

The boundary fluxes are smooth and integrable. This follows from the reflected Gaussian series for the Dirichlet heat kernel: the middle segment stays a fixed positive distance from the walls, so its boundary normal derivative and every fixed derivative in the age extend by zero at age zero and are bounded on the finite time interval uniformly over the source segment. Both fluxes are supported in the closed time hull of supp h (gaps within supp h's hull need not persist at a wall). Formula (1) is a compact signed measure and is nonzero by its interior component h. Testing it against exp(zx−z²t/2) gives zero by distributional integration by parts against the compactly supported u, since the test solves the backward heat equation.

## Projective and Galilean changes of coordinates

Time and space translations preserve the zero-transform property by multiplication by a nonzero exponential scalar. Under a Galilean map (t,x)↦(t,x−dt), use the positive weight exp(dx−d²t/2) before pushing forward the measure. The new test transform is the old transform at z+d. This transformation is invertible by the parameter −d, whose pullback weight is exactly the reciprocal of the original weight.

For a real c and compact support where 1+ct>0 define

    T_c(t,x)=(t/(1+ct), x/(1+ct)),
    G_c(t,x)=(1+ct)^(-1/2) exp[cx²/(2(1+ct))].

The measure (T_c)_*(G_c μ) is null whenever μ is. For c≥0 take V=z+√c N with N standard normal and z real. The Gaussian expectation of exp(Vx−V²t/2) equals G_c times the new exponential test at T_c(t,x), by completing the square. For c<0 take V=z+i√(−c)N instead; the same Gaussian integral converges absolutely and gives the same identity. On the compact window 1+ct is bounded strictly away from zero, which controls Fubini in both cases. The original transform vanishes for complex V too. The inverse is T_{−c} and G_{−c}(T_c(t,x))=1/G_c(t,x) exactly. Thus the equivalence acts in either direction, never produces zero from a nonzero source, and preserves smoothness and finiteness of line densities on compact windows. No time Jacobian is inserted in G_c; the ordinary Jacobian for a new time-density is already part of the measure pushforward.

Translate an interior time of the common subinterval to zero, and index the two outer lines by their spatial order so that a₂>a₁. Shrink the time window to avoid all intersections if necessary. Take c=(v₂−v₁)/(a₂−a₁). Under T_c, a line x=a_i+v_i t becomes x'=a_i+(v_i−ca_i)t', so the two outer lines become parallel. Moreover 1+ct is their original separation divided by a₂−a₁ and is strictly positive throughout the chosen interval. After shifting space and subtracting their common slope by the Galilean map, the outer lines are 0,L with L=a₂−a₁. Time ordering and the strict ordering of the three spatial branches at each time are preserved. Construct (1) on the transformed window and pull it back by the respective inverses. This proves the three-line statement for arbitrary offsets and slopes, without a concurrency hypothesis.

If one wants the three segments to lie on a single continuous graph t=f(x), one may take nonzero slopes and shrink the time window further until their spatial projections are disjoint, then interpolate across the unused gaps. This is a signed nullspace statement for that graph; it does not prescribe a Gaussian spatial marginal.
