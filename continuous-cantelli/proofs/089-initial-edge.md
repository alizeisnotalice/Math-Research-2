# The universal initial Stefan edge and the large-tilt saddle

Assume the genuine early survival representation in round 84 on a time interval (0,t₀), with initial δ₀ and deposited spatial density g=p₁. By enthalpy uniqueness the interval is symmetric; denote its radius by b(t). Put L=log(1/t). Then

    b(t)²/t = L + 2/L + 1/L² + o(L⁻²).                               (1)

**Integrated equation.** For U=∫₀ᵗ q_s ds, the identity for the past and surviving mass gives U_t−U_xx/2=δ₀−g 1_{|x|<b(t)} in distributions. U is positive exactly in the interval interior: Brownian paths can remain inside and reach compact interior neighborhoods with positive density; this also follows from interior heat positivity for the nonzero survivor once that interval is open. Boundary points have zero time integral. U has zero initial data. Local spatial derivatives in the distributional equation suffice throughout.

**Caloric barriers.** Fix any real c. For small r>0 let M(r)>1 be the large root of r²=M exp(−M). Extend k_c(r)=g(0)[1−1/M+c/M²] to a bounded positive even function, bounded away from zero and smooth off zero. Set w(t,x)=P_t k_c(x) and q_c=(p_t−w)_+. Its positive region for small t is an interval about 0: w/p_t is a strictly increasing function for x>0 by pairing ±y in the heat integral and using the strict increase of cosh(xy/t); it starts below 1 and eventually exceeds 1. At its unique endpoint b_c(t) put s=b_c²/t. Since w is bounded above and below, s=L+O(1).

We give details of the needed uniform local expansions. At x=b_c the local window |y−b_c|≤b_c/2 has M(y)=L+O(1) by its defining equation and s=L+O(1). Derivatives there satisfy |k'_c(y)|≤C/(b_c L²), |k''_c(y)|≤C/(b_c² L²) by implicit differentiation. Taylor's formula about b_c in the symmetric window yields w(t,b_c)=k_c(b_c)+O(t/(b_c² L²))+O(exp(−c₁L)). For derivatives of w, integrate by parts with a smooth cutoff identically one in the half-size window and zero outside the window; derivatives of the cutoff cost at most C_j b_c^{-j}, and the Gaussian tail and its derivatives outside the inner window are exponentially small in L after multiplication by t^{j/2}. Thus |w_x|≤C/(b_c L²), |w_xx|≤C/(b_c² L²). These estimates avoid differentiating k across its origin singularity. The root gives s=L−2log(w/g(0))=L+O(1/L), hence M(b_c)=L+O(L⁻²) from M−log M=L−log s. It follows that w/g(0)=1−1/L+c/L²+O(L⁻³), and

    s=L+2/L+(1−2c)/L²+O(L⁻³).

Implicit heat differentiation gives b_c'=−(p_xx−w_xx)/(2(p_x−w_x))>0 for small t by the estimates above. The boundary density in the integrated equation of q_c, in spatial units, is therefore

    g_c(b_c(t))=(p_x−w_x)²/(p_xx−w_xx)
                =w(t,b_c) s/(s−1) [1+O(L⁻³)]
                =g(0)[1+c/L²+O(L⁻³)].

Indeed the heat operator applied to the clipped caloric function q_c has at each boundary the mass (p_x−w_x)/2 on the right and its reflection on the left. Its initial measure is δ₀ (b_c/√t→∞, and w is bounded). Time integration thus gives (U_c)_t−(U_c)_xx/2=δ₀−g_c 1_{U_c>0}, for U_c=∫₀ᵗ q_c. This defines g_c locally by the inverse of the strictly increasing boundary.

Since g(b_c)=g(0)[1+O(tL)], c>0 gives g_c≥g on all sufficiently early swept positions. The comparison U_c≤U follows by testing their difference equation against (U_c−U)_+: wherever that test is positive its forcing is nonpositive, as both functions are nonnegative. Singular initial sources cancel; the remaining forcing is bounded with compact spatial support on a sufficiently short time slab. Hence the elementary heat energy test is valid by approximation, with zero initial difference. Reversing inequalities gives U≤U_c for c<0. Support inclusion sandwiches b(t) between the two radii. First let t↓0 with c fixed, then let positive and negative c approach 0. This proves (1).

**Saturation of the global entropy bound.** Retain the hypotheses of round 84, including liminf f>1 on both tails and its unique zero at zero. Write t hereafter for a large real tilt, u for the positive solution of u+log u=2log|t|, and Z=tX under the tilted law. The tilt has density g(z/t) exp(z−log|t|−t²f(z/t)/2) in the z coordinate. Then

    A(t)=t²E_t f=E_t Z=u−1+o(1).

In particular A+log A=2log|t|−1+o(1). This is a conditional necessity for exact graphs in the stated early regime.

Here are error details for the saddle calculation. By early symmetry f is even near zero. Write T_z=t²f(z/t), s_z=z²/T_z; then (1) gives s_z=log(1/f(z/t))+O(1/log(1/f(z/t))) at small nonzero |z/t|. On z∈[u/4,2u], first s_z=u+O(1): indeed s_z−log s_z=2log|t|−2log z+O(1/log(1/f)), with large s_z, and the right hand side equals u−log u+O(1). If z=u+√u w, |w|≤u^(1/20), then log T_z=log u+2w/√u+O((1+w²)/u). Consequently s_z=u−2w/√u+O((1+w²)/u), and T_z=u+2√u w+w²+2w/√u+O((1+|w|⁴)/u). After the Jacobian √u the density in w equals g(w)[1−w/√u+O((1+|w|⁸)/u)] uniformly on this window. On the rest of [u/4,2u], the preliminary bound s_z=u+O(1) gives a constant times exp(−w²/2) as a dominating density. On positive z≤u/4 use the direct bound g(0)exp(u/4−log|t|); for all negative z the factor exp(z−log|t|) suffices. On z>2u in a fixed small origin neighborhood use monotonicity there to bound s_z≤log(1/f(2u/t))+O(1)≤u+C, again giving an integrable Gaussian tail centered near u. Compact physical annuli give exponentially small mass and first absolute moment; on the favorable far tail use f≥1+δ and complete the Gaussian square, which beats all powers of |t|. These estimates give E[(Z−u)/√u]=−1/√u+O(1/u) by integration, proving the assertion.
