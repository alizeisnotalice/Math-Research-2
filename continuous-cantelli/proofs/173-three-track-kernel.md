# Two affine tracks are identifiable, three parallel tracks have a smooth null measure

Parameterize the space-time lines by (t,x)=(t,a_j+vt). For a compact time measure η_j write H_j(w)=∫exp(wt)dη_j(t). Its contribution to the stopped Gaussian transform is exp(a_j z)H_j(w) with w=vz−z²/2. For two distinct offsets the transforms of a vanishing sum evaluated at z and 2v−z give an invertible 2 by 2 system off a discrete set. Both H_j vanish by analytic continuation and Fourier uniqueness.

For three offsets 0,L,2L, L>0, take a nonzero smooth compactly supported bump h of Gevrey order 1<s<2 in any chosen open time interval. For example a translated and rescaled exp(−1/(1−t²)²) on |t|<1, extended by zero, has order 3/2: on a complex disk of radius c(1−|t|) at each real t its modulus is ≤exp(−c'/(1−|t|)²); Cauchy differentiation and maximization give C A^n(n!)^(3/2). Put D=d/dt and define three real signed densities
$$
 \eta_0=e^{2Lv}h,\qquad
 \eta_1=-2e^{Lv}\sum_{n=0}^\infty\frac{L^{2n}}{(2n)!}(v^2+2D)^nh,\qquad
 \eta_2=h .
$$
The series and each fixed derivative converge uniformly: the binomial expansion is bounded by C_j A_j^n((n+j)!)^s before division by (2n)!, for derivative order j. All terms vanish outside the support of h. The transform of D is multiplication by −w. The transformed middle coefficient is consequently −2exp(Lv)cosh(L(z−v))H(w). With Y=exp(Lz) the entire sum vanishes by
$$
 Y^2-2e^{Lv}\cosh(L(z-v))Y+e^{2Lv}=0.
$$
Conversely any compactly supported null triple on these tracks has the two outer transforms in the ratio exp(2Lv): evaluate at both roots as above, subtract the quadratic polynomials and use their common two roots. This also gives the stated middle transform. With v≠0 and sufficiently short time window the spatial projections of the tracks are disjoint, so the three pieces lie on a single continuous piecewise affine graph. They define a nonzero signed measure there with zero transform at every complex z, and give infinitely many independent such measures by using disjoint time windows. Translation of all times or all spatial coordinates multiplies the transform by a common entire factor and preserves the cancellation.

For completeness two nonparallel affine tracks are also identifiable as spacetime measures. Denote their intersection by (t*,x*) and their distinct speeds by v_i. After multiplication by exp(−zx*+t*z²/2), a common transform of the two tracks has the form G(z)=∫ exp((t−t*)w_i(z))dη_i(t) for each i (absorbing the sign into one measure). Each representation is invariant under reflection in v_i, so G is periodic on the real line and bounded there. Each w_i covers all sufficiently negative reals; boundedness implies each η_i is supported on [t*,∞), by the signed moment support lemma of round 163 applied to exp(−(t−t*)). Dominated convergence along real z→∞ gives a finite limit for G, so periodicity makes it constant. The two measures can then only be matching point masses at the intersection, which cancel as spacetime measures.

This exhibits failure of injectivity for a *fixed continuous graph* and even for smooth compactly supported perturbations of spatial measures on that graph. Prescribing the full Gaussian spatial marginal is an additional constraint.
