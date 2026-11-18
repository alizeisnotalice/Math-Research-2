# Spatial oscillation forced by exponentially weighting late times

Suppose the real transform is exact for a bounded continuous graph 0≤f≤H with Gaussian spatial marginal. Fix any nonempty compact interval J of positive length in {f<1}. Take τ<1 with f<τ on a neighborhood of J. On this cylinder use the smooth signed heat density q and Q_H=0 as in round 161. Choose x₀ in the interior of J and t₀>1 with q(t₀,x₀)<0, possible by that round. Put
$$
 F_n(x)=\int_\tau^H e^{n(t-\tau)}q(t,x)\,dt,\qquad k_n=\sqrt{2n}.
$$
There exist a number d>0 and infinitely many positive integers n for which F_n alternates in sign at at least floor(k_n |J|/π) ordered points in J, with absolute value ≥ exp(d n) at each point. One can choose the points on a single equally spaced lattice of spacing π/k_n. Weighting the entire time interval [0,H] gives the same signs and an exponentially large margin on a possibly smaller d scale when tested as measures against nonnegative smooth unit-integral tests in short intervals around these extrema. All the statements concern a common time weight for the entire spatial interval; the sign changes need not occur on a single time slice.

**The scalar growth lemma.** For a finite real signed measure ν on a compact interval in (0,∞), with b=sup supp |ν|, its moments m_n obey limsup |m_n|^(1/n)=b (for ν≠0). The upper bound follows by total variation. For the lower, if the limsup were <c<b, the function Σ m_n/z^(n+1) would be holomorphic on |z|>c and real on the positive real axis there. It agrees with the Cauchy integral ∫(z−u)⁻¹ν(du) first for |z|>b and then in the upper half-plane outside the disk by uniqueness. The negative imaginary part divided by π is convolution with the Poisson kernel at height Im z. Integrating against a smooth test supported in (c,∞) and letting the height decrease to zero gives both zero (analytic real continuation) and ν applied to the test (approximate identity). Contradiction. Apply this to the pushforward of q(t,x₀)dt on [τ,H] under u=exp(t−τ). Since q is nonzero near t₀>1, b>exp(1−τ). In particular some d>0 and infinitely many n have |F_n(x₀)|>2exp(dn).

**Spatial comparison.** Integration by parts in time in the heat equation gives
$$
 F_n''+k_n^2 F_n=-2q(\tau,\cdot).
$$
On J and for x on either side of x₀ the solution differs from
$$
 F_n(x₀)\cos(k_n(x-x₀))+k_n^{-1} F_n'(x₀)\sin(k_n(x-x₀))
$$
by at most 2|J| sup_J |q(τ,·)|/k_n, by the sine kernel for the initial-value problem. The sinusoidal amplitude is at least |F_n(x₀)|. Its extrema on J give the asserted signs. The same reasoning gives stable signs and a lower bound equal to a fixed fraction of the amplitude on intervals of width proportional to 1/k_n around extrema inside J. For the contribution of times [0,τ] the spatial total variation of the weighted measure with weight exp(n(t−τ)) is ≤2τ, since each slice has variation ≤2. Smooth nonnegative unit-integral bumps in the half-amplitude intervals may be chosen with supremum O(k_n), hence this earlier contribution is negligible. At the two ends one may drop the outermost extrema when selecting supports. No pointwise regularity before τ is required.

This supplies global-transform tests even for early components without boundary smoothness, a fixed modulus or a finite order of contact. Their time weights increase toward the terminal time, in contrast to the decaying Laplace weights previously considered. The divergence of spatial oscillation concerns hypothetical exact continuous counterexamples and does not establish their existence.
