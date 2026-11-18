# Repeated loss and recovery of heat positivity for a single deterministic graph

Choose in advance strictly interlaced times s₁<t₁<s₂<t₂<⋯ increasing to s_*<1, all in [7/8,15/16]. There is a bounded deterministic exact graph, continuous off one compact nowhere-dense subset of [2,2.1], with a strictly negative smooth heat test at a time in (s_j,t_j) for every j, yet

    Q_t≥(1/2)p_t dx            whenever t_j≤t<s_{j+1}, and at s_*≤t<1.

The negative tests may depend on j. The early time distribution need not contain a single atom. After time 1 the remainder is nonnegative by the actual residual Brownian embedding. Thus positivity at arbitrarily late stages of the early history does not propagate through the next onset. The construction does not establish this behavior for continuous graphs.

Start with the thick-complement Cantor set K₀ used in round 154. Inside it choose disjoint compact positive-measure cylinder sets K_j, ordered so they accumulate only at one endpoint a of K₀. They may be chosen arbitrarily small individually: for example reserve one distinct right daughter branching off each node of the leftmost infinite path, then select a sufficiently deep cylinder inside that daughter. Put K={a}∪⋃K_j. As a closed subset of K₀ it inherits the complement thickness lower bound. Set η_j=1_{K_j}γ with masses m_j. For s∈[7/8,15/16], x∈[2,2.1], and s<t≤1, completing the Gaussian square gives

    sup_y p_{t−s}(y−x)/p_t(y)=√(t/(t−s)) exp(x²/(2s))≤C/√(t−s).

Choose m_j so C√2 m_j/√(t_j−s_j)≤2^(−j−2); their total can also be arbitrarily small.

At the hypothetical common time s_j the local slice p_{s_j}dx−η_j is strictly negative on K_j, since p_{s_j}/g<1 there. Select a bounded smooth nonnegative test supported near K_j with negative integral against this signed measure. Choose δ_j>0 so 2δ_j<(t_j−s_j)/2, and assign an affine nonconstant time h_j(x)∈[s_j,s_j+δ_j] on K_j. At t=s_j+2δ_j all of this source has been subtracted. Its pairing against the chosen smooth test tends to its pairing at s_j for the frozen time as δ_j↓0: heat ages are at most 2δ_j and the test has bounded second derivative. The free test integral also converges. Choose δ_j small enough to retain strict negativity. Earlier sources subtract additional nonnegative mass and no later ones have yet acted. Since each h_j is injective on its cylinder, the early time marginal has no atoms. The single point a may be given the value s_*.

On [t_j,s_{j+1}) only the first j cylinders contribute, and their combined propagated density is at most p_t Σ_{i≤j}2^(−i−2) by the square-completion bound. At and above s_* the same argument sums all cylinders. At time 1 put μ=γ−Σ∫_{K_j}p_{1−h_j(x)}(·−x)γ(dx), λ=1_{Kᶜ}γ. Then μ≥γ/2 and λ−μ is the difference between the propagated early submass and its original spatial law. These measures have the same mass and mean and are in convex order by conditional Gaussian smoothing. The occupation potential for their Root completion is

    U(y)=Σ∫_{K_j}∫₀^{1−h_j(x)} p_v(y−x)dv γ(dx).

It obeys the same upper bound as round 154 with ε=1/8 and m=Σm_j. The exterior trap and the compact component-halving argument of that round give a bounded barrier r: in the compact region the thickness of the fixed K₀ is enough even if the cylinders have very unequal sizes. Local continuity holds off K. The full graph equals h_j on K_j and 1+r on the complement. Its exact real transform follows by bounded optional sampling and the variable-duration Gaussian identity, cylinder by cylinder. Any endpoint value on the null singleton leaves the law unaffected. This verifies all intervals for one fixed graph and one strictly positive terminal time bound.
