# A common last active time and simultaneous signs on every interior interval

Assume R161's continuous bounded exact graph and let C be any connected component of {f<1}. On (1,H]×C write q(t,x) for the density of the signed slice. There is a time σ=σ(C)∈(1,H] such that q vanishes on [σ,H]×C. For every nonempty open subinterval J of C and every 1≤s<σ there exist a single t∈(s,σ) and x₁,x₂∈J such that

    q(t,x₁)>0 and q(t,x₂)<0.

In particular both signs, and times of actual spatial sign change, have the same supremum σ in every J.

**Analyticity in space.** Fix x₀∈C. Choose d>0 and τ<1 with f≤τ on [x₀−2d,x₀+2d]. On the complex disk |w−x₀|<d/2, the positive-age Gaussian kernels from this segment are holomorphic in w and locally bounded uniformly when 1≤t≤H, since t−f≥1−τ. For exterior source positions x', |Re w−x'|≥3d/2 while |Im w|<d/2. Consequently the real part of (w−x')² is bounded below by 2d². Each kernel of positive age v is bounded uniformly by (2πv)^(−1/2)exp(−d²/v) there, which remains bounded down to v=0. Sources of exactly zero age are spatially remote and have no local measure; negative ages are omitted. Dominated holomorphic integration proves spatial real analyticity of q(t,·) on C for each t∈[1,H]. The same local separation gives ordinary spacetime smoothness around these real points. No time analyticity is used.

**Extinction and the common time.** By R161 any finite interval compactly inside C has at least one endpoint with strictly negative heat density at some time after 1. Indeed it admits a neighborhood below a fixed τ<1. Thus the set of times t∈(1,H) with q(t,·) nonzero on C is nonempty; let σ be its supremum. Continuity and Q_H=0 imply q(σ,·)=0 and all later local slices vanish. Analyticity shows that on every J the nonzero-time set is exactly the same as on C.

**Simultaneous signs.** Fix J and s. If there are no mixed times on (s,σ), every slice there has a weak spatial sign on J. At least one slice is nonzero on J. Replacing q by −q if needed, there is t₀∈(s,σ) such that some point of J has a strictly positive value. The set of times with a positive value in J is open. Let (u,v) be its component relative to (s,σ) that contains t₀. Then q≥0 throughout (u,v)×J. At v it is identically zero on J, by continuity and absence of positive values at a finite endpoint v<σ, or by definition when v=σ. Take a bounded interval J' compactly inside J containing a positive point at t₀. On its two endpoints q stays nonnegative through time v. Comparison with the Dirichlet heat flow started at t₀ from a nonzero nonnegative smooth bump dominated by q(t₀,·) in J' shows q(v,·)>0 throughout J'. The strict positivity of the killed heat kernel for all positive durations gives this even at the final time v. The contradiction proves the theorem.

## Arbitrarily many spatial changes of sign for exponentially weighted late slices

There is a stronger assertion for time averages. On every bounded open interval J compactly contained in C and every 1≤s<σ, the functions

    F_λ(x)=∫_s^σ exp(λ(t−s))q(t,x)dt

have unbounded numbers of strict sign alternations in J along some sequence λ→∞. In fact there is a sequence along which there are at least c_J√λ alternations on any fixed smaller nonempty interval, where c_J can be replaced by a positive constant depending on the smaller interval.

Choose x₀ in J with q(t₀,x₀)≠0 for some t₀∈(s,σ). The signed time measure q(t,x₀)dt on [s,σ] has upper support endpoint d>s. The signed compact Laplace-support lemma of R083 §3, applied after time reversal, yields limsup_{λ→∞} λ⁻¹ log|F_λ(x₀)|=d−s>0. In particular there is a sequence with exponentially large absolute value. For spatial points in a compact interval containing x₀ and any specified smaller interval in J, integration of the heat equation gives

    F_λ''+2λ F_λ=−2q(s,·).

Writing k=√(2λ), the solution is A cos(k(x−x₀))+B sin(k(x−x₀)) plus the forced convolution −(2/k)∫_{x₀}^x sin(k(x−u))q(s,u)du. The latter is bounded by a fixed constant divided by k, whereas the amplitude of the homogeneous wave is at least |A|=|F_λ(x₀)|. Alternating extrema spaced π/k inside the smaller interval give the claimed number of strict alternations. This result concerns weighted time integrals; such an integral alone does not exhibit a single time slice with the same number of sign changes.
