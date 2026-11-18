# Heat return to a cluster of terminal maxima

Let X have standard normal density g. Assume 0≤f≤H=sup f, v=H−f, and the exact Gaussian output law for X+√v(X)G, with G an independent standard normal. Let Z={v=0}. For every s>0 the smoothed output density is p_{H+s}=∫p_{v(x)+s}(·−x)g(x)dx. The results below are local near a chosen compact nonempty A⊂Z and use no assumptions at infinity. Assume f continuous for the terminal-set interpretation. The basic estimate also holds for measurable v with the stated output law.

## Positive measure tests and an irregular cluster

For every probability μ supported on A,

    ∫_{v(x)>0} (p_{v(x)}*μ)(x) g(x)dx ≤ ∫ p_H(a) μ(da).                   (1)

Indeed at fixed a apply Fatou on {v>0}, then integrate in a by Tonelli. This bound does not depend on a regularity property of μ. More precisely its right side can be decreased by ∫ ζ(a) μ(da), where ζ(a)=limsup_{s↓0} P_s(g 1_Z)(a). At each fixed a take a sequence realizing that limsup and apply Fatou to the remaining positive-variance contributions. The limsup is Borel by continuity in s>0 and restriction to rational s. The zero contribution can, for example, exhaust g(a) at a Lebesgue density point of Z. If |Z|=0 it vanishes. R021, proofs/021-contact-phase.md lines 3–12, gives a precise pointwise return identity by also accounting for positive variances below the parabolic scale.

Write r(x)=dist(x,A). Let h_-(u)=inf_{a∈A} μ([a−u,a+u]). For each k≥0 and x with v(x)>0 and r(x)≤k√v(x), take a nearest point of A. The points within √v(x) of that point are all within (k+1)√v(x) of x. Consequently

    (p_{v(x)}*μ)(x) ≥ e^{−(k+1)²/2} h_-(√v(x))/√(2π v(x)).              (2)

One can use the actual μ-mass about a nearest point instead of its infimum. Equivalently, without a nearest-point choice, the integrand in (1) is exactly

    (2πv)^(−1/2) ∫₀^∞ t e^(−t²/2) μ([x−t√v,x+t√v])dt.

This identity follows by integrating t e^(−t²/2) from |a−x|/√v to infinity. Thus (1) retains the mass variations of an irregular cluster even when its uniform lower mass bound is uninformative.

## Ahlfors-regular clusters and partial shells

Suppose 0<d<1 and c u^d≤μ([a−u,a+u])≤C u^d for every a∈A and all sufficiently small u. With any fixed k>0, (1)–(2) imply

    ∫_{0<v(x)<v₀, r(x)≤k√v(x)} v(x)^{−(1−d)/2} g(x)dx < ∞.               (3)

The numerical constant is bounded in terms of k, the lower regularity constant and the exact right side of (1). The upper regularity constant is only needed in the volume calculation below.

Let A_t={x:r(x)<t}. A maximal t-separated subset covers A with t-balls, and its disjoint balls of radius t/3 have μ-mass at least c' t^d apiece. Thus it has at most C't^(−d) and, by the upper mass bound and covering, at least c't^(−d) points. Applying the same disjoint and covering intervals to Lebesgue length gives c₁t^(1−d)≤|A_t|≤C₁t^(1−d). Choose fixed M>1 large enough that C₁M^{−(1−d)}<c₁/2. The shells S_n=A_{M^(−n)} minus A_{M^(−n−1)} then have length comparable to M^{−n(1−d)}. In particular |A|=0.

Set t_n=M^(−n), B_n≥1 with t_n²B_n→0. Suppose a measurable subset E_n⊂S_n of relative Lebesgue length θ_n satisfies

    0<c t_n² B_n≤v(x)≤C t_n² B_n    (x∈E_n),

with constants independent of n. Then r(x)≤k√v on these subsets for fixed k, g is bounded below near A, and (3) gives the scale budget

    Σ_n θ_n B_n^{−(1−d)/2}<∞.                                                (4)

Additional maxima and other variance profiles throughout the complement of the E_n are allowed. Assuming the profile throughout a full punctured neighborhood, up to null sets, allows θ_n=1. Ascribing full-shell estimates to a distance profile that holds only on the positive-variance portion requires a separate lower bound on that portion's relative measure. Taking A strictly smaller than Z is harmless for (1)–(4) with the E_n as specified.

## Logarithms and the critical coefficient

If v is comparable on a fixed positive fraction of each shell to r²L^β(log L)^γ, where L=log(1/r) and B=L^β(log L)^γ stays bounded below away from zero, (4) requires

    β>2/(1−d), or β=2/(1−d) and γ>2/(1−d).

In particular the doubly critical choice of β and γ is excluded. For a further factor (log log L)^δ at this doubly critical choice one must have δ>2/(1−d). These statements follow from the ordinary logarithmic series criterion. If each θ_n may decrease then (4) gives its required summability with the indicated weight. For example at the single critical power B_n≈n^{2/(1−d)}, the active relative lengths must satisfy Σ θ_n/n<∞.

## Non-power mass gauges

More generally suppose c h(u)≤μ([a−u,a+u])≤C h(u), with h(u)=u^d L₀(1/u) a positive regularly varying mass gauge of index 0<d<1 and L₀ slowly varying. The preceding packing calculation gives |A_t| comparable to t/h(t). For sufficiently large fixed M the shell lengths are comparable to t_n/h(t_n), since [M^(−1)t/h(t/M)]/[t/h(t)]→M^(d−1). For the same good subsets, (2) and regular variation under fixed multiplicative comparability give

    Σ_n θ_n h(t_n√B_n) / [√B_n h(t_n)] < ∞,                                 (5)

or equivalently Σ θ_n B_n^{−(1−d)/2} L₀(1/(t_n√B_n))/L₀(1/t_n)<∞. Here t_n√B_n→0. The ratio of slowly varying factors should be kept for a general gauge, because the dilation √B_n may be unbounded. For L₀(T) that is a finite product of powers of iterated logarithms and log B_n=o(log(1/t_n)), that ratio tends to one. For an arbitrary compact set with only the lower mass function h_-, apply (2) directly with actual shell volumes.

These inequalities are necessary. Convergence of the return-energy series by itself gives no construction of an exact graph.
