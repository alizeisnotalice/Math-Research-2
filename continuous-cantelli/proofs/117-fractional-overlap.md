# A universal fractional likelihood envelope and its equality criterion

Let γ have standard normal density g. Suppose f≥0 is measurable and R_s(y)=exp(sy−s²f(y)/2) integrates to one under γ for all real s. Fix 0<k<1 and let s→+∞. Define L by s=exp(L)√(2L), a=√(2L)/s=exp(−L), y₀=2L/s. Then

    limsup a^{k−1} ∫R_s^k dγ ≤ 1/√k.

Along any sequence, convergence to the displayed upper bound is equivalent to convergence of the law of (X−y₀)/a under R_s dγ to N(0,1) in total variation. It is also equivalent to the local profile
s² f(y₀+a z)=(√(2L)+z)²+o(1) in Lebesgue measure on each bounded z interval. Thus saturation at any one power implies saturation at every fixed power in (0,1) at the same base parameters. Reflection covers negative s. The bound applies even to a random nonnegative mark f on an enlarged space; the equality assertions here are for a deterministic graph. (For a random mark the upper bound follows by the same two-scale argument on the joint space, or conditioning in the Hölder steps; see clarification below.)

**Tail localization.** Write g₀=g(0). For y>y₀+Ma take u=1+1/√(2L); for y<y₀−Ma take u=1−1/√(2L). Set r=ku²<1 for large s. The identity R_s^k=R_{s/u}^r exp(ks(1−u)y) and Hölder show that on either half-line the integral is at most
exp(ks(1−u)y₀−kM)[g₀(1−r)a/k]^{1−r}.
Dividing by a^{1−k}, the L-dependent exponential at y₀ becomes exp(kL(u−1)²)=exp(k/2). Thus the normalized tails are O_k(exp(−kM)), uniformly for large s and any M≥0.

**Sharp bound on the central window.** In the coordinate z=(y−y₀)/a let p(z)=a g(y) R_s(y), a probability density. The normalized integral is ∫p^k g(y)^{1−k}dz. On |z|≤M discard p outside [exp(−L^{1/4}),exp(L^{1/4})]. Its contribution vanishes at fixed M since ∫p=1. On the retained set G, log R_s=L+O_M(L^{1/4}) and sy=2L+√(2L)z. From
log(R_{s(1+e)}/R_s)=(2e+e²)log R_s−(e+e²)sy
with e=λ/√(2L) one obtains −λz−λ²/2+o(1) uniformly there for bounded λ. Therefore limsup ∫_G exp(−λz)p(z)dz ≤ exp(λ²/2). Integrate in λ against N(0,1−k). At fixed M the integrands are bounded by exp(M|λ|), so reverse Fatou gives limsup ∫_G p exp((1−k)z²/2) ≤ k^(−1/2). Weighted Hölder now gives
∫_G p^k ≤ (∫_G p exp((1−k)z²/2))^k (∫exp(−kz²/2)dz)^(1−k).
Multiplying by g₀^{1−k} and using tail localization then M→∞ proves the bound.

**Equality.** Along any candidate saturating sequence take M→∞ sufficiently slowly by a diagonal choice for all preceding fixed-M estimates. With A=∫_G p exp((1−k)z²/2), saturation forces A→1/√k and equality asymptotically in the weighted Hölder bound. The order-k affinity between p1_G exp((1−k)z²/2)/A and the normalized density proportional to exp(−kz²/2) thus tends to one. This implies total variation convergence: splitting on a set realizing the variation distance and applying Hölder on each part bounds the affinity by the binary affinity, which is uniformly below one on pairs of probabilities a fixed positive distance apart. Multiplication by the bounded factor A exp(−(1−k)z²/2) shows p1_G→g in L¹. The remaining mass tends to zero. Conversely p→g in L¹ gives local convergence of the fractional-power integrals by |p^k−g^k|≤|p−g|^k; uniform tail localization proves saturation. Finally the exact density identity is
s² f(y₀+a z)=2L+2√(2L)z+2log g(y₀+a z)−2log p(z).
On compact z sets, 2log g(y₀+a z)→2log g₀, while 2log g(z)=2log g₀−z². This gives the local profile. Conversely this profile implies local convergence in measure p→g; any subsequence has a further subsequence with almost everywhere convergence, and Scheffé's lemma for probability densities gives L¹ convergence.

We claim the deterministic-graph theorem above; for the optional random-mark upper bound replace p(z) by the joint likelihood a g(y)R_s(y,F) against dz times the conditional mark law at y. The high/low truncation bounds are identical since the conditional mark measure has mass one; the same normalization and Hölder argument applies without a conditional density in the mark coordinate. No existence claim for a continuous nonconstant saturator is made.
