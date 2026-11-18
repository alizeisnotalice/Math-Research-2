# A deterministic exact graph with a strictly negative heat remainder

Let ε=1/8, s=1−ε, and let K be the compact nowhere-dense fat Cantor set in round 153 (or a smaller affine copy inside the same J). There exists a bounded nonnegative Borel graph f such that

    ∫exp(zx−z²f(x)/2)γ(dx)=1   (z real),       f=s on K,       f≥1 on ℝ∖K.

It is continuous at every point of the open complement and discontinuous at every point of K. At t=s+10⁻¹² the explicit unscaled K has the same strict negative smooth test as in round 153. In particular its time-position law cannot come from a Brownian stopping time starting at zero, even with independent randomization. For any H>sup f the measurable coefficient √(H−f) gives an exactly Gaussian output of variance H. This last measurable phenomenon was known from Kleptsyn and Kurtzmann; the construction here forces failure of the signed-heat positivity criterion. [Kleptsyn and Kurtzmann](https://arxiv.org/abs/1202.2250).

**Remaining transport.** Put η=1_Kγ, m=η(ℝ), μ=γ−P_εη, λ=γ−η. The Gaussian ratio bound in round 128 uses only η≤1_Jγ and gives μ≥γ/2. Both residual measures have mass 1−m and the same mean after normalization; their potential difference ∫|x−y|(λ−μ)(dy) is strictly positive everywhere by Gaussian smoothing of η. Equivalently μ precedes λ in convex order. Take the regular Root embedding of λ from μ after normalization and common translation. Its stopping time τ is uniformly integrable in the sense of the stopped Brownian martingale, with finite expected duration by the two finite second moments. Write r for the lower semicontinuous barrier and D_t={r>t}. On ℝ∖K the target has a strictly positive continuous density locally, so r is continuous and finite there. Since λ(K)=0 the stop satisfies τ=r(B_τ) almost surely: at any terminal point in the complement this follows by continuity in a neighborhood and continuity of the path. [Root embedding for general starting laws](https://arxiv.org/abs/1505.03169) · [local barrier continuity](https://arxiv.org/abs/2010.14695).

**Uniform time bound.** Use unnormalized laws in all mass computations. Tanaka's formula gives the total occupation density
[
 U(x)=int |x-y|(P_εη−η)(dy)=int_0^ε (P_uη)(x),du .
]
The second identity uses the Brownian heat equation and the second distributional derivative of the absolute value. By compact support in [2,2.1] and ε<1,
[
 U(x)≤C m sqrt ε,exp[-mathrm{dist}(x,K)^2/(2ε)]
       ≤ C_1 m,g(x)/(1+x^2).
]
The constants may depend on ε and the fixed containing interval, but are independent of the copy size. Here and below g denotes the standard Gaussian density.

For the small-copy bound put ℓ=A√m/(1+|y|) with a fixed sufficiently large A. On [y−2ℓ,y+2ℓ] the Gaussian density is comparable to g(y) by absolute constants once m is sufficiently small. Each of [y−ℓ,y] and [y,y+ℓ] contains a point with r≤ℓ². Otherwise Brownian paths starting inside that interval and remaining there have no stops through time ℓ². The initial density ≥g/2 dominates a constant times g(y) times the positive first Dirichlet sine eigenfunction there. Heat evolution killed at the two endpoints gives a surviving density lower bound by that sine times exp(−π²t/(2ℓ²)); integration for t≤ℓ² gives U≥c g(y)ℓ² on the middle third, contradicting the upper bound for large A. This comparison also works as an inequality of occupation measures, without invoking a boundary density.

At time ℓ², if y is still in D_{ℓ²}, its open component I is thus bounded and has length D≤2ℓ. At any time t, a bounded component I of D_t contains exactly λ(I) of surviving mass. No earlier stop lies there. Paths inside cannot cross the two closed endpoints afterward and paths outside cannot enter; endpoint stops have probability zero since λ has no atoms.

Write δ for the diameter of K. Every interval L satisfies
[
 |L∖K|≥cmin(|L|,|L|²/δ).
]
For verification, if the central half contains no K the linear bound holds. Otherwise follow the nested Cantor intervals about a point of its central half until their length first drops below |L|/4. For |L|≤4δ the selected interval has length comparable to |L|, because successive ratios are bounded above by 1/2 and below by a positive constant. Stage j lengths are comparable to δ2^(−j) since the retained infinite product is positive. The next middle gap consequently has length at least c|L|²/δ and lies inside L. Larger intervals use |K|≤δ. Constants persist under affine scaling.

On I the density g is comparable throughout. Every subinterval of length D/2 therefore has λ-mass ≥c min(1,D/δ)λ(I). Brownian survival inside a static interval of width D over a duration h is bounded by C exp(−c'h/D²), uniformly over starting points, for example by the Markov property at successive times D². Any surviving subcomponent of length >D/2 at time t+h would have at least the preceding mass, but could have at most this survival fraction of λ(I). Thus after a duration C D²(1+log_+(δ/D)) all subcomponents have width ≤D/2. This comparison can be iterated because the density comparability and thickness constants hold throughout the original interval. The function D²(1+log_+(δ/D)) increases with D, and summing over successive halvings costs at most C D²(1+log_+(δ/D)). If a point still belonged to D_t after these times, its open neighborhood would have survived all the width bounds. This proves
[
 sup_y r(y)≤C(m+δ²)
]
for sufficiently small affine copies. For the fixed unscaled set, finite boundedness follows by the same argument: in both far tails use ℓ=1/(1+|y|) and the first, exponentially decaying bound on U to get the two comparison points. Trap the remaining compact region between one such point on each side at a common finite time. Its components have uniformly bounded width and a positive Gaussian density ratio, so the same halving proof gives a finite bound. This proof of bounded duration does not require μ≥η, which fails on K.

**Assembling the single graph.** Assign the point x the value s on K and 1+r(x) outside. These branches are carried by the mutually singular measures η and λ. Bounded-time optional sampling for the residual transport gives ∫exp(zx−z²r(x)/2)λ(dx)=∫exp(zx)μ(dx). Multiplying by exp(−z²/2) and adding the early term cancels P_εη exactly. The Gaussian spatial marginal and all real identities follow. At s the heat remainder has density p_s−1_K g, strictly negative on K. Up to any t<1 it evolves freely from this signed slice. For the unscaled K the smooth pairing from round 153 is at least as negative because the early sink has increased from (15/16)η to η and the heat kernel is positive. A genuine stopped Brownian law would leave the nonnegative distribution of its surviving paths at such t, a contradiction. Finally the open complement accumulates at every point of K while the time gap is at least ε, proving the discontinuity claim.
