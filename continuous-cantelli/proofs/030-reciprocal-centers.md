# Raw-center discrepancy for the reciprocal Gaussian representation

Let f be a continuous, bounded nonnegative exact graph normalized by ∫exp(tx−t²f(x)/2)γ(dx)=1 (t∈R). The zero set has Lebesgue measure zero: any positive-measure subset in x≥δ>0 or x≤−δ makes the identity diverge along one tilt. On f>0 set

    c=x/f(x), s=1/f(x), w=f(x)^(−1/2)exp[(1/f(x)−1)x²/2],
    λ=c_#(w dx).

Then, for T→∞ on either reflected tail,

    |λ([0,T])−T| = O(T/√log T).                                  (1)

In fact the same rate holds uniformly for all intervals I contained in [−T,T], with |λ(I)−|I|| on the left. If f is bounded below on compact subsets of R (e.g. f>0 everywhere), the bound improves to o(√T). If inf_R f>0, interval discrepancy is O(1) uniformly over all intervals.

The result concerns the raw (unsmoothed) center measure, including arbitrarily small times and origin zeros. The reciprocal identity itself is immediate by completing the square in t.

**Transport and packing.** With p_s the Gaussian density, direct algebra gives

    ∫ w(x)p_s(t−c)dx=1.

Thus for every measurable interval J the mass delivered into J is precisely |J|. Whenever an input set has c in an interval J of length L and √s≤L, its raw weight is ≤3|J|/P(|Z|≤1), since its probability of landing in the concentric interval 3J is ≥P(|Z|≤1). All integrations with nonnegative integrands use Tonelli. No a priori local finiteness of λ is needed.

**Small source coordinates.** Set a=(log T)/(4T), for T>e. For all intervals I⊂[−T,T], the total delivered weight from |x|≤a is bounded, using the original kernel expression, by

    C T a exp(Ta) = O(T^(1/4)log T).

The raw weight from that same region whose c∈I is at most

    2∫₀ᵃ √(T/u)exp(Tu/2)du = O(T^(1/8)√log T),

since |c|≤T entails f≥|x|/T. The single point x=0 has zero source measure. These estimates are uniform in I.

**Central centers.** On |x|>a, |c|≤2T the standard deviation is ≤H=√(2T/a). Tile the real line by intervals of length H. By packing, every tile carries at most C H of raw mass from this group. Conditional on center c, the absolute difference between the raw indicator 1_I(c) and the conditional probability to land in I is bounded by the probability of crossing either endpoint b of I, hence by ∑_{b∈∂I} exp(−(c−b)²/(2H²)). Summing over tiles costs C H since the Gaussian tail sum about each endpoint is uniformly bounded. This also works if an endpoint carries raw mass. The estimate establishes local finiteness.

**Distant centers.** On |x|>a and R≤|c|≤2R, R=2^j T (j≥1), one has s≤2R/a. For T large the standard deviation ≤R, so the two shells have raw mass ≤C R by packing. Their distance from I is ≥R/2. Their total smoothed mass leaking into I is at most C R exp(−aR/16). Summing this geometric dyadic bound gives O(T exp(−Ta/8)), since for z→∞ one has ∑_{j≥1}2^j e^(−2^j z)≤C e^(−2z). The raw indicators are zero in these shells, up to endpoint conventions of a disjoint partition. Combining the three groups gives the uniform (1).

**Strictly positive graphs.** Fix R₀≥1. On the compact |x|≤R₀, positive f has a positive minimum, so total raw weight there is finite and bounds the contribution to the discrepancy uniformly in T. On the remaining groups repeat the preceding proof with a=R₀, obtaining O(√(T/R₀))+O(T exp(−TR₀/8)), with universal implied constants. Divide by √T, take T→∞, then R₀→∞. If instead inf f>0, use the uniform scale H=(inf f)^(−1/2) globally; the central-center tiling proof then gives O(H) on every interval. This also applies to a strictly positive globally Hölder exact graph, since the first-round mass bound implies (1−f)_+=O(log|x|/x²) on both tails. A globally Hölder graph positive only off the origin may still have unbounded reciprocal scales near zero; the O(1) conclusion here does not cover that case.

**A related integrated sign identity.** Tonelli and integration in t of the differences against the constant graph give the parameter-free budget

    ∫_{f<1} (w−1) dx = ∫_{f>1} (1−w) dx,

where both sides may be infinite. Indeed for each t the signed difference between the two kernels has sign 1−f at each source x, and the difference integrated over all t at fixed x equals w(x)−1. This identity does not require an a priori absolute-integrability assumption.
