# Nonradial contact with arithmetic Cantor sets: a real frequency exclusion

Write v=H−f. Let E be a compact null set of highest points, isolated from other zeros in a neighborhood. Assume that there are 0<A_-≤A_+<∞ and β>1 such that for all x at sufficiently small positive distance r from E,
$$
 A_-[\log(1/r)]^{-\beta}\le v(x)\le A_+[\log(1/r)]^{-\beta}.
$$
There is no equality assumption across gaps or addresses. Two exclusions follow.

1. Suppose there are k_m→∞ and real phases θ_m with cos(k_m x−θ_m)≥c>0 on E, uniformly in m. Then no exact graph has this contact. This includes every affine copy of the two-end-digit Cantor set in integer base B≥4: for the normalized set x=(B−1)Σ e_n B^(−n), e_n∈{0,1}, take k_m=2πB^m/(B−1). Modulo 2π the phases lie in [0,2π/(B−1)], a strictly shorter interval than a semicircle.
2. For the ordinary middle-thirds set (or any affine copy), put d=log 2/log 3. The same exclusion holds whenever
$$
 A_+/A_- < [2/(1-d)]^\beta .
$$
In particular an arbitrary address-dependent relative error tending uniformly to zero in the leading logarithmic contact is allowed. The same result applies with constant multiplicative bounds whose ratio satisfies this strict condition.

**Localization.** Take a bounded open source neighborhood U={dist(x,E)<r₀} with closure in the indicated collar, choosing r₀ sufficiently small. Its positive-variance output has an analytic density D_U on a common horizontal strip with integrable traces on a smaller strip. Indeed near E the omitted sources are separated in space, so their heat kernels are uniformly bounded and holomorphic on a fixed complex neighborhood, even for arbitrarily small variances. Zero-age omitted mass is absent there. Subtract those sources from the analytic target p_H. Away from E the included sources sufficiently close to the observation position have a uniform lower variance bound, and the other included sources are spatially separated. Compactness gives a common strip over bounded real parts; outside a sufficiently large compact set the included means are bounded and all variances ≤H, giving a uniform integrable Gaussian tail on the strip. The extensions agree by uniqueness. Consequently for some c₀>0,
$$
 \left|\int_U e^{ikx}e^{-k^2v(x)/2}g(x)\,dx\right|
 =|\widehat D_U(-k)|\le C e^{-c_0 k}\quad(k>0).
$$
This uses the same local continuation principle as round 23; the maximum set may be extended in space. Throughout let α=2/(1+β)<1. Constant affine changes of spatial scale add O(1) inside log(1/r) and divide the chosen frequencies by the scale; strict coefficient margins may be slightly widened without affecting the hypotheses.

**Strict phase margin.** For k=k_m the rotated real part is nonnegative and ≥c/2 for all r<c/(2k). The negative part can only come from r≥c/(2k), where damping bounds its absolute integral by C exp(−c' k²/(log k)^β). Fix any point of E. In an interval of radius exp(−k^α) about it the phase is uniformly positive, while v≤A_+ k^(−αβ) off E. The Gaussian input is bounded below there. This contributes at least exp(−C' k^α), contradicting localization.

**Middle-thirds, including the borderline semicircle geometry.** Use coordinates [0,1] and k=π3^m. The phase sin(kx) is nonnegative throughout the 2^m surviving intervals at depth m, because their left endpoints are even integer multiples of 3^(−m). Negative phase points are in the coarser removed gaps or outside the hull. On these gaps and the two exterior collars, |sin(kx)|≤k r whenever its sign is negative: the nearest endpoint of the gap is a zero of sin(kx) for every gap removed through level m. There are O(2^m) such intervals. Up to polynomial factors in k their total negative contribution is at most
$$
 \int_0^{r_0} r\,\exp\!\left[-\frac{A_- k^2}{2(\log(1/r))^\beta}\right]dr
 \le \exp\big[-(J(2,A_-)-o(1))k^\alpha\big],
 \qquad J(q,A)=\min_{u>0}\left(qu+\frac{A}{2u^\beta}\right).
$$
Here the Gaussian spatial weight is bounded. Substitute log(1/r)=k^α u: on a fixed bounded u-interval use the minimum, and for u beyond a sufficiently large fixed cutoff integrate the tail exp(−q k^α u) directly. Polynomial factors do not change the rate.

For the positive part, choose in each depth-m surviving interval the same fixed subcylinder in its normalized coordinates, for instance [2/9,1/3]; there sin(kx)≥sin(2π/9). At any much deeper depth n select, in all descendant intervals of these subcylinders, central portions of their next removed gaps. Their total length is ≥c (2/3)^n with a fixed c>0, and log(1/r)=n log 3+O(1) uniformly there. The input density has a positive lower bound on the hull. Given u>0 choose n log 3=u k^α+O(1), so n≫m. This gives a positive contribution ≥exp(−[(1−d)u+A_+/(2u^β)+o(1)]k^α). Optimize in u. The strict condition is precisely J(1−d,A_+)<J(2,A_-), since J(q,A) is a common constant times (q^β A)^(1/(1+β)). The positive part therefore dominates the negative part and contradicts the exponential strip bound. This real-frequency proof does not require control of complex phase errors from comparing different variances at different addresses.
