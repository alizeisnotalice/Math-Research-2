# Exponential height budget inside a bounded interval

Assume a bounded exact continuous graph f in the normalization of round 169. Write g for the unit Gaussian density and p_s for variance s. Let s>0 and let I=(a,b) have length L. Assume either E(a),E(b)≤s, where E is the last activity height, or that the whole strict future {f>s} is contained in [a,b]. The second condition is stated only in terms of f. Then, writing φ(x)=sin(π(x−a)/L) on I,
$$
 \int_{I\cap\{f>s\}}\!\!\phi(x)\exp\!\left(\frac{\pi^2(f(x)-s)}{2L^2}\right)g(x)\,dx
 \le \int_I \phi(x)p_s(x)\,dx .
$$
For the first hypothesis use the local future identity of round 169 with z=iπ/L, multiply by exp(−iπa/L) and take the imaginary part. The identity holds at complex z by entire continuation (bounded time and Gaussian tails). The resulting left side is exactly ∫_I φ dQ_s; for the inequality use Q_s≤p_s dx and φ≥0 there. For the second hypothesis apply the whole-line future identity; the signed convex-hull support theorem of round 58 places Q_s on [a,b]. Its endpoints have no atoms for s>0. This proves the same assertion. The argument permits a time-s platform.

**Width at a peak.** Suppose for s₀≤s<M=sup f, with s₀>0, every strict superlevel {f>s} is a nonempty bounded interval I_s. The graph is unrestricted below s₀. Write L_s=|I_s|, and let R be the supremum of p_s(x)/g(y) for s₀≤s≤M and x,y in the closure of I_{s₀}. This is finite. Also R≥1, since the preceding inequality at s₀ already forces such a bound. For s<t<M set θ=L_t/L_s. Every subinterval of length θL_s inside I_s captures at least L_s(1−cos(πθ))/π of the first Dirichlet sine's spatial integral: the minimum over translations is attained at an endpoint. Hence
$$
 t-s\le \frac{2L_s^2}{\pi^2}
       \log\!\frac{2R}{1-\cos(\pi L_t/L_s)} .
$$
In particular for every s₀≤s<M one has
$$
 M-s\le \frac{8\log(2R)}{3\pi^2}\,L_s^2 .
$$
Indeed split the levels above s by their widths into bins (2⁻ⁿ⁻¹L_s,2⁻ⁿL_s]. Each bin is an interval in the level variable by monotonicity; any two levels in it have ratio of widths at least 1/2. Its length is therefore at most 2·4⁻ⁿL_s² log(2R)/π². Their sum gives the result even in the presence of plateaus. Thus a peak with connected upper levels and L_s=o(√(M−s)) along a sequence is excluded, regardless of how many branches, extrema or oscillations occur below s₀.

The exponential budget also obstructs disconnected upper levels. If p_s(x)≤ρ g(x) throughout I with ρ<1, the set {f≤s} must carry at least (1−ρ) of the φ-weighted Gaussian mass of I. This follows on dropping the exponential. In particular no connected bounded strict superlevel can be contained in {p_s<g}, since then the strict integral comparison already contradicts the budget. These restrictions need no asymptotic expansion of f at a peak.
