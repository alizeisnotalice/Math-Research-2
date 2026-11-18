# Arbitrarily deep time-difference obstruction for deterministic time marks

For each N≥1 there is a nonconstant bounded Borel function f, continuous off a compact nowhere-dense set of arbitrarily small positive Gaussian mass, such that the spatial marginal is Gaussian, every real tilt is exact, and

    (−∂q)^n R(q,y)>0   for all n≤N, q≥0 and y∈ℝ,

while a positive smooth spatial test gives a strictly negative heat remainder at some fixed time for this graph. In particular no integer difference of order k≤N at any base j can detect this graph; some finite higher order does detect it. The conditional time law is a point mass at every spatial location.

Here are the dependencies and parameter order. Round 142 gives a fixed large-q cutoff Q from the mark range [19/20,21/20] and the finite list of orders. Round 128 gives a positive upper-mark tolerance at that Q. Thus there is d_N>0 such that every exact model with |F−1|≤d_N passes these signs simultaneously over the whole unbounded spatial and parameter domains. First choose ε∈(0,min(d_N,1/8)). The construction of round 154 works with this ε held fixed: put s=1−ε, pick a short interval inside [2,2+ε^(3/2)/100] and an affine copy K of the same thick-complement fat Cantor set. The Gaussian ratio on this interval satisfies p_s/g<1, while for small diameter δ the residual μ=γ−P_ε(1_Kγ) is at least γ/2 by the positive-variance Gaussian kernel bound. All constants in the occupation and component-halving estimates may now depend on this fixed ε. Those estimates give sup r≤C_ε(m+δ²), m=γ(K), for all sufficiently small copies. Choose δ next so sup r≤d_N and, if desired, m lies below any prescribed positive tolerance. The resulting deterministic graph has the required range. This argument does not exchange the ε and δ limits.

On K the time-s slice is strictly negative; smooth testing and right continuity give a negative test at a time in (s,1). The compactly supported time trace against that test has a negative interval in time, so the Hausdorff difference criterion from round 123 yields some negative integer difference. The positive derivative bounds through N force its order to exceed N. This transfers the all-parameter obstruction previously obtained for mixed marks to deterministic graphs. Continuity at points of K is still missing.

There is a sharp order of magnitude for the graph's collapse off K as δ→0 at fixed ε. By the second moment of the bounded residual stopping time,

    ∫ r(x)λ(dx) = ∫x²(λ−μ)(dx)= εm.

Thus sup r≥εm/(1−m), while the upper bound above is O_ε(m), since scaled copies inside the fixed interval have m comparable to δ. The exact graph satisfies E|f−1|=2εm. Consequently shrinking K to a single point at fixed ε cannot yield a nonconstant continuous limit with local uniform convergence: the entire law tends to the unit mark in L¹ while the prescribed downward displacement on K stays equal to ε.
