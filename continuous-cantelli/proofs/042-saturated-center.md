# Origin saturation forces an unbounded signed discrepancy, conditionally on exactness

Assume the hypotheses of an exact graph, global Hölder continuity, positivity away from zero, and the following explicit positive-side contact for small u>0:

    f(u)=u²/(2D(ℓ)),   ℓ=log(1/u),
    D(ℓ)=ℓ+(1/2)log ℓ+d+(j log ℓ+q)/ℓ,   d=log√2.

Then for T→∞,

    λ([0,T])−T ∼ −T/(2 log T).                                   (2)

In particular, replacing the uniform positive lower bound in the O(1) theorem by positivity only away from zero would force exclusion of every such explicit saturated contact in an exact graph. This statement does not claim existence of such a graph.

Indeed, on the positive origin branch c=2D/u is decreasing in u for sufficiently small u. Changing variables gives its raw-center density at large c,

    k(c)=exp(D−ℓ−u²/2)√(2D)/(2(D+D')).

With L=log c, we have ℓ=L+O(log L). Expanding yields

    k(c)=1+(j−1/4)(log L)/L+[q−1−d/2]/L+o(1/L).

To see this, the product equals √(ℓ/D)·exp[(j log ℓ+q)/ℓ−u²/2]·(1+D'/D)^(−1), with D'/D=1/ℓ+o(1/ℓ). Round 29 gives the delivered density of this same positive core as 1+(j−1/4)(log L)/L+[q−(1+d)/2]/L+o(1/L). Their difference is −1/(2L)+o(1/L), which integrates to −T/(2log T)(1+o(1)).

For x≥r>0 and x≤−r, global Hölder continuity forces the first-round deficit bound, hence f bounded below on tails; compact positivity covers the rest. These sources in each half have bounded reciprocal variance, so each component's raw versus delivered discrepancy is O(1) by the tile argument applied to the submeasure (its delivered density is at most 1). On −r<x<0 raw centers are negative. Their delivered mass in [0,T] is ≤C∫₀ʳ min(T,1/u)du=O(log T), using exp(−tu−t²f/2)≤exp(−tu). These errors are negligible on the T/log T scale, yielding (2). For the same calculation the exact smooth model can be replaced by a differentiable expansion having an o(1/ℓ) remainder and derivative o(1) added to D'.
