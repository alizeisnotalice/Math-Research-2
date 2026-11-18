# Fourier comparison with a radial template in averaged error, without arithmetic resonances

Let E be any strongly separated homogeneous Cantor set of dimension 0<d<1 on the line, with equal contraction ratio ρ and m pieces, as in round 185. Fix a bounded source collar U={dist(x,E)<r₀}, containing no other zero of the continuous variance deficit v. Write r=dist(x,E), Q=1−d, and ℓ=log(1/r). Suppose V(r)ℓ^β→A>0 and v(x)≥(A_-−o(1))/ℓ^β uniformly as r↓0 with 0<A_-≤A and β>0. Assume that for some η>0 and all large integers n,
$$
 \int_{n\le\ell<n+1}|v(x)-V(r)|\,dx\le C e^{-(Q+\eta)n},
 \qquad A_-(Q+\eta)^\beta>AQ^\beta .
$$
Take V positive continuous away from r=0 in the collar, or restrict to a smaller collar where the asymptotic gives the needed positive lower bound on each separated compact distance interval. These hypotheses allow address-dependent changes on a set of relative volume O(exp(−ηn)) in the nth shell with bounded amplitude; they also follow from a uniform power error O(r^η) and the stated lower bound. Uniform relative o(1) alone does not imply the displayed shell inequality.

For β>1 this local configuration is impossible for an exact graph. For β=1 it is impossible when the distance to all other zeros exceeds 2J(Q,A), including when E is the only zero set; for any β>0 a uniform exterior positive variance bound suffices. Here J(q,A)=min_{t>0}(qt+A/(2t^β)), α=2/(1+β). With arbitrary tails one obtains the exterior mass requirement of round 187 along a sequence of real frequencies.

**Shell rates.** A shell of weight at most a polynomial in n times exp(−qn), with damping exp(−k²(A'−o(1))/(2(n+1)^β)), sums to ≤exp(−[J(q,A')−o(1)]k^α) as k→∞. For a direct upper bound split at n=B k^α, sum the exponential tail beyond it and use the minimum for each of O(k^α) smaller terms, first shrinking A' by an arbitrarily small fixed factor and then sending that factor to one. B can be chosen arbitrarily large. Fixed bounded n contribute exp(−c k²). The shell hypothesis and the mean value bound for exp(−k²v/2) show that replacing the true damping by the template costs at most a polynomial k² times the above estimate with (q,A')=(Q+η,A_-). Its exponential rate is strictly smaller in magnitude of contribution than the main term below.

**Address comparison using real frequencies.** Let μ be the uniform self-similar law and σ=gμ. On a basic gap G, positions x=S_ω(u) at level j have r=ρ^j dist(u,∂G), independent of the word ω. The average of g(S_ω(u))exp(ik S_ω(u)) over all words differs from ∫g e^{ikx}dμ by at most C(1+k)ρ^j: couple S_ω(u) to S_ω(Y) for Y∼μ. Thus the template integral on the hull inside U equals σ̂(−k)Z(k) plus an error, with
$$
 Z(k)=\sum_{G,j}(m\rho)^j\int_{u\in G:\,r<r_0}e^{-k^2V(r)/2}\,du>0 .
$$
Its error has an extra factor C(1+k)ρ^j inside this sum. The geometric weight (including this factor apart from C(1+k)) in an r-shell of logarithmic index n is O(e^(-n)): put s=j log(1/ρ); at fixed j the allowed basic-gap endpoint strip has length ≤C e^(−n+s) and the coefficient is exp(−(2−d)s), giving ≤C e^(−n)exp(−Qs). Sum over j. The two outside-hull parts of U have shell length O(e^(−n)) as well. Their contribution and the address error therefore have decay rate at least J(1,A)>J(Q,A). Meanwhile by choosing j log(1/ρ)=t k^α+O(1) and u in a fixed middle part of a basic gap, one obtains Z(k)≥exp(−[J(Q,A)+o(1)]k^α) after optimizing t.

The nonzero compact singular measure σ has |σ̂(−k)|≥1/k along some sequence k→∞ by Plancherel. Both errors have a strictly better exponential rate, so the true compact source has absolute Fourier transform at least exp(−[J(Q,A)+o(1)]k^α) on that sequence. The strip localization and exterior alternatives in rounds 186–187 now give all conclusions. The argument uses neither an arithmetic ratio nor positivity of a chosen spatial phase. The exponent cost is obtained from the shell size and a one-dimensional minimum.
