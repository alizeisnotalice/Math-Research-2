# Uniform quadratic flatness in measure at the support of a positive-mass top layer

Let 0≤f≤H=sup f be a continuous exact graph, v=H−f, A={v=0}, and K=supp(1_A dx), the closed essential support of A. Let g=p₁. For each compact K₀⊂K, R>0 and η>0,
$$
 \lim_{h\downarrow0}\ \sup_{a\in K_0}
 h^{-1}|\{x:|x-a|\le Rh,\ v(x)\ge\eta h^2\}|=0 .
$$
The center is allowed to vary with the spatial scale. If A has no positive-measure part, K is empty and the assertion gives no restriction.

**Continuity of the positive-variance contribution on K.** Set
$$
 J(a)=\int_{\{v>0\}}p_{v(x)}(a-x)g(x)\,dx .
$$
The output density identity at variance H gives J+1_A g=p_H almost everywhere. The function J is lower semicontinuous by Fatou. Approaching any point a∈K through full-measure good points of A yields J(a)≤p_H(a)−g(a). Conversely the contact-phase identity in round 21 gives J(a)=p_H(a)−θ(a)g(a) with 0≤θ(a)≤1 at every a∈A. Hence equality holds on all of K, and J restricted to K is continuous. On K₀ the continuous truncations integrating only over v≥1/n increase to J, so convergence is uniform by Dini. Each truncated integrand is bounded by √(n/(2π))g(x). It follows that
$$
 \sup_{a\in K_0}\int_{0<v(x),\,|x-a|\le r}p_{v(x)}(a-x)g(x)\,dx\longrightarrow0 .
$$

**Uniform blow-up estimate.** Fix B>R+1. In the exact smoothed density at a, the part |x−a|>Bh differs from the corresponding J integral there by ≤C/B². Indeed the variance derivative of p_w(u) has absolute value ≤C/|u|³ uniformly for u≠0, also down to w=0 with the off-diagonal zero convention. Integrate this bound times h²g(x). Thus the inner smoothed integral on |x−a|≤Bh is at least g(a)−C/B²−o(1), uniformly in a, by J(a)=p_H(a)−g(a) and continuity of the target. On that inner interval split the variances into v<δh², δh²≤v≤Mh², and v>Mh², with 0<δ<η<M. The middle set has length o(h) uniformly for fixed B,δ,M, since p_v(a−x)≥c(B,δ,M)/h there and the unsmoothed near integral tends uniformly to zero. The final group contributes at most CB/√M to the smoothed density. On the first group the smoothed kernel differs from p_{h²}(a−x) by at most Cδ/h. The full p_{h²} integral over the line against g is g(a)+o(1). In comparing with that full integral, the omitted group v≥ηh² on |x−a|≤Rh costs at least c_R/h times its Lebesgue measure (g has a positive uniform lower bound there). Consequently the limsup in the theorem is ≤C_R(B⁻²+B/√M+Bδ). First send M to infinity and δ to zero at fixed B, and then send B to infinity.

**Consequences at gaps.** For any sequence of bounded complementary gaps (a_j,b_j) of A with a_j∈K₀ and lengths ℓ_j→0, the distributions of v/ℓ_j² under uniform measure on the gaps converge in probability to zero. The same holds with right endpoints in K₀. For a fixed gap with such an endpoint this also applies on shrinking segments adjacent to that endpoint. If A equals its essential support near a gap endpoint, a lower bound v(x)≥c dist(x,A)² throughout the adjacent gap segment is impossible: on half of [a+h,a+2h] the ratio to h² already exceeds a fixed positive number, for all small h. More generally the uniform statement forbids a fixed positive fraction of each of a sequence of such shrinking gaps from carrying a fixed positive quadratic-scale deficit, even when there is no single fixed endpoint at which these fractions recur. Exceptional sharp spikes on fractions tending to zero are not excluded by this estimate.
