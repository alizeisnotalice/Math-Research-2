# Positive critical Gevrey type and endpoint refinements for the quartic diagonal

Use the normalization of round035-general-quartic-controlled-jet.md: k even, m=k/2, c=w(0)>0, b=w/c, g=p₁(x₀+·), β_i=c^{-i}2^i(m)_i/(2k)_{2i}, h=(4ck)^{-1/2}, λ=|b₁|/(4√c). The diagonal of the potential divided by its factor L_k is ∑_{0≤i<m}β_i d_i with d_i=[u^{2i}]g b^{m−i}. Define
τ=limsup (|[u^j]w|/√(j!))^{1/j}.
Let Q be the unique root in (0,1) of 1−Q=(1−Q/2)^4, Q≈0.912621974616, and K=2/[√Q(1−Q/2)]≈3.85064940161.

**Claim.** The original positive diagonal bound of order exp(λ√k) still holds when τ<K√c. At τ=K√c it also holds under either
(i) |[u^j]w|≤(K√c)^j √(j!) exp(−D√j) eventually for every fixed D>0, or
(ii) λ>0 and |[u^j]w|≤(K√c)^j √(j!) exp(ε√j) eventually for every ε>0.
All constants here concern Taylor coefficients, not derivatives. Combined with the independently positive exterior coefficient, these are exclusions for an exact graph.

## Weight geometry

With j=2i, q=j/k and v_i=β_i^{1/j}, factorial ratios give for q in any fixed [q₀,1], uniformly,
log(√c v_i √(j/e))=log F(q)+O_{q₀}(log k/k),
F(q)=√q/2 · exp{[−(1−q)log(1−q)/2+2(1−q/2)log(1−q/2)]/q}.
The convention is 0 log0=0. This follows by applying log(L!)=L log L−L+O(log(L+2)) to m!, (m−i)!, (2k)! and (2k−2i)!, including m−i bounded. Differentiation gives q²(log F)'=(1/2)log(1−q)−2log(1−q/2). The sign changes exactly once from + to −: on 0<q<1 the derivative of log[(1−q)/(1−q/2)^4] has sign 1−3q/2. The asserted value at the maximum follows by substituting its defining equation. At q→0 use β_i/h^{2i}≤C uniformly: each factor (1−2l/k)/[(1−l/k)(1−(l+1/2)/k)]≤(1−1/[2(k−l)])^{−1}. Their product is bounded uniformly for l<i≤k/2 by the harmonic sum. Thus v_i ≤h exp(O(1/j)), so √c v_i√(j/e)≤C√q. This also gives v_i/h≤1+O(1/k) on any fixed macroscopic band q≥q₀.

## Low indices

Fix T large as specified below and let P be b truncated to J=min(k,ceil(2T√k)+2). Then d_i is unchanged for i≤T√k. For every fixed C₁, P^+(C₁h)=1+|b₁|C₁h+O(h²), with uniform derivative bound on the same disk. Indeed |b_l|≤C A^l√(l!) for a fixed finite A; here A C₁h√J→0. Even the endpoint conditions satisfy this hypothesis with enlarged A. For i≤T√k, β_i/h^{2i}=1+O_T(k^−1/2) by taking logarithms of the product. The sum of h^{2i}|d_i| over this range is ≤g^+(h)P^+(h)^m≤C exp(λ√k). Thus one may replace β_i by h^{2i} here with o(exp(λ√k)) absolute error.

The inverse at zero of z=u√P(u), defined on |z|≤2h, gives A_k(z)=g(u)(zu'(z)/u)P(u)^m. Rouché on |u|≤4h and the derivative bound imply a unique analytic inverse; [z^{2i}]A_k=[u^{2i}]gP^{m−i} by residue. On |z|≤2h one has |A_k|≤C exp(C₀√k), where C₀ may be chosen independent of any fixed T once k is sufficiently large depending on T: the first two Taylor terms determine the √k exponent. The tail i>T√k of the even series tested at h is therefore O(exp(C₀√k)4^−(T√k)). The low sum equals (A_k(h)+A_k(−h))/2 up to o(exp(λ√k)) for sufficiently large fixed T. The real Jacobian prefactors tend to g(0)>0, and m log P(u(±h))=±m b₁h+O(1), proving the claimed positive main lower bound for these indices.

## Mid and upper indices

For i>T√k, analyze d_i with the separate truncation b_{≤j}, j=2i. With r_i>0,
β_i |d_i| ≤ (v_i/r_i)^j · g^+(r_i)(b^+_{≤j}(r_i))^{m−i}.

In the strict case choose ρ>1, A>τ so that ρ A<K√c. Use r_i=ρv_i. Weight geometry yields A r_i√(j/e)≤1−δ uniformly for all sufficiently large k and i>T√k. The factorial estimate √(l!)≤C l^{5/4}(l/e)^{(l−2)/2} implies
∑_{l=2}^j |b_l|r_i^l ≤ C r_i² ∑_{l=2}^j l^{5/4}(A r_i√(l/e))^{l−2} = O(r_i²).
Thus |d_i|β_i≤C exp(C₁√k)ρ^{-2i}; C₁ can be taken independent of T. Choose T large enough to make the whole sum o(exp(λ√k)).

For the endpoint, fix a closed small band around Q inside (0,1), say q∈[q₁,q₂]. On its complement there is a strict gap in the maximum of F. The previous argument still applies there with A>K√c and ρ>1 sufficiently close to one. Within the band write τ*=K√c. Uniformly τ* v_i √(j/e)≤exp(C log k/k). In case (i) take r_i=v_i exp(D/√k) for a fixed large D. For l≥2, excess of (τ*r_i√(l/e))^{l−2} over one is bounded by exp[(D+1)√l] for sufficiently large k (use l≤k and (C log k)/√k→0). The assumed damping in the coefficients then makes the previous displayed sum O(r_i²), choosing damping exponent D+3. Here (v_i/r_i)^j≤exp(−D q₁√k). Choosing D large dominates the bounded exponential C₁√k coming from (b^+)^{m−i}, as r_i/h≤1+o(1) for every fixed D and the leading constant in C₁ is independent of D.

In case (ii), choose 0<η<λ and r_i=v_i exp(−η/√k). For l≤j/2, τ*r_i√(l/e)≤(1+o(1))/√2, so the l^{5/4}exp(o(√l)) Taylor tails are summable. For l≥j/2≥q₁k/2 the same factor raised to l−2 is ≤exp[−(η/2)(l−2)/√k] and hence ≤C exp(−c_η√l); choose the allowed coefficient upper exponent smaller than c_η/2. Thus in this case too b^+_{≤j}(r_i)≤1+|b₁|r_i+O(r_i²). Since r_i/h≤1+o(1) and m−i=(1−q)m, restoring (v_i/r_i)^j gives
β_i |d_i| ≤ C exp{[(1−q)λ+ηq+o(1)]√k},
uniformly on the band. The right side, summed over at most k indices, is o(exp(λ√k)). This proves both endpoint refinements. The constant K measures the maximum in this absolute high-index estimate; optimality as a rigidity threshold is open.

## A polynomial endpoint with zero linear term, and its limitation

If λ=0 and |w_l|_{coef}=o((τ*)^l√(l!) l^−7/4), the same positive diagonal conclusion holds at τ=τ*. Normalize w by c in the estimates. In the fixed interior band q∈[q₁,q₂] about Q the factorial ratio actually satisfies log(√c v_i√(j/e))=log F(q)+O(1/k), j=2i. The half-log prefactors from each factorial ratio are bounded because their arguments have ratios uniformly bounded above and below. The maximum of log F is nondegenerate since the derivative of log[(1−q)/(1−q/2)^4] is strictly negative at Q>2/3. Use
r_i=v_i F(Q)/F(q) · (√k M_k)^{1/j},
where M_k→∞, M_k≤log(k+1), is chosen so slowly that M_k times sup_{l≥q₁k/2}|b_l|l^{7/4}/[(τ*)^l√(l!)] tends to zero. Then r_i=O(k^−1/2), and
(v_i/r_i)^j ≤ exp[−c₁k(q−Q)²]/(√k M_k)
after narrowing the band if necessary. For l≤j/2 the geometric tail in the Taylor absolute majorant is unchanged. On j/2≤l≤j we have
(τ*r_i√(l/e))^{l−2}≤C√k M_k exp[−c₂(j−l)]
using l/j≤exp[−(j−l)/j]. When estimating the Taylor tail divided by r_i², the factor l^(5/4) from Stirling is multiplied by l^−7/4 and the small coefficient envelope; l^−1/2√k=O(1). Summing the resulting geometric series in j−l bounds this tail uniformly. As b₁=0, (b^+(r_i))^{m−i}=O(1). The summed band contribution is O(1/M_k), since a Gaussian in q−Q has lattice sum O(√k). The rest of the proof is unchanged.

The strict small-o at this endpoint has content for a uniform positive-diagonal proof. To see this fix c=1 and g=p₁(x₀+·); choose b=1 plus sparse negative blocks of even-order Taylor coefficients. A block at a large even k puts
b_l=−C K^l √(l!) l^−7/4   if l even and |l−Qk|≤√k.
Separate the chosen k so that each new block begins beyond all previous observed orders, and choose its k after fixing its finite predecessor polynomial P=1+O(u³). For the normalized diagonal at this k all block coefficients enter linearly because their indices exceed k/2. For each l in the block its multiplier is
Σ_{i≥l/2}^{m−1} β_i(m−i)[u^{2i−l}]g P^{m−i−1}.
The leading i=l/2 term is β_{l/2}(m−l/2)g(0). Consecutive β ratios are ≤C₀ h² for all i (the explicit product factors have a uniform upper bound). After division by β_{l/2}(m−l/2) the absolute sum of the remaining terms is ≤g^+(√C₀ h)(P^+(√C₀ h))^m−g(0)=o(1); the convergence is uniform in l for fixed P since m h³→0. The earlier low-index proof gives normalized diagonal of P→g(0), since P=1+O(u³). Throughout this central block β_{l/2}K^l√(l!)≥c₃ k^1/4 by interior-uniform Stirling and log F(q)≥log F(Q)−C₄(q−Q)². Also m−l/2≥c₅ k. Each negative term contributes at most −c₆ C k^−1/2; there are ≥c₇√k terms. Choose a fixed sufficiently large C, then at each stage choose k sufficiently large depending on the previous polynomial. This gives an infinite sequence with strictly negative normalized diagonal. Future blocks do not affect these orders. Borel's theorem on smooth jets produces a smooth b with those coefficients, positive near zero since b(0)=1. This disproves a universal eventual positive diagonal bound for arbitrary coefficients satisfying only the endpoint big-O with unspecified constant. It does not imply existence of a Gaussian-output variance graph; other coefficients and the exterior constraints still apply.
