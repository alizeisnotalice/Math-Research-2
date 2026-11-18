# Phase concentration and harmonic switching at a middle-thirds contact

Let E be the middle-thirds Cantor set in [0,1]. Let the input have Gaussian density g, let 0≤v≤H be the remaining heat time in an exact graph with analytic Gaussian output p_H, and assume v=0 on E. For some r₀∈(0,1), β>1 and fixed 0<A₋≤A₊<∞ assume
$$
 A₋[\log(1/r)]^{-β}≤v(x)≤A₊[\log(1/r)]^{-β},\qquad 0<r=\operatorname{dist}(x,E)<2r₀ .
$$
We can decrease r₀. Set U={r<r₀}, α=2/(1+β), k=π3^m, and δ=3^{-m}. For a source set V⊂U write
$$
 Z_V(s)=\int_V e^{-s²v/2}g,\quad dμ_{s,V}=Z_V(s)^{-1}1_V e^{-s²v/2}g\,dx,\quad
 F_V(s)=\int_V e^{isx-s²v/2}g\,dx .
$$
All constants below may depend on the graph, collar, and fixed bounds. No restriction on their ratio is imposed.

**Theorem.** There exist b,η>0 and m₀ such that for every m≥m₀ the following conclusions hold. Let P_m and N_m be the left and right endpoints of the depth-m surviving intervals. Write H_m for the union of their neighborhoods of radius k^{-1}exp(−b k^α). Uniformly for |s/k−1|≤η,
$$
 μ_{s,U}(H_m^c)≤e^{-b k^α}.
$$
At s=k, the mass in each endpoint class is 1/2+O(exp(−b k^α)), and their unnormalized first spatial moments under μ_{k,U} differ by O(exp(−b k^α)). On I_m={|sin(kx)|≥1/2},
$$
 μ_{k,U}(I_m)≤e^{-b k^α},\qquad μ_{2k,U}(I_m)≥1/3-O(e^{-b k}).
$$
In particular, writing q=k²/2 and ν_t(dx) proportional on U to exp(−tv)g dx,
$$
 D(ν_{4q}\Vert ν_q)=\int_q^{4q}(u-q)\operatorname{Var}_{ν_u}(v)du≥b k^α .
$$
It follows that E_{ν_q}v−E_{ν_{4q}}v≥b k^{α-2} and some u∈[q,4q] satisfies Var_{ν_u}(v)≥b k^{α-4}, with constants decreased when needed.

The same conclusions hold *conditionally in each of many individual cylinders*, uniformly in m and in the cylinder, with error O(exp(−b k^α)) in the second-harmonic lower bound. More precisely, there are j₀ and M<∞ such that every depth-j Cantor cylinder J=[a,a+L], j≥j₀, with L=3^{-j}≥M k^{α-1} qualifies. Use V=V_J=(a−L/8,a+9L/8) and only the depth-m endpoint classes descending from J. In the asserted range j≤m for large m. Within each J, let P_i,N_i be the conditional μ_{k,V_J} masses in the two endpoint neighborhoods of its i-th depth-m interval, and a_i its left endpoint. Then
$$
 \sum_i(a_i-a)(P_i-N_i)=\frac{δ}{2}+O(e^{-b k^α}),\qquad
 \frac12\sum_i|P_i-N_i|≥\frac{δ}{2L}-O(L^{-1}e^{-b k^α}).
$$
In particular at the finest guaranteed depth, L comparable to M k^{α-1}, the conditional total variation of the two cell-mass vectors is at least c k^{-α}. These are necessary conditions on all fixed coefficient ratios. They leave open whether R186's ratio hypothesis can be dropped in its exclusion.

**Uniform analytic localization.** R186 gives |F_U(s)|≤Cexp(−c s); the same contour shift with the extra weight 1+|x| gives |F'_U(s)|≤Cexp(−c s). Here F'_U is also the derivative, on the real axis, of the source integral including its Gaussian factor. The unweighted and first-moment trace bounds needed for shrinking J admit the following uniform version. Take j₀ so L<r₀/4. Distinct depth-j cylinders have gap at least L. Put E_J=E∩J. At real observation points y with dist(y,E_J)<L/32, the omitted source positions have distance at least 3L/32 from y. Their positive-age kernel at complex argument y+iz, |z|≤L/128, is bounded by C/L, since
$$
 t^{-1/2}\exp[-((y-x)^2-z^2)/(2t)]≤C/|y-x| .
$$
Their zero-age contribution has no support in this neighborhood. Subtract the omitted contribution from the analytic output there. For dist(y,E_J)>L/64, split included means into those with |x−y|≥L/128 and the rest. The former obey the same bound on the narrower strip |z|≤L/512, and among the latter dist(x,E)≥L/128. Indeed E outside J is separated from V_J by ≥7L/8. For such nearby means v≥A₋[log(128/L)]^{-β} and v≤H, so their complex Gaussian kernel is bounded by C(1+log(128/L))^{β/2}; the factor exp(z²/(2v)) is bounded uniformly in j. These two real neighborhoods overlap. The real densities agree almost everywhere there by exactness; their holomorphic extensions agree throughout the overlap. On real tails |y|>4 the included means lie in a fixed bounded interval, and v≤H gives a uniform Gaussian majorant, by retaining half of the spatial Gaussian for the tail and using the other half to bound the prefactor. Thus the analytic continuation D_V of the local output satisfies
$$
 \sup_{|z|≤L/1024}\int(1+|y|)|D_V(y+iz)|dy≤C L^{-1}(1+|\log L|)^{β/2}.
$$
All holomorphic assertions follow by integrating on compact substrips with the same bounds. Contour shifting both the transform and its derivative yields |F_V(s)|+|F'_V(s)|≤C L^{-1}(1+|log L|)^{β/2} exp(−sL/1024), s>0. For V=U use the fixed R186 strip throughout.

**Tilt mass.** For s∈[k/2,2k] one has uniformly Z_V(s)≥exp(−C₀ k^α). Indeed take an E_J endpoint (any E point for U) and integrate on its exp(−k^α)-radius interval, contained in V. The variance there is ≤A₊k^{-αβ}, and 2−αβ=α. Enlarge M in the shrinking-cylinder condition so |F_V(s)|/Z_V(s)+|F'_V(s)|/Z_V(s)≤exp(−2 k^α) on this whole window. For V=U the quotient is O(exp(−c k)). For s=k choose B₀>0 small so A₋/(2B₀^β)>C₀+4. Then
$$
 μ_{k,V}\{r>e^{-B₀ k^α}\}≤e^{-2k^α}.
$$
For any sufficiently large C₁ and s∈[k/2,2k], μ_{s,V}{v>C₁ k^{α-2}}≤exp(−2k^α), directly from the tilt and Z lower bound.

**Phase and centroid.** The negative part of sin(kx) is bounded by kr on the collar. Indeed the phase is nonnegative on the surviving intervals, and elsewhere the nearest coarse gap endpoint (including the two exterior endpoints) has zero sine. Consequently E_{μ_{k,V}}(sin(kX))₋≤exp(−c k^α). Using the transform bound on its signed mean gives E|sin(kX)|≤exp(−c' k^α). If r≤exp(−B₀ k^α) and |sin(kx)|≤exp(−c'k^α/2), then x is within C exp(−c' k^α/2)/k of a surviving endpoint, for the sine zeros spaced δ in a surviving interval are just its endpoints, and all other grid points have distance ≥δ from E. For V_J such an endpoint must descend from J. Markov's inequality and the r bound prove the neighborhood estimate at s=k with a positive exponent. The estimates remain valid if the exponent in the neighborhood size is decreased. To transfer to nearby s, first discard {v>C₁ k^{α−2}} under both tilts. On the rest the unnormalized weights have ratio between exp(−C₂η k^α) and exp(C₂η k^α), so Z_V(s)/Z_V(k)≥(1−exp(−2k^α))exp(−C₂η k^α). Take η sufficiently small, then decrease b if necessary.

On these endpoint neighborhoods, cos(kx) equals +1 or −1 up to exponentially small error. Re F_V(k) therefore gives the class masses. Also
$$
 \Im F'_V(k)/Z_V(k)=E_{μ_{k,V}}[X\cos(kX)-k v(X)\sin(kX)] .
$$
Since v≤H and X is bounded on V, the sine estimate transfers to the signed first spatial moment. Rounding captured points to their endpoints changes this moment by ≤exp(−b k^α)/k. The total signed mass is exponentially small, so one may subtract a from both signed spatial moments. The right endpoint in cell i is a_i+δ. The stated weighted imbalance identity follows. For d_i=P_i−N_i, the sum of their positive parts is (Σ|d_i|+Σd_i)/2. Using 0≤a_i−a≤L proves the asserted total variation bound.

**Second harmonic.** The exact identity E_{μ_{2k,V}}sin²(kX)=(1−Re F_V(2k)/Z_V(2k))/2 implies μ_{2k,V}(I_m)≥1/3−O(|F_V(2k)|/Z_V(2k)), since sin²≤1/4 outside I_m. For shrinking cylinders the error here is O(exp(−2k^α)); for U it is O(exp(−c k)). Coarse graining relative entropy over I_m and its complement yields D≥p log(1/q₀)−log2, where p=μ_{2k,V}(I_m), q₀=μ_{k,V}(I_m). For Λ(t)=log ∫_V exp(−tv)g dx, bounded v gives Λ''(t)=Var_{ν_t}(v). Substituting the two tilt densities gives D=Λ(q)−Λ(4q)+3qΛ'(4q), the stated integral. Since 0≤u−q≤3q and ∫_q^{4q}(u−q)du=9q²/2, both consequences follow.

**Scope.** Between k and 2k an angular sector must gain an exponentially large factor in probability even in each indicated shrinking cylinder. The variance lower bound is of the same order as that of a regularly varying log partition function proportional to −q^{α/2}. By itself it causes no contradiction with the two-sided logarithmic contact bounds.
