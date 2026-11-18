# Audit of the Gaussian-multiplier lower bound, all contact orders 2m≥6

There is a strict interior entropy excess for every integer m≥3. The phase, parity and torus arguments below use N even; for the explicit test take N a multiple of 2m(m−1).

Write d=m−1, a=(cN)^(−1/(2d)), S=Nv(az)/a²=z^(2m)(1+o(1)) locally. For λS<1 the Hermite sign is positive if S(1−λS)/z²≤1/4. The supremum of S^((m−1)/m)(1−λS) is m/(2m−1) times [(m−1)/((2m−1)λ)]^((m−1)/m). Thus choose λ strictly above

    λ_m=(m−1)/(2m−1) [4m/(2m−1)]^(m/(m−1)).

The strict inequality handles the uniform 1+o(1) for the actual graph on a fixed sufficiently small core. All off-diagonal Gaussian coefficients after multiplication by exp(λNW²/2) are then nonnegative.

For constant source and pure jet the kth diagonal group is

    D_{N,k}=binom(2mk,2k)(2k−1)!! N^(−k)
            (λN/2)^(N/2−d k)/(N/2−d k)!.

Analyticity of the original density forces its complete weighted coefficient to have exponent at most f(0)=(1+log λ)/2. Writing ν=k/N and 0 log 0=0, the normalized log D is f(ν)+O(log(N+2)/N) uniformly on 0≤ν≤1/(2d), where

    f(ν)=ν[2m log m−2d log d+log(2ν)−1]+(1/2−dν)[1+log(λ/(1−2dν))].

At λ=λ_m and ν=1/(m d), simplification gives the excess f(ν)−f(0)=A(1/m)/d with

    A(x)=−(3−2x)log(1−x)+(2−x)log(1−x/2)−1
         −(1−x)(1−2x)/(2x)log(1−2x).

I independently expanded the geometric logarithm series: a₂=1/12, a₃=1/24, and the coefficients a_j at j≥4 are

    [(j−3)(1−2^(j−1)/(j+1))+2^(1−j)]/[j(j−1)].

These are negative and bounded in absolute value by 2^(j−1)/(j(j−1)). For x≤1/4 the summed negative tail divided by x² is ≤(2/3)x²/(1−2x)≤1/12, leaving a positive cubic term. For m=3 one obtains 3A(1/3)=log(84375/4096)−3>0 via e<49/18. Thus each fixed m has a λ>λ_m close enough to keep a strict interior excess. Constants need not be uniform in m. Also f'(ν)=2m log m−2d log d+log(2ν)−d log[λ/(1−2dν)], so f' tends to minus infinity at the right endpoint and the global maximum is strictly above that endpoint as well. Hence a compact interior range and a positive gap above both endpoint neighborhoods exist. This works regardless of which endpoint is larger.

The weighted pure sextic note's density-jet identity extends to v₀=c u^(2m)H(u)² by the same finite Taylor argument (near variable u=w+|w|^m t, producing O(|w|^((m−1)(L+1))) and the arbitrary flat kernel difference contributes only its off-diagonal coefficients. For a general analytic H the kth diagonal group is the corresponding positive Gaussian moment times the bivariate coefficient with (W+T)^(2mk) H(a(W+T))^(2k)g(a(W+T)). The circular majorant used to bound its endpoints remains valid, replacing 3/2 by m/(m−1) in the binomial bound. In the middle the positive saddle radii are R²=(1−2dν)/λ, U=R/d. The log-radius Hessian is the sum of diag(2λR²,0) and 2mνRU/(R+U)² times [[1,−1],[-1,1]], positive definite uniformly on compact interior ranges. The two real saddles have equal parity signs; the small real analytic 2ν log H perturbation preserves both. The two local radial deformations of the Cauchy torus in disjoint angular patches are homologous to the original, with uniformly bounded Jacobians and exponentially small transitions. Thus every middle group is eventually positive uniformly and has log-rate f(ν)+O(a)+o(1); a maximizing interior group dominates the absolute contributions from the endpoint neighborhoods. No off-diagonal saddle lower bound is required.

For the pure sextic independently of the series argument, at λ²=121/54, N=20l, k=3l, the normalized block ratio has (1/l)th power approaching K/e⁹, K=1029455660473245/116101021696. The rational bound K/(68/25)^9>1 and e<68/25 give a distinct strict lower bound. Code and output in r043_weighted_hermite_independent_check.py/.json verify this certificate, the Hermite sign margin and the stated real saddle margin. The pure sextic one-dimensional saddle in the source also passes the circular-contour audit; it is unnecessary once the entropy excess is available.

**Divergent jets in the final manuscript.** I also checked r043-all-orders-weighted-jets.md after it became available. Its degree-N replacement H_N suffices because the kth group's largest requested residual degree is N+2k−2mk≤N. With s=1/(2d) and Rmax=m/(d√λ), τ<(ec)^s/Rmax indeed gives an absolute-majorant bound H_N(aZ)=1+O(a) on a strictly larger fixed Z-disc. Stirling gives (j!)^s≤C j^(3s/2)(j/e)^(s(j−1)), so after extracting the first factor a the remaining sum is dominated by a convergent geometric series times a polynomial in j, uniformly in N. The same bound in complex modulus yields O(a) for each fixed derivative on a smaller disc. Consequently all perturbed Cauchy-torus estimates remain uniform without analyticity of the infinite formal series. In particular λ=λ_m exp(m g_m/2) retains test excess g_m/2 while allowing equality at the written type threshold (ec)^s d√λ_m/m. The stated open extension up to type (ec)^s(d/m)√λ_m exp(mg_m/2) is valid by choosing λ<λ_m exp(mg_m). No claim for arbitrary divergent jets follows.
