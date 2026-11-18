# Infinite-order flat contact with an integrated mass hypothesis

Consider the normalized exact Gaussian variance graph and the local positive potential S_R of round 34. Translate the contact to u=0. Assume v is continuous, nonnegative and U(r)=sup_{|u|≤r}v(u)=O_N(r^N) for every integer N. No Hölder bound is assumed. Write A(r)=∫_{−r}^r v(u)du. If v is positive in every neighborhood, choose any small attained level m₀>0, set m_j=2^(−j)m₀, and define the first radius r_j attaining maximum m_j. For each j set R=r_{j−1}. Then for every fixed 0<s<1,

    (1−r_j/R)^p ≤ C_s m_j^s log(C_s r_j/A(r_j)),                    (1)

where p is the absolute wedge exponent of round 31 and j is sufficiently large. Therefore for every fixed 0<s<1,

    ∑_j [m_j^s log(C_s r_j/A(r_j))]^(1/p) = ∞.                      (2)

In particular, for each 0<β<1 there are arbitrarily large j with A(r_j)<exp(−m_j^(−β)). More precisely the same statement holds for A(r_j)/r_j and any fixed multiplicative constants in the exponent. Hence a uniform lower bound of stretched-exponential form exp(−K m_j^(−β)), β<1, on the integrated first-hit masses forces a zero neighborhood. This criterion allows strictly non-Hölder moduli on the graph.

**Uniform central estimate.** At the fixed level m=m_j the core satisfies v≤2m. The entire jet of S_R at zero is given by its integrable off-diagonal jets, as in round 34. This identity requires only flat pointwise upper bounds and holds independently of a Hölder modulus. For q(u)=u²/v(u) on v>0 and any N>2,

    q≥max(u²/(2m),1/(K_N |u|^(N−2))) ≥c_N m^(−1+2/N).

The last inequality is obtained by taking powers 1−2/N and 2/N, respectively. Fix s<1 and choose N large with t=1−2/N>s. Interpolate q≥c_N m^(−t) with q≥1/(K_N|u|^(N−2)) to obtain q≥c m^(−s)|u|^(−a), a=(N−2)(1−s/t)>0. Thus q≥c' m^(−s) log(4R/|u|) uniformly for 0<|u|<R≤1. For n≤c''m^(−s) the factor (4R/|u|)^n in the normalized Taylor bound is absorbed by exp(q/8). The resulting coefficients satisfy R^n|[y^n]S_R|≤C exp(−c_N m^(−t)). The rest of the Taylor series is estimated on |y|≤R/32 using the fixed holomorphic bound on radius 3R/8, producing

    sup_{|y|≤R/32} S_R(y) ≤ C_s exp(−c_s m^(−s)).

By the established wedge chain, for |y|≤r_j one has S_R(y)≤C_s exp[−c_s(1−r_j/R)^p m^(−s)].

**Integrated testing without a modulus.** For each |u|≤r_j all the variance scales 0<t≤v(u) obey t≤m_j≪r_j² by flatness at the level. Among intervals [−r_j−u,r_j−u] for |u|≤r_j, the centered Gaussian N(0,t) has least mass at u=±r_j. Alternatively each window contains either [0,r_j] or [−r_j,0]. Therefore its probability is bounded below by Φ(1)−1/2 for large j. Fubini and the strictly positive lower bound of the Gaussian source near zero give

    ∫_{−r_j}^{r_j} S_R(y)dy ≥ c A(r_j).

Integration of the wedge upper bound now gives (1). The telescoping product ∏r_j/r_{j−1} vanishes, so ∑(1−r_j/r_{j−1}) diverges, proving (2). If A(r_j)/r_j≥exp(−K m_j^(−β)) eventually, choose s∈(β,1); the upper bound in (1) is summable after taking pth roots, a contradiction.

**An explicit regularity class.** Suppose in the neighborhood |v(u)−v(u')|≤K[log(e/|u−u'|)]^(−γ) with γ>1. At a first hitting point choose an inward interval of length ℓ=exp[−(2K/m_j)^(1/γ)]. It lies in the core for large j since r_j≫√m_j; throughout v≥m_j/2. Hence A(r_j)≥m_jℓ/2. This contradicts (2) and proves the zero-neighborhood assertion without local Hölder regularity. The parameter-free version is a lower bound on the area itself and requires no modulus hypothesis.

**Why continuity alone does not supply this area.** Take triangular positive bumps of height 2^(−j) centered at 1/j (j≥10), with halfwidth exp(−4^j), and take v=0 elsewhere. These disjoint bumps form a continuous nonnegative field, flat to every order in value at 0 and positive in each neighborhood. At its dyadic first hits r_j=1/j one has A(r_j)≤C 2^(−j) exp(−4^j). Indeed higher peaks are outside r_j and the sum of all later areas is smaller than the current one. Thus continuity and infinite-order pointwise flatness cannot by themselves yield the area bound used above. This example has no exactness assertion and does not refute flat-zero rigidity under the full Gaussian identity.

## Local output refilling and ultra-thin changes

The central smallness bound above holds on the **complex** disc |z|≤R/32 by the same summation of absolute Taylor coefficients. Take an arbitrary small R and m=U(R)>0. The same proof, with v≤m instead of 2m, gives by Cauchy's estimate

    sup_{|y|≤R/64}|S_R''(y)| ≤ C_s R^(−2)exp(−c_s U(R)^(−s)),   0<s<1.

Hence the truncated output measure has density g(y)+S_R''(y) on this interval, since H_v''=P_vδ₀−δ₀ in the distributional sense, also for v=0. At any point y with v(y)>0 the usual positive-variance density integral from the core equals this analytic representative pointwise: near y the variance has a positive lower bound, while means farther from y have a uniform Gaussian off-diagonal bound. Any zero-variance source is locally absent there.

For example, the isolated triangular bumps above cannot themselves produce exact output, although their areas are too small for (2) to exclude. At y=1/j set R=128/j, small for large j. The home bump contributes at most Cℓ_j/√m_j to the core output density there, integrating 1/√v over its triangular support (ℓ_j=exp(−4^j)). Every remaining included center 1/k has k≥j/129 for large j, and at distance at least c/j² from y. For large j the widths do not reduce this bound by more than one half. Since sup_{s>0}p_s(z)≤C/|z|, these other contributions are bounded by Cj²∑_{k≥j/129}ℓ_k=o(1). Zero-variance sources do not contribute at this positive point. This contradicts the required g(1/j)+o(1)→g(0)>0. This uses no endpoint regularity or analysis of an individual cusp.

There is a separate **analytic stability criterion**. Suppose a nonzero bilateral flat reference v₀ obeys the local Hölder condition of round 34, or the area criterion above. Let a continuous graph v agree with it in a fixed core except on a measurable E. If v,v₀>0 on E up to null sets and some δ>0 satisfies

    ∫_E ∑_{s∈{v(u),v₀(u)}} s^(−1/2) exp(δ²/(2s))du <∞,

then v cannot be exact. Indeed the output difference from E is analytic on |Im z|<δ by the complex Gaussian bound and dominated integration. If v were exact, its fixed-core output has an analytic density near the origin by off-diagonal analyticity of exterior input. The fixed-core output of v₀ would therefore also be analytic. Its fixed-core potential is analytic there because its second distributional derivative is the analytic output density minus g. Removing the annular inputs down to any smaller R preserves all uniform wedge bounds needed by round 34: the fixed analytic potential is bounded on a fixed disc and off-diagonal H has a common bound in the wedge discs. Thus the local proof applies to v₀, contradicting its nonzero bilateral flat germ.

This last result includes continuous graphs with no zero except the contact, no positive logarithmic Hölder exponent, and insufficient first-hit area: start with v₀(u)=exp(−exp(1/|u|)) (value zero at 0), and add positive triangular bumps at u=1/j of halfwidth ℓ_j=exp[−exp(3e^j)], so the peak value is 2^(−j). For large j, v≥v₀≥exp(−2e^j) on each bump and the stated exponential integral is bounded by a constant times ∑ℓ_j exp(C exp(2e^j))<∞. The first-hit radii at these levels equal 1/j for j large, since the bounded derivative of v₀ is dominated by each triangular slope. Moreover A(1/j)≤2j^(−1)exp(−e^j)+C∑_{k≥j}ℓ_k 2^(−k), which is smaller than exp(−(2^j)^β) for every fixed 0<β<1 and all large j. The field violates every log-power modulus at these bump widths and is flatter than every power in value at zero. This is an exclusion of a specified non-Hölder class, not a proof for every continuous bilateral flat graph.
