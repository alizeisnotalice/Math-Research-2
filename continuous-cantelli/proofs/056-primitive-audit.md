# Independent audit of the spatial potential primitive

The following finite-order conclusion is valid. If 0≤v(x₀+u)≤C|u|^α, α>2, on a bilateral neighborhood in an exact graph, then for some single polynomial Q,

    ∫₀ʰ p₁(x₀+u)v(x₀+u)du = Q(h)+O(|h|^(3α/2)),

simultaneously as h↓0 and h↑0. Its coefficients of order <α+1 vanish. Terms of degree exactly 3α/2, when integral, are not identified by this error alone.

Indeed use L(y)=∫g(u)∫₀^{v(u)}p_t(y−u)dt du, analytic at zero by exactness and the Gaussian call identity. For I between 0 and h and |u|≤8|h|, every t≤v(u) satisfies t≤U=C'|h|^α. Away from the endpoints where indicators are unambiguous, the difference between P(u+√tZ∈I) and 1_I(u) is bounded by the sum of the two tails crossing the nearest endpoint, hence by a constant times exp(−dist(u,∂I)²/(2t)). Integration over u is O(√t), and integration over time yields O(U^(3/2)). Orientation is reversed on both sides of the identity when h<0. Choosing the split 8|h| avoids any contour radius issue.

For |u|>8|h| and y∈I the mean-specific potential has a Cauchy expansion on |z|≤|u|/4 with bound C exp(−c|u|^(2−α)) near zero. At every fixed finite order its remainder is absolutely integrable as O(|h|^(k+1)), and the coefficients integrated over the omitted inner means are negligible to any finite order. The far means have an analytic exterior contribution by fixed separation. Thus the integrated outer contribution has a single two-sided polynomial through the chosen order, and subtracting from the analytic primitive of L gives the claim. On an interval [r,r(1+δ)] with 0<δ≤1 subtract the primitives at the endpoints. After division by the leading average scale δr^(α+1), the remainder is O(r^(α/2−1)/δ). This permits all windows δ≫r^(α/2−1).

**Signed residuals.** Let P≥0 be a smooth O(u⁴) reference, |v−P|≤C|u|^α and v+P≤C|u|^p with p>2. There is a single polynomial Q_res such that

    ∫₀ʰ g(u)(v(u)−P(u))du = Q_res(h)+O(|h|^(α+p/2)).

In the inner region the time interval between P(u) and v(u) has length O(|h|^α). At every intermediate variance its boundary leakage is ≤C exp(−dist(u,∂I)²/(C|h|^p)). Integrating over u now costs O(|h|^(p/2)), proving the displayed error. All off-diagonal Taylor integrals for the two variances are treated as above. For the reference potential primitive, first replace P by its Taylor polynomial of a sufficiently high fixed degree. If this polynomial is nonzero its first nonzero coefficient is positive and has even degree, so after shrinking the core it is a nonnegative finite-type reference. Its local output has two-sided Peano expansions through every fixed order by the standard near-diagonal change of variables at its first nonzero order; its potential therefore does too by integration twice against the density difference from g. If the Taylor polynomial vanishes, use zero instead. The replacement affects ∫g(v−P) and the hypotheses by errors of arbitrarily high prescribed order. Thus the signed conclusion applies also when the original reference has zeros accumulating at the contact. The relative window floor is r^(p/2−1). An averaged claim at a comparable or smaller window needs additional estimates; it does not follow just from the primitive error.
