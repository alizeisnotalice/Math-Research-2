# Distance-power variance at an arbitrary minimum set

Let \(Z\) be a nonempty closed subset of the circle \(\mathbb T_L\), let \(v\ge0\) be bounded, and suppose that on a neighborhood of \(Z\), almost everywhere,
\[
 c\,d(x,Z)^\alpha\le v(x)\le C\,d(x,Z)^\alpha
\]
for \(c,C>0\). We study the stationary condition that \(X+\sqrt{v(X)}G\pmod L\) is Haar when \(X\) is Haar and \(G\) is an independent standard Gaussian. The graph is assumed nonconstant. One may take \(Z\) to be the zero set of a continuous graph after subtraction of its minimum. In the measurable formulation, it suffices that at the particular complementary gap endpoint used below the indicated comparisons hold and that the graph has essential infimum zero.

A positive measure of zero variances is impossible. The output on that set already receives its full Haar measure from the deterministic zero-variance sources, and any positive-variance source has strictly positive wrapped transition density. Thus \(Z\) has measure zero and its complement has a gap. For \(\alpha\ge2\), the upper comparison gives a one-sided \(O(u^2)\) variance bound at every endpoint, excluded by the R65 measurable quadratic-average theorem. The following arguments concern \(0<\alpha<2\).

## All geometries are excluded for \(1<\alpha<2\)

Only one complementary gap is used. Put its endpoint at 0 in local circular coordinates, with gap on the positive side, and take \(r_0>0\) small enough that \(v(u)\) lies between \(c u^\alpha\) and \(C u^\alpha\) almost everywhere on \((0,r_0)\). No assumption on the opposite-side minima is used. For \(y>0\) small, split the output density into the principal Gaussian from \(0<u<r_0\), denoted \(q_+(y)\), the principal Gaussian from \(-r_0<u<0\), denoted \(q_-(y)\), and the remainder \(q_{
m ext}(y)\). Normalize the source as Lebesgue measure so the total target density is 1.

The local gap density is continuous on the open gap. For observation points there, small source variances are separated from the observation or else bounded positively below by the distance comparison. The remainder \(q_{\rm ext}\) is continuously differentiable near zero, with uniformly bounded derivative by separated Gaussian derivative estimates and summability over wraps. The positive wing density extends continuously from the right to 0 by \(p_{v(u)}(y-u)\le C'u^{-\alpha/2}\), which is integrable. At the negative wing, the contribution at 0 is finite by Fatou applied to the exactly uniform wrapped outputs after adding a constant \(\zeta>0\) to all variances. For \(y\ge0\), \(p_{v(u)}(y+|u|)\le p_{v(u)}(|u|)\) when \(u<0\) and \(v(u)>0\). This gives continuity from the right by domination and \(q_-(y)\le q_-(0)\). Consequently the sum of the three boundary values equals 1.

Let \(\beta=2/\alpha-1\in(0,1)\). The integrand of \(q_+(y)-q_+(0)\) is nonpositive when \(0<u\le y/2\). On the source interval \(u\in[y^{2/\alpha},2y^{2/\alpha}]\), contained there for small \(y\), we have \(v(u)\asymp y^2\). Thus \(p_{v(u)}(u)\gtrsim y^{-1}\) and the quotient \(p_{v(u)}(y-u)/p_{v(u)}(u)=e^{-(y^2-2yu)/(2v(u))}\) is bounded below 1 by a fixed gap. These sources yield a decrease of at least \(c_1y^\beta\). The only possible positive part is \(u>y/2\). Its absolute integral is bounded by
\[
 C_1y\int_{y/2}^{r_0}u^{1-3\alpha/2}\,du
 =\begin{cases}O(y),&\alpha<4/3,\\O(y\log(1/y)),&\alpha=4/3,\\O(y^{3-3\alpha/2}),&\alpha>4/3.\end{cases}
\]
Indeed the spatial derivative at arguments \(s-u\), \(0\le s\le y\), is bounded by \(C u^{1-3\alpha/2}\), since \(v(u)\gtrsim u^\alpha\) and \(|s-u|\le3u\). Every displayed bound is \(o(y^\beta)\); for the last use \(3-3\alpha/2-\beta=(2-\alpha)(3\alpha-2)/(2\alpha)>0\). Hence \(q_+(y)-q_+(0)\le-c_1y^\beta+o(y^\beta)\). The opposite wing is nonincreasing and the remainder changes by \(O(y)=o(y^\beta)\). The sum contradicts stationarity. Positive bounded source density near the contact and a differentiable target suffice, so the same local proof also applies to the ordinary Gaussian source and target on the line.

## Even derivatives force opposite-wing crowding for every \(0<\alpha<2\)

At each complementary gap endpoint, choose coordinates as before. Fix any \(0<\theta<1\). There is \(a_\theta>0\) such that for all sufficiently small \(r>0\), the opposite wing must contain another point of \(Z\) with distance from the endpoint in \([a_\theta r,r^\theta]\).

To prove this, choose \(n\ge1\) so \(2n(1-\theta)>2/\alpha-1\). Let \(\delta_n>0\) be small enough that \((-1)^n He_{2n}\) is strictly positive on \([-\delta_n,\delta_n]\). Add \(r^2\) to the variance, preserving stationarity, and differentiate the smooth output at target 0 to order \(2n\). For \(|u|\le\delta_n r\), the principal Gaussian term has the sign \((-1)^n\). Its integral after multiplication by \((-1)^n\), restricted to the inward wing \(u\asymp r^{2/\alpha}\), is at least \(c_n r^{-2n-1+2/\alpha}\), since the variance with smoothing is comparable to \(r^2\) and the source length is comparable to \(r^{2/\alpha}\). The rest of the very near part has the same sign.

Now suppose the opposite wing has no zero at distances in \([\delta_n r/2,r^\theta]\). For \(\delta_n r\le|u|\le r^\theta/2\) on that wing, the nearest possible minimum is at distance at least \(|u|/2\). The principal kernel argument \(|u|/\sqrt{v(u)+r^2}\) is then bounded by \(C |u|^{1-\alpha/2}\), so on these shrinking annuli the even derivative has the same sign. On the inward wing the same conclusion follows from \(d(u,Z)=u\) for sufficiently small \(u>0\). The only possible adverse principal contribution is thus at \(|u|>r^\theta/2\). The uniform separated Gaussian derivative bound \(|p_s^{(2n)}(-u)|\le C_n|u|^{-2n-1}\) gives a total for it of at most \(C_n r^{-2n\theta}\); nonprincipal wraps give \(O(1)\). Since \(2n(1-\theta)>2/\alpha-1\), neither can cancel the inward contribution, a contradiction. Take \(a_\theta=\delta_n/2\).

A more quantitative condition is available from the same argument. For each fixed \(n\) with \(2n>2/\alpha-1\), the adverse principal contribution in the opposite near wing can occur only where \(|u|\ge\delta_n r\) and \(v(u)\le\delta_n^{-2}u^2\). On the dyadic shell \(s\le|u|\le2s\), these sources lie within distance \(C_n s^{2/\alpha}\) of \(Z\). Write \(T_n(s)\) for the ordinary length of this tube within the shell on the opposite wing. With a sufficiently small fixed upper cutoff \(r_0\), stationarity requires
\[
\sum_{j\ge0:\, 2^j\delta_n r<r_0}2^{-j(2n+1)}T_n(2^j\delta_n r)
\ge c'_n r^{2/\alpha}
\]
for all sufficiently small \(r\). The bound follows by using the separated derivative upper bound on each adverse shell and absorbing the bounded outer and wrapped derivatives into the divergent positive term. An opposite wing with power-sized logarithmic holes in its zero distances along a sequence is ruled out even for \(\alpha\le1\). The estimate allows a sufficiently crowded minimum set and gives no hypothesis-free exclusion in that remaining range.

## Checks

`gap_check.py` verifies the exponent comparison and example choices of even derivative orders exactly in SymPy. It also runs 45- and 70-digit mpmath quadrature on a pure inward wing \(v(u)=u^\alpha\) for \(0<u<1/4\). For \(1<\alpha<2\), the synthetic normalized density difference converges to
\[
\frac{2^{1/2-1/\alpha}}{\alpha\sqrt{2\pi}}
   \Gamma\left(-\tfrac12(2/\alpha-1)\right),
\]
obtained by scaling \(u=y^{2/\alpha}z\) in the variance-only difference. The expression is negative. These numerical checks calibrate the one-sided mechanism and are not a construction of any exact stationary graph.
