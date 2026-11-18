# Directional drift integrability forces stretched-exponential inverse variance at a contact

At a highest point a write u=x−a, v=H−f, T=1/(2v), k=u/v on v>0. Restrict to a fixed small bounded interval U around a. The finite positive measure
$$
 d\nu=p_{v(x)}(a-x)g(x)\,1_{U\cap\{v>0\}}dx
$$
has finite mass by Fatou in the output identity. Suppose ∫exp(d|k|^p)dν<∞ for some d>0 and 1≤p<∞. Then for some c>0
$$
 \int \exp(c T^{\,p/(p+2)})\,d\nu<\infty .
$$
If k is essentially bounded, the conclusion holds with exponent one. Neither isolation of the highest point nor a pointwise lower variance bound on both wings is required.

To justify local analyticity, the entire output near a minus source positions off U is analytic, using spatial separation for small variances. The positive-variance source within U has a continuous density in a small interval around a: at a+y it is exactly ∫exp(ky−Ty²)dν and the hypothesis dominates by an integrable exponential in |k|. The term remaining from U is g times the zero-set indicator almost everywhere near a. It therefore has a continuous representative, which divided by g takes values in {0,1} almost everywhere. On a connected smaller interval it is identically 0 or 1 by continuity. In both cases the positive-variance density is analytic near a, after subtracting either zero or g. Consequently its even part
$$
 S(r)=\int e^{-Tr}\cosh(k\sqrt r)\,d\nu
$$
extends analytically as a function of r near zero (its even power series in y has r=y²). Its derivatives at all sufficiently small positive r are bounded by C Aⁿ n!.

For r>0 one may differentiate under the integral: an exponential in −T controls each finite power of T near that fixed r and the assumed exponential controls the drift powers for small enough r. The nonnegative power series of cosh gives 0≤∂ᵣʲcosh(k√r)≤|k|^(2j)cosh(k√r), by termwise comparison. The absolute difference between (−1)ⁿS⁽ⁿ⁾(r) and the leading integral ∫Tⁿe^(−Tr)cosh(k√r)dν is at most
$$
 \int [(T+k^2)^n-T^n]\,e^{-Tr}\cosh(k\sqrt r)\,d\nu .
$$
Where T≥4nk² the bracket is ≤(exp(1/4)−1)Tⁿ. On the complement it is ≤(5n)ⁿ|k|^(2n). Absorb the first part and bound the remaining cosh by exp(|k|√r). The exponential drift moment, after reserving half its exponent, bounds the integral of |k|^(2n)exp(|k|√r) by C Bⁿ n^(2n/p) uniformly for small positive r (maximize t^(2n)exp(−d t^p/2)). For bounded k the bound is C Bⁿ. Fatou at zero gives
$$
 \int T^n d\nu \le C B_1^n n^{(1+2/p)n}.
$$
Markov's inequality with n of order R^(p/(p+2)) and a sufficiently small proportionality constant gives an exponential tail for T>R at the claimed power, hence the asserted moment with a smaller exponent coefficient. The bounded-drift case uses n of order R.

**One-sided implication.** If on a right wing c₁u^α≤v(a+u)≤c₂u^α for 0<α≤1, then for every d>0 and every sufficiently small δ the left integral
$$
 \int_{-\delta<u<0,\ v(a+u)>0}
 v(a+u)^{-1/2}\exp[(d|u|-u^2/2)/v(a+u)]\,du
$$
must diverge. On the right wing the drift is bounded and the input is locally bounded above and below. If the left integral were finite, the p=1 theorem would force a stretched exponential inverse moment under ν; this diverges on the right because u²/v→0 and v≤c₂u^α. This identifies a quantitative necessary contribution from the opposite wing even if that wing has no pointwise contact order.
