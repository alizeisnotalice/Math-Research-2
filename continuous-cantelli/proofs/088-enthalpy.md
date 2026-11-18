# Early enthalpy uniqueness and the first unequal-tail branch

Let \(g=p_1\), let \(f\) be bounded continuous and nonnegative, and assume the exact Gaussian identity
\[
 \int e^{\lambda x-\lambda^2 f(x)/2}g(x)\,dx=1
                  \qquad(\lambda\in\mathbb R).
\]
Throughout the early band \(0<t<T\), assume the hypotheses of round 84 have yielded connected past intervals \(I_t=\{f\le t\}\) containing zero and the nonnegative surviving density \(q_t\) of \(Q_t\) supported on \(I_t\). In particular \(0\le q_t\le p_t\). This is available for connected early intervals from round 84, or from the record-birth exclusion hypotheses of round 87. We do not assume the sign of \(Q_t\) above \(T\).

## 1. Energy proof of early enthalpy uniqueness

Put \(\mathcal H_t(x)=q_t(x)+g(x){\bf1}_{I_t}(x)\). On \(I_t\), \(\mathcal H_t\ge g\) and \(q_t=\mathcal H_t-g\), whereas outside \(I_t\) both are zero. Thus
\[
 q_t=(\mathcal H_t-g)_+,\qquad
 \partial_t\mathcal H_t=\tfrac12\partial_x^2 q_t .                \tag{1}
\]
The distributional equation follows by adding the activated graph density to \(Q\), cancelling the graph sink in the heat equation. The initial measure is \(\delta_0\). More explicitly,
\[
 \mathcal H_t=P_t\delta_0-\int_{f(x)\le t}
                  (P_{t-f(x)}-I)\delta_x\,g(x)dx.              \tag{2}
\]
The zero-age difference vanishes. Bounded continuous spatial tests and dominated convergence give weak continuity in \(t\), including at zero, from (2). The total enthalpy is one.

To compare two early enthalpies with the same \(g\) and initial \(\delta_0\), let
\[
 d_t=\mathcal H^{(1)}_t-\mathcal H^{(2)}_t,\qquad
 b_t=q^{(1)}_t-q^{(2)}_t,\qquad
 D(t,x)=\int_{-\infty}^x d_t(y)dy .
\]
On \(0\le t\le T'<T\), there is a common compact spatial support for the two enthalpies and the two \(q\)'s, since the early sublevel intervals are bounded and increase. The total masses agree, so \(D\) also has compact support and \(|D|\le2\). It is continuous in \(t\) with values in \(L^2(dx)\): for positive \(t\), the measures in (2) have no atoms, so their distribution functions converge pointwise and are uniformly bounded; at \(t=0\), the two limits are the same \(\delta_0\), so their difference tends to zero for every \(x\ne0\).

Equation (1) gives \(\partial_t D=\tfrac12\partial_x b\) distributionally. The integration ambiguity is spatially constant and vanishes by the common compact support. For \(0<s<t\le T'\), mollify in space with a compactly supported smooth approximate identity. On \([s,t]\), \(d\) and \(b\) are bounded and square-integrable by \(q_i\le p_u\), where \(u\) is the time variable. The mollified equation has an \(L^2\)-valued time derivative and gives
\[
 \frac12\|D_\varepsilon(t)\|_2^2-\frac12\|D_\varepsilon(s)\|_2^2
 =-\frac12\int_s^t\!\int
                     (d_u)_\varepsilon(x)(b_u)_\varepsilon(x)\,dxdu.
\]
This follows by pairing with \(D_\varepsilon\) and integrating by parts in space; a further time mollification gives the same identity if desired. At the two endpoints use strong spatial \(L^2\) mollification. In the last integral both factors approximate strongly in spacetime \(L^2\). By the pointwise monotonicity of \(h\mapsto(h-g(x))_+\), \(d_ub_u\ge0\). Hence \(\|D(t)\|_2\le\|D(s)\|_2\). Sending \(s\downarrow0\) gives \(D(t)=0\), and therefore both enthalpies agree at each early time.

Apply this to reflection about zero, using evenness of \(g\) and invariance of the initial point source. Thus \(\mathcal H_t\), \(q_t\) and \(I_t\) are symmetric for every \(0<t<T\). For the interval, use that \(\mathcal H_t\ge g>0\) almost everywhere inside and vanishes outside. Consequently \(f(r)=f(-r)\) whenever at least one of these two values is strictly less than \(T\). This last statement deliberately leaves the equality level \(T\) itself open to a one-sided plateau.

## 2. Unequal limits force the compensating branch at the critical level

Assume distinct actual tail limits, both greater than one, and reflect if necessary to arrange
\[
 H_+=\lim_{r\to\infty}f(r)>H_-=\lim_{r\to\infty}f(-r)>1.
\]
Take the early band above to extend through every \(t<T_*:=H_--1\). Write \(f_+(r)=f(r)\), \(f_-(r)=f(-r)\) for \(r\ge0\). The exact odd identity is, for \(\lambda>0\),
\[
 \int_0^\infty \sinh(\lambda r)
         \big[e^{-\lambda^2 f_+(r)/2}
                 -e^{-\lambda^2 f_-(r)/2}\big]g(r)\,dr=0.        \tag{3}
\]
Its tail integrand is negative: the higher-tail wing is the one with the larger \(f\), which produces the smaller exponential weight. Consequently \(S=\{r>0:f_+(r)<f_-(r)\}\) is a nonempty bounded open set and supplies the positive compensating contribution. The exact critical level is
\[
 \inf_{r\in S} f_+(r)=T_*,
 \qquad
 \inf_{\{r:f_+(r)\ne f_-(r)\}}\min(f_+(r),f_-(r))=T_* .          \tag{4}
\]
Early symmetry gives both lower bounds. To prove the first upper bound, pick \(\delta>0\) sufficiently smaller than \(H_+-H_-\). On a fixed far ray, \(f_-(r)\le H_-+\delta\) and \(f_+(r)-f_-(r)\ge c>0\). The negative contribution from that ray has magnitude at least
\(c' e^{-(H_--1+\delta)\lambda^2/2}\) for large \(\lambda\), since the Gaussian-weighted integral of \(\sinh(\lambda r)\) on the ray is at least \(c''e^{\lambda^2/2}\). Meanwhile the total positive contribution, supported in \([0,R]\) for some finite \(R\), is at most \(C e^{-\lambda^2\inf_S f_+/2+R\lambda}\). If the infimum exceeded \(T_*\), choose \(\delta\) smaller than the gap, contradicting (3). This proves (4).

In particular arbitrarily near the level \(T_*\) there are points at which the higher-tail \(f\)-wing is lower than the other. If speaking in terms of the original conditional variance \(a=C-f\), both inequalities are reversed. Statements (3)–(4) do not assert a single sign change, or that an infimum of mismatch radii is attained as an unequal point.

## 3. First radius and time-level plateaus under monotone wings

If both outward profiles are nondecreasing on the entire half-line, let \(a_*=\lim_{t\uparrow T_*}a(t)\) be the common early endpoint. Let \(r_\pm=\sup\{r:f_\pm(r)=T_*\}\); both plateau intervals are \([a_*,r_\pm]\), allowing zero length. When \(H_+>H_-\), one must have
\[
 r_+\ge r_- .
\]
Indeed if \(r_+<r_-\), no compensating point \(f_+<f_-\) can be at \(r\le r_-\), since there the lower-tail wing is at or below \(T_*\), and the two wings agree below that level. Beyond \(r_-\), \(f_+(r)\ge f_+(r_-)>T_*\), contradicting (4).

The first mismatched radius, as an infimum, is \(r_-\). If \(r_+>r_-\), the whole intervening open strip \(r_-<r<r_+\) has \(f_+<f_-\). If the two plateau ends coincide, (4) forces compensating points approaching their common end from the right, although no full neighborhood sign is asserted. Both values at the first radius equal \(T_*\) by continuity. A common plateau may make the radius strictly larger than the entry endpoint \(a_*\). Without global wing monotonicity, (4) remains true but this radius argument is unavailable.

The one-sided early support argument of round 87 actually excludes the high-wing \(T_*\) plateau here. The right boundary support continues to all \(t<H_+-1\), whereas \(T_*<H_+-1\). If \(r_+>a_*\), the right endpoint \(B(t)\) jumps from \(a_*\) to \(r_+\) at \(t=T_*\), and the new extreme has an entire positive-length level interval on its inner side. The reflected-kernel bound of round 87, applied with this one-sided time band, would make the local mass of that level exponentially small as the excess parameter tends to zero. This is impossible. Thus \(r_+=a_*\), and then \(r_-=a_*\) as well. The actual first mismatch radius is \(a_*\). If global wing monotonicity is omitted but the right record at \(T_*\) satisfies the local regularity required by the same jump theorem, its exclusion still gives \(B(T_*)=a_*\). The compact compensating set \(S\) can then approach the lowest level \(T_*\) only at \(a_*\), yielding the same first radius.

There is also a quantitative restriction on the lower wing:
\[
 \limsup_{r\to\infty} r[H_- - f_-(r)]\le 2r_+.                   \tag{5}
\]
For large \(\lambda\), use the negative tail contribution from \(r\in[\lambda-1,\lambda]\) in (3). Monotonicity gives \(f_-(r)\le f_-(\lambda)\); the Gaussian shift gives the lower bound \(c e^{-(f_-(\lambda)-1)\lambda^2/2}\). For any \(a>r_+\), the positive compact branch on \(r\le a\) is bounded by \(C e^{-T_*\lambda^2/2+a\lambda}\); on \(a<r\le R\) there is a fixed stronger time exponent since \(f_+(a)>T_*\). Comparing and sending \(a\downarrow r_+\) proves (5).

## 4. Symmetric early endpoint tail budget and local checks

Early symmetry gives the common entry radius \(a_*\) in the round-84 quadratic tail budget. Thus, for monotone wings in the orientation above,
\[
 H_+-H_-\ \ge\
 2a_*\left(\frac{g(a_*)}{\bar\Phi(a_*)}-a_*\right)
   + \bigl[H_+-E(f_+(|X|)\mid X>a_*)\bigr]
   + \bigl[H_--E(f_-(|X|)\mid X>a_*)\bigr].
\]
Here both smaller conditional mean terms are on the same Gaussian half-line, and \(a_*\) comes from the interval approached through times \(t<T_*\).

The companion Python 3.13.18 script uses SymPy 1.14.0 to check the constitutive product in the mixed regime, the Gaussian shift identity, and the branch crossing factorization of an explicit synthetic example. mpmath 1.3.0 at 80 digits evaluates the model odd integrals for \(\lambda=8,\ldots,128\). After removing the common quadratic exponent, their logarithms divided by \(\lambda\) approach 1 for the compact positive piece and \(1/2\) for the negative tail piece, as predicted by the respective onset radius and lower-tail deficit. The two values at \(\lambda=128\) are approximately \(0.9134697950\) and \(0.4956000797\); the model serves to check scales rather than satisfy the full identity.
