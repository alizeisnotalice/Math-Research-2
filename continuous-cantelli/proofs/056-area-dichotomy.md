# Weak contact tangents and local averages under a superquadratic envelope

Translate an attained contact to \(u=0\), and write \(g(u)=p_1(x_0+u)\).
Suppose the exact output law after shifting by \(x_0\) has Gaussian density \(P(y)\), and \(v\ge0\) is bounded continuous with
\[
 v(u)=O(|u|^\beta),\qquad \beta>2 ,
\]
near zero. No derivative of \(v\) and no convergence of normalized values is assumed. For \(s>0\), let \(p_s\) denote the Gaussian density of variance \(s\). Define
\[
 K_\ell(z)=\int_0^\ell p_s(z)\,ds,\qquad
 S(y)=\int g(u)K_{v(u)}(y-u)\,du .
\]
As a locally integrable function, \(S''=2(P-g)\) distributionally: for a smooth compactly supported test, integrate the heat equation in \(s\) from zero to \(v(u)\), then apply Fubini. It follows that \(S\) has an analytic representative. Write \(s_j=[y^j]S(y)\).

For every integer \(j\ge0\), define the finite off-diagonal jet
\[
 a_j=\frac1{j!}\int_{u\ne0}g(u)\,\partial_y^j K_{v(u)}(-u)\,du
\]
and let \(b_j=s_j-a_j\). Finiteness follows from
\[
 \sup_{|y|\le |u|/2,\ 0<s\le v(u)}|\partial_y^j p_s(y-u)|
       \le C_j |u|^{-j-1}\exp(-c |u|^{2-\beta})         \tag{1}
\]
near zero. Indeed, with \(z=|y-u|/\sqrt s\), the normalized Gaussian derivative is a polynomial in \(z\) times \(e^{-z^2/2}\); extract \(e^{-z^2/4}\) and use \(|y-u|\ge |u|/2\). Away from zero the usual Gaussian off-diagonal derivative bound and bounded variance suffice.

## 1. Dichotomy and uniqueness of the weak tangent

Fix \(L>0\), choose any \(C>2L\), and split \(S=S_{\mathrm{in},r}+S_{\mathrm{out},r}\) at the source radius \(Cr\). Taylor expansion in \(y=r z\) using (1) shows, for every fixed \(N\), uniformly on \(|z|\le L\),
\[
 S_{\mathrm{out},r}(rz)=\sum_{j=0}^N a_j r^j z^j+O_{N,L,C}(r^{N+1}).
\]
To see that the coefficients do not depend on \(C,r\), extend the off-diagonal jet integrals to the full line. The missing core integrals are smaller than every power of \(r\) by (1). Thus the analytic representative of the inner potential satisfies
\[
 S_{\mathrm{in},r}(rz)=\sum_{j=0}^N b_j r^j z^j+O_{N,L,C}(r^{N+1}).   \tag{2}
\]
This is an equality of analytic representatives on the stated output interval, which suffices for all integrals below.

On testing (2) against any \(\psi\in C_c^\infty((-L,L))\), the inner term is exactly
\[
 \int_{|u|\le C} g(ru)\int_0^{v(ru)}
                 (P_{s/r^2}\psi)(u)\,ds\,du ,                    \tag{3}
\]
where \(P_a\) is the Gaussian smoothing operator. It is \(O(r^\beta)\) by the envelope. Testing a possible first nonzero coefficient in (2) proves \(b_j=0\) for all integer \(j<\beta\).

If \(m\) is the first nonzero index of \(b_j\), then for every finite \(L\)
\[
 \int_{-L}^L v(ru)\,du = O_L(r^m).                               \tag{4}
\]
For this, choose a nonnegative test at least 1 on an enlarged interval around \([-L,L]\), and choose a proportionally larger \(C\) for that test. Since all smoothing variances \(s/r^2\) in (3) are \(O(r^{\beta-2})\), its smoothed values are at least \(1/2\) on \([-L,L]\) for small \(r\). The positive lower bound for \(g\) and (2) give (4). Returning to any fixed smooth test, \(|P_{s/r^2}\psi-\psi|\le (s/2r^2)\|\psi''\|_\infty\). In (3), divide by \(r^m\), use (4), and replace \(g(ru)\) by \(g(0)\). This gives
\[
 \int \psi(u)\,\frac{v(ru)}{r^m}\,du
   = \frac{b_m}{g(0)}\int u^m\psi(u)\,du
           +O_\psi(r+r^{\beta-2}).                              \tag{5}
\]
In particular \(m\) is even and \(b_m>0\), because the limiting distribution is nonnegative on both half-lines and nonzero. The limit holds for all \(r\downarrow0\). Consequently any locally \(L^1\) or uniform value tangent at an exponent \(\alpha>2\) is zero if \(\alpha<m\), is this single monomial if \(\alpha=m\), and cannot exist as a locally finite limit if \(\alpha>m\).

If all \(b_j\) vanish, the same nonnegative test in (3) gives, for every \(N\),
\[
 \int_{-r}^r v(u)\,du=O_N(r^N).                                 \tag{6}
\]
This alternative is infinite-order flatness of area. No pointwise infinite-order flatness follows from (6) alone.

The convergence assumption on value tangents has a role away from a same-order envelope. Along a sequence, mere pointwise convergence under a weaker \(\beta\)-envelope may miss mass concentrated on shrinking sets. At \(\alpha=m\), Fatou alone yields that any pointwise limit is at most \((b_m/g(0))u^m\) almost everywhere. If the convergence instead holds as \(r\downarrow0\) without restricting \(r\) to a sequence, evaluation at \(u=1,-1\) supplies a same-order \(O(|x|^\alpha)\) envelope, so pointwise convergence implies the required local dominated convergence.

The argument still applies for envelopes \(v=o(u^2/\log(1/|u|))\). In that case (1) has a superpolynomial upper bound in \(|u|\) in place of the displayed exponential; \(s/r^2=o(1)\) uniformly in the core and the rate term \(r^{\beta-2}\) in (5) becomes \(O(\sup_{|x|\le Cr}v(x)/r^2)\).

## 2. Local interval error and shrinking averages

In the finite-order alternative, fix \(L\), now choose \(C>2(L+2)\), and put
\[
 Q_r(du)=r^{-m}g(ru)v(ru){\bf1}_{|u|\le C}\,du,\qquad
 \rho_r=\frac{\sqrt{\sup_{|x|\le Cr}v(x)}}r
        =O(r^{\beta/2-1}) .
\]
For all sufficiently small \(r\), \(\rho_r>0\) because of (5), and (4) bounds \(Q_r(\mathbb R)\) uniformly. The density \(r^{-m}S_{\mathrm{in},r}(rz)\) is obtained by averaging each source \(u\) of \(Q_r\) with the normals \(N(u,s/r^2)\), \(s\) uniform on \([0,v(ru)]\) where \(v(ru)>0\). All these Gaussians have standard deviation at most \(\rho_r\). Equation (2) with \(N=m\) gives
\[
  r^{-m}S_{\mathrm{in},r}(rz)=b_m z^m+O_L(r)
                      \quad (|z|\le L+2).                      \tag{7}
\]
On every source interval of length at most \(\rho_r\) contained in \([-L-1,L+1]\), enlarge the interval by \(\rho_r\) at each end. At least \(\mathbb P(|G|\le1)\) of each such source is retained; since (7) bounds the smoothed density, the interval's \(Q_r\)-mass is \(O_L(\rho_r)\). For any \(J\subset[-L,L]\), near each endpoint sum strips of width \(\rho_r\) and use Gaussian crossing bounds \(O(e^{-j^2/2})\) on the \(j\)-th strip. Sources at distance at least 1 from the endpoints are controlled by the total \(Q_r\)-mass and \(O(e^{-1/(2\rho_r^2)})\). The total difference between the sharp source window mass and the smoothed mass is \(O_L(\rho_r)\), uniformly in \(J\). Integrating (7) yields
\[
 \left|
 r^{-m-1}\int_{rJ}g(x)v(x)\,dx - b_m\int_Ju^m\,du
 \right|
                  \le C_L(\rho_r+r|J|).                       \tag{8}
\]
In particular, for the physical window \([r,r+\ell]\) with \(0<\ell=o(r)\),
\[
 \frac1{\ell}\int_r^{r+\ell} g(x)v(x)\,dx
 = b_m r^m+
 O\!\left(r^m\left[\frac{\ell}{r}+r+
       \frac{\sqrt{\sup_{|x|\le Cr}v(x)}}{\ell}\right]\right).
\]
Choosing \(\ell=r^{1+\theta}\), \(0<\theta<\beta/2-1\), gives a power-saving relative error \(O(r^\theta+r+r^{\beta/2-1-\theta})\). Both signs of the contact are covered by the same coefficient.

## 3. Signed reference primitive

Let \(w\ge0\) be a fixed smooth local reference with \(w(u)=O(u^4)\), extended by zero outside a fixed core, and put \(e=v-w=O(|u|^\alpha)\), \(\alpha>2\), \(q=\min(4,\alpha)\). The difference potential is
\(\int g(u)\int_{w(u)}^{v(u)}p_s(y-u)\,ds\,du\).
On the oriented window \(0\) to \(h\), split the input at \(|u|=2|h|\). Inside, replace the Gaussian window probability by the sharp window indicator. The error is at most
\[
 C|h|^\alpha\int_{|u|\le2|h|}
           e^{-\operatorname{dist}(u,\{0,h\})^2/(C|h|^q)}du
       =O(|h|^{\alpha+q/2}),
\]
because \(s\le C|h|^q\); use the absolute length \(|e(u)|\) even when the time interval is reversed. Outside, each fixed derivative has the analogue of (1) with \(\beta=q\); the signed time integral contributes a factor \(|e(u)|\). Taylor expansion and extension of the off-diagonal coefficients across the omitted core yield consistent polynomial jets.

The primitive of the reference potential has all-order Peano jets at zero. To show this, split at \(|u|=2|h|\) again. The off-diagonal part has the same Taylor bounds. If \(w\) is flat in value, the integrated inner contribution is bounded by \(\int_{|u|\le2|h|}g(u)w(u)du=O_N(|h|^N)\). Otherwise write \(w(u)=u^{2m_0}b(u)\), \(m_0\ge2\), \(b\) smooth positive, and choose the smooth signed square root \(\sigma(u)=u^{m_0}\sqrt{b(u)}\). The reference potential equals the output density of the mixture with source weight \(g w\), output map \(u+\sqrt{\tau}\sigma(u)G\) and \(\tau\) uniform on \([0,1]\). Restrict \(|G|\le|h|^{-\theta}\) for \(0<\theta<1\). On the inner core this map has derivative at least \(1/2\) and its image contains the oriented window. It fixes zero. If \(u_*(h,G,\tau)\) denotes the inverse of \(h\), the inner primitive is the normal-and-\(\tau\) average of \(\int_0^{u_*}g(u)w(u)du\). The inverse and its fixed-order derivatives on the image have bounds polynomial in \(|G|\), by implicit differentiation with the denominator bounded below by \(1/2\). Taylor expansion to any fixed order has an integrable Gaussian-weighted remainder of the same order in \(|h|\). Its coefficients at zero are the same for both signs of \(h\). The removed Gaussian tails are smaller than every power of \(|h|\), both for the integral and for the Taylor coefficients.

Hence there is a single bilateral polynomial \(P_*(h)\), with no terms below degree \(\alpha+1\), such that
\[
 \int_0^h g(u)e(u)\,du=P_*(h)+O(|h|^{\alpha+q/2}).
\]
Keep only polynomial terms of degree strictly below the requested error exponent. For \(\ell=h^{1+\delta}\), \(0<\delta<q/2-1\), it follows that the mean signed residual on \([h,h+\ell]\) is
\[
 P_*'(h)+O(h^{\alpha+\delta}+h^{\alpha+q/2-1-\delta})
\]
when weighted with \(g\). In particular any leading signed tangent at order \(\alpha\), under locally \(L^1\) or locally uniform convergence, is zero unless \(\alpha\) is an integer, and is a single monomial of degree \(\alpha\) in the integer case. No sign condition on the residual is needed.

## Local verification

The companion script ran under Python 3.13.18. SymPy 1.14.0 verifies the heat derivative identity and, for the quartic signed-root prototype with constant source weight, the Gaussian-averaged inverse-map primitive \(h^5/5+2h^7+55h^9\) through degree nine. mpmath 1.3.0 at 70 digits calibrates endpoint leakage. A unit uniform source with Gaussian standard deviation \(\rho\) loses \(2\rho/\sqrt{2\pi}+o(\rho)\) from its original unit window. For a continuous signed residual bump \(\pm h^5 z(1-z)\) at \(x=h+h^2 z\), \(0<z<1\), over the reference \(w=x^4\), the leakage from the outside into \([0,h]\) divided by \(h^7\) converges to \(\pm0.0208003547537603\ldots\). This illustrates the orders of the two comparison bounds. The checks do not assert that those source profiles give exact Gaussian output.
