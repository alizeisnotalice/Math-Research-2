# Fractional variance contact detected by weak local tests

Let \(\phi=p_1\), \(X,Y\) be independent standard normal random variables, \(v\ge0\) be bounded and continuous, and put \(Z=X+\sqrt{v(X)}Y\). At a point \(a\), assume for \(\alpha>2\)
\[
 v(a\pm u)=u^\alpha(c_\pm+o(1)),\qquad u\downarrow0,
 \qquad c_\pm\ge0.
\]
Write \(V(z)=c_+(z_+)^\alpha+c_-((-z)_+)^\alpha\). The following expansion holds with no derivative assumptions on \(v\). For each compactly supported smooth real test \(\psi\) and each integer \(m>\alpha-2\),
\[
\begin{split}
 &\mathbb E\psi((Z-a)/h)-\int_{\mathbb R}\psi(y/h)\phi(a+y)dy\\
 &\qquad-\sum_{j=0}^{m-1}\frac{h^{j+1}}{j!}C_j\int_{\mathbb R} z^j\psi(z)dz
 =\frac{h^{\alpha-1}\phi(a)}{2}\int_{\mathbb R} V(z)\psi''(z)dz
       +o(h^{\alpha-1}),
\end{split}
\tag{1}
\]
where
\[
 C_j=\int_{v(a+u)>0}p_{v(a+u)}^{(j)}(-u)\phi(a+u)du.
\]
All the displayed background coefficients are finite. The density term on the right is the ordinary locally integrable function
\(V''=\alpha(\alpha-1)[c_+(z_+)^{\alpha-2}+c_-((-z)_+)^{\alpha-2}]\), interpreted under the integral. No boundary distribution occurs because \(V\) and its first derivative vanish at zero.

## Proof of the expansion

Choose \(K\ge1\) so that \(\psi\) is supported on \([-K,K]\), let \(u=X-a\), and split sources at \(|u|=2Kh\). In the near field put \(u=hz\) and \(s=\sqrt{v(a+hz)}/h\). Taylor expansion of the test and vanishing odd Gaussian moments give, uniformly in \(|z|\le 2K\),
\[
\mathbb E\psi(z+sY)-\psi(z)=\tfrac12s^2\psi''(z)+O_\psi(s^4).
\]
Only bounded derivatives of \(\psi\) are used. The variance is \(O(h^\alpha)\) uniformly there, so integrating the remainder costs \(O(h^{2\alpha-3})=o(h^{\alpha-1})\). After division by \(h^{\alpha-1}\), the main term converges to \(\phi(a)\int V\psi''/2\) by the variance asymptotics and uniformly bounded source density.

In the far field, a source with \(v=0\) makes no contribution. For the others the tested kernel is \(h\int\psi(z)p_v(hz-u)dz\). Taylor expand the kernel argument about \(-u\) to degree \(m-1\). For every derivative used in its remainder, \(|y-u|\ge |u|/2\) when \(|y|\le Kh\). For \(w\ne0\) and \(t>0\), the Gaussian derivative formula implies
\[
 |p_t^{(m)}(w)|\le C_m |w|^{-m-1}e^{-w^2/(4t)}.
\]
Indeed extract \(|w|^{-m-1}\) and bound the Hermite polynomial times \(|w/\sqrt t|^{m+1}e^{-w^2/(4t)}\). For \(u\) near zero, \(v(a+u)\le C|u|^\alpha\), so the last bound is at most \(C|u|^{-m-1}e^{-c/|u|^{\alpha-2}}\), integrable independently of \(h\). Away from zero the same derivative is uniformly bounded by the distance estimate, regardless of how small \(v\) is. Hence the integrated far remainder is \(O(h^{m+1})=o(h^{\alpha-1})\). The far Taylor coefficients can be extended from \(|u|>2Kh\) to all \(u\) with \(v>0\), because the added near-source derivatives at \(-u\) have exponentially small integral by the same Gaussian bound. This proves (1).

## Exclusion under a smooth output and a one-sided version

If \(Z\) has a smooth density near \(a\), fix a nonnegative smooth nonzero \(\chi\) supported in \((1,2)\) and take \(\psi=\chi^{(m)}\). All its moments of degree below \(m\) vanish, so the smooth output, input, and background together are \(O(h^{m+1})\) on this test. If \(c_+>0\) and \(\alpha\) is noninteger, integration by parts gives the contradictory nonzero coefficient
\[
 \int V(z)\psi''(z)dz=(-1)^{m+2}c_+\prod_{j=0}^{m+1}(\alpha-j)
       \int_1^2 z^{\alpha-m-2}\chi(z)dz\ne0.
\]
Use a reflected test if \(c_->0\). The result is local and applies in the exact Cantelli heat identity, where \(Z\) has a Gaussian output distribution. If \(\alpha\) is integral, the general necessary condition from (1) says that \(V''\) is a polynomial across zero, equivalent to \(c_-=(-1)^\alpha c_+\). Indeed tests with vanishing degree-\(<m\) moments annihilate precisely such polynomials, and on each open half-line the formula for \(V''\) is explicit.

For this one-sided test the other variance exponent need not match. If the tested side has \(v(a+u)\sim c u^\alpha\) for \(c>0\), noninteger \(\alpha>2\), while the opposite side satisfies \(v(a-u)=O(u^\beta)\) for any \(\beta>2\), retain \(\psi\) supported inside \((1,2)\). In the opposite near field \(u<0\), the baseline is zero and reaching the support has probability at most \(C e^{-c/h^{\beta-2}}\). In the far field use \(\min(\alpha,\beta)>2\) in the integrable remainder. Thus the same contradiction holds without a leading term on the other side. This variant also shows that the coefficient on the tested side cannot be canceled by a different, stronger polynomially flat contact on the opposite side.

## Why a pointwise expansion would need another condition

A local synthetic construction illustrates the weaker topology. Work at \(a=0\) with \(v_0(u)=|u|^\alpha\) near zero and a bounded extension. Let \(h_n=2^{-n}\), \(s_n=h_n^{\alpha/2}\), and \(\varepsilon_n=h_n^\gamma\) for \(0<\gamma<\alpha-2\). Choose nonnegative nonzero \(q\in C_c^\infty((-1/2,1/2))\), and set \(v_1(u)=v_0(u)[1+\varepsilon_n q((u-h_n)/s_n)]\) on these disjoint bump intervals for large \(n\), with \(v_1=v_0\) elsewhere. Both variances have the same asymptotic on both sides.

At \(y=h_n\), substitute \(u=h_n+s_n z\) in the contribution to the output density from the \(n\)th bump. The rescaled base variance is \((1+s_n z/h_n)^\alpha\to1\) uniformly on the bump. After dividing the density difference of the two models by \(\varepsilon_n\), a single differentiation of \(p_T(z)\) in \(T\), confined to a compact positive range of \(T\), gives the limit
\[
 \frac{\phi(0)}2\int (z^2-1)\phi(z)q(z)dz<0.
\]
Other bumps at \(u<0.6h_n\) contribute exponentially little at \(y=h_n\) by the separated Gaussian bound. The bumps at \(u>1.5h_n\) have a smooth-background expansion in \(y\) with remainder \(O(h_n^m)\), for arbitrary fixed \(m\), by the same derivative bound used in (1). Their background coefficients can include the near bumps at exponentially small cost. Hence after subtracting that polynomial background the pointwise density difference is of order \(\varepsilon_n\), which can dominate \(h_n^{\alpha-2}\). These nonexact models have the same leading weak coefficient in (1). They illustrate why no shared pointwise fractional correction follows from the given variance asymptotics alone.

## Executed local checks

`fractional_audit.py` verifies three vanishing moments and the nonzero integration-by-parts coefficient exactly with SymPy for \(\alpha=5/2,7/2\), using a finite-degree smooth-enough compact test. The script also runs 50- and 75-digit mpmath checks of the local weak coefficient for a Gaussian derivative test on a fixed scaled core, and of the synthetic narrow-bump correction. The hard cutoff in those numerical weak checks is in the source core only. Equation (1) covers arbitrary compact tests by its analytic near- and far-source estimates.
