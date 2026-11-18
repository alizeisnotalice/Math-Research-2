# Leakage and common jets at an infinitely flat seam

Fix a center \(a\) and a sufficiently small core \((-r,r)\) in the source coordinate \(u=X-a\). On each half-core, let the variance deficit be either zero or
\[
v_\sigma(u)=c_\sigma e^{-d_\sigma/u^{p_\sigma}},\qquad
0<u<r,\quad c_\sigma,d_\sigma,p_\sigma>0,\quad \sigma=\pm1.
\]
Let \(g=p_1\). For \(y\ne0\) sufficiently close to zero, write \(\mathcal D(y)\) for the local output density at \(a+y\), summing the two heat integrals over the positive-variance half-cores and the pushforward density over any zero-variance half-core. The latter contributes \(g(a+y)\) if \(y\) is on that side. For \(j\ge0\), define the normalized leakage jets
\[
J_j=\frac1{g(a)}
 \sum_{\sigma:\ v_\sigma>0}\int_0^r g(a+\sigma u)\,
 p_{v_\sigma(u)}^{(j)}(-\sigma u)\,du .
\]
Then for every fixed integer \(N\ge0\), from both signs of \(y\),
\[
\frac{\mathcal D(y)-g(a+y)}{g(a)}
=\sum_{j=0}^N \frac{J_j}{j!}y^j+O(|y|^{N+1}).
\]
In particular \(J_0>0\) if either branch is positive.

## Fixed-order asymptotics

For the integration region \(u\ge2|y|\), Taylor expand each heat kernel in the output coordinate. For intermediate arguments the displacement has magnitude at least \(u/2\). Gaussian differentiation supplies a uniform integrable bound of the form \(C_N v_\sigma(u)^{-(N+2)/2}e^{-u^2/(16v_\sigma(u))}\). The omitted portions of the integrals defining the coefficients, over \(0<u<2|y|\), decay faster than every power of \(|y|\).

For \(u<2|y|\), the branch with sign opposite to \(y\) has displacement at least \(|y|\). On the matching positive-variance branch, put \(w=|y|\), \(\tau=\sqrt{v_\sigma(w)}\), and choose a neighborhood \(|u-w|\le\rho\) with \(\rho\) a sufficiently small fixed multiple of \(w^{p_\sigma+1}\). Outside that neighborhood the displacement is at least \(\rho\); the largest variance is at most \(v_\sigma(2w)\). For small \(w\), maximizing the kernel over these allowed variances still gives a bound smaller than every power of \(w\). Inside the neighborhood, use \(u=w+\tau t\) and the uniform bound \(\big|\log[v_\sigma(w+\tau t)/v_\sigma(w)]\big|\le C\). The contribution over \(|t|>w^{-1}\) is smaller than every power of \(w\) by a uniform Gaussian tail. For \(|t|\le w^{-1}\), the log variance ratio and the variation in \(g\) are both smaller than every power of \(w\), because \(\tau\) is. Expanding the scaled kernel and using the same Gaussian bound gives \(g(a+y)+O(w^N)\) for every \(N\). For a matching zero-variance branch the identity contribution is exactly \(g(a+y)\).

At \(y=0\), the local density defined by positive smoothing variance \(\varepsilon\) tends to \(g(a)(1+J_0)\). On \(u<\varepsilon^{1/4}\), the rescaling \(u=\sqrt\varepsilon\,t\) shows that each half-core contributes \(g(a)/2+o(1)\), since \(v_\sigma=o(\varepsilon)\) uniformly there. Outside that interval, \(v_\sigma(u)+\varepsilon\le C u^4\). The Gaussian kernel is then bounded uniformly near zero by \(C' u^{-2}e^{-c'/u^2}\), so dominated convergence gives the leakage on each positive-variance side.

If an exact solution has such a core, the output from \(|X-a|\ge r\) is analytic in the output coordinate near \(a\). Its positive-variance contribution has uniform complex Gaussian bounds because of the fixed displacement, and its zero-variance contribution has no mass in that interval. Thus analyticity of the target Gaussian forces the local density to have a convergent Taylor expansion there. Its leakage coefficients must obey \(\limsup_j |J_j/j!|^{1/j}<\infty\). Common finite-order jets across the two sides supply necessary coefficients. Following finitely many coefficients does not resolve this convergence condition.

## A rational center witness outside the unit interval

Take \(a=51/50,\ r=1/4\), the left plateau, and
\[
v_+(u)=\frac1{64}\exp[-3(r/u-1)].
\]
On \(15/64\le u\le1/4\), we have \(v_+\ge1/80\) since \(e^{-1/5}\ge4/5\), and \(v_+\le1/64\). Consequently \(u^2/(2v_+)\le5/2\) and \(g(a+u)/g(a)\ge571/800\). The rational bounds \(e^{5/2}<13\) and \(q_0>39/100\) give the strict lower bound
\[
J_0>\frac1{64}\frac{571}{800}\frac{78}{25}\frac1{13}.
\]
For \(M\ge1\), the normalized target density at \(a>1\) is maximized at \(M=a^2\), and
\[
\sup_{M\ge1}\frac{p_M(a)}{g(a)}
=\frac{e^{(a^2-1)/2}}a
\le \frac{250000}{249849}.
\]
The lower bound for \(1+J_0\) is strictly greater than this upper bound. Since the remaining contribution to the full smoothed density is nonnegative, this core is incompatible with every \(M\ge1\).
