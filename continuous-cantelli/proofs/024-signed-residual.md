# Signed power residuals and powers of a logarithm

Let \(p_t\) be the Gaussian density with variance \(t\), write \(\phi=p_1\), and suppose the exact output satisfies
\[
D_\varepsilon(y)=\int p_{v(x)+\varepsilon}(y-x)\phi(x)\,dx=p_{M+\varepsilon}(y),
\qquad v\ge0,\quad M>0,\quad\varepsilon>0.
\]
Fix \(x_0\) and a smooth nonnegative local reference \(w(u)=O(u^4)\) around \(u=x-x_0=0\). Choose a common symmetric core \(|u|<r<1\) and set \(w=0\) outside it when defining the reference output \(W_\varepsilon\). Let \(r_\pm(u)=v(x_0\pm u)-w(\pm u)\) for \(0<u<r\). For \(\alpha>2\), \(\beta\in\mathbb R\), assume
\[
r_\pm(u)=u^\alpha (\log(1/u))^\beta(c_\pm+o(1)).
\]
The coefficients \(c_\pm\) may have either sign. If either is nonzero, then \(\beta=0\), \(\alpha\) is an integer, and \(c_-=(-1)^\alpha c_+\). In particular the nonzero leading residual then has the coefficient of one smooth monomial on both sides.

## Gaussian derivative bounds and the reference output

Let \(He_j\) be the probabilists' Hermite polynomials and
\[
F_k(z,T)=T^{-(k+1)/2}He_k(z/\sqrt T)e^{-z^2/(2T)}.
\]
The Gaussian heat equation gives \(\partial_T^j F_k=2^{-j}F_{k+2j}\), and \(|F_m(z,T)|\le C_m\min(1,|z|^{-m-1})\) for \(T\ge1\). This follows from the uniform bounds on \(He_m(s)e^{-s^2/2}\) and on its product with \(|s|^{m+1}\), where \(s=z/\sqrt T\). The same second bound holds for all \(T>0\) and \(z\ne0\). For fixed \(\varepsilon>0\), every spatial derivative can pass through the integral by its uniform bound for variances at least \(\varepsilon\). With signed \(u=\delta z\), \(\delta=\sqrt\varepsilon\), the Hermite signs from spatial differentiation and from \(y-x=-u\) cancel at \(y=x_0\).

The derivatives \(W_\varepsilon^{(k)}(x_0)\) are uniformly bounded as \(\varepsilon\downarrow0\). In the core, expand \(F_k(z,1+w(\delta z)/\delta^2)\) in the second variable to degree \(J-1\) for an integer \(2J>k\). The remainder, before the prefactor \(\delta^{-k}\), is bounded by
\[
C\delta^{2J}|z|^{4J}\min(1,|z|^{-k-2J-1}).
\]
The integral to \(|z|<r/\delta\), after restoring \(\delta^{-k}\), is bounded by \(C'\delta^{2J-k}(1+(r/\delta)^{2J-k})=O(1)\). Taylor term \(j\) is a fixed scalar times
\[
\delta^{-k-2j}\int_{|z|<r/\delta} \phi(x_0+\delta z)w(\delta z)^j He_{k+2j}(z)e^{-z^2/2}\,dz.
\]
For \(k+2j>0\) expand the smooth weight to order \(k+2j-1\) in \(u=\delta z\). The polynomial cancels by full-line Hermite orthogonality after extending to the whole line, at only an exponentially small truncation cost. The Taylor remainder \(O(|\delta z|^{k+2j})\) gives \(O(1)\). The \(k=j=0\) term is directly bounded. Outside the core, the unscaled derivative kernel is bounded by \(C|u|^{-k-1}\) for all variances; the Gaussian source measure is finite. This proves the assertion.

## Reduction to a kernel on logarithmic scale

Choose \(\theta>0\) small enough that \(\theta<1\) and \(p:=\min(4,\alpha-\theta)>2+\theta\). The residuals are \(O(u^{\alpha-\theta})\), and \(v(x_0+u)+w(u)=O(|u|^p)\). Choose \(k\) as large as needed in each parity so that \(k+2>\alpha\) and \(k+4>\alpha-\theta+p\). As \(\delta\downarrow0\), the derivative comparison gives
\[
\delta^{k+2-\alpha}\big(D_{\delta^2}^{(k)}(x_0)-W_{\delta^2}^{(k)}(x_0)\big)
=\frac{\phi(x_0)}{2\sqrt{2\pi}}
\int_0^{r/\delta} z^\alpha He_{k+2}(z)e^{-z^2/2}
  g_k(\log(1/\delta)-\log z)\,dz+O(\delta^\eta)
\]
for some \(\eta>0\), where \(g_k(t)=e^{\alpha t}[r_+(e^{-t})+(-1)^kr_-(e^{-t})]\) for \(t>\log(1/r)\), and is zero elsewhere. To verify the error, interpolate the two variance parameters and subtract the first-order Taylor term at \(T=1\). The scaled pointwise second-order error is bounded by
\[
C\delta^{-2-\alpha}|r_{\operatorname{sgn}(z)}(\delta|z|)|[v(x_0+\delta z)+w(\delta z)]\min(1,|z|^{-k-5}).
\]
Its integral is \(O(\delta^{p-2-\theta})\) by the choice of \(k\). Replacing the source density by \(\phi(x_0)\) in the linear term gives \(O(\delta^{1-\theta})\) from the Gaussian tail and \(|\phi(x_0+\delta z)-\phi(x_0)|\le C\delta|z|\). The derivative integral outside the core is \(O(\delta^{k+2-\alpha})\). Exactness and the uniform reference bound make the left side \(O(\delta^{k+2-\alpha})\).

Set \(L=\log(1/\delta)\) and change variables \(s=\log z\). The scaled equation is therefore
\[
(K_k*g_k)(L)=O(e^{-\eta L}),\qquad
K_k(s)=e^{(\alpha+1)s}He_{k+2}(e^s)e^{-e^{2s}/2}.
\]
The factor \(\phi(x_0)/(2\sqrt{2\pi})\) is positive. The core cutoff is accounted for by the zero extension in \(g_k\).

## Mellin moments, parity, and the logarithm

Write \(d=k\bmod2\) and \(m_0(a)=\int_0^\infty z^a e^{-z^2/2}dz\). For \(a>2\), the Mellin moment is
\[
J_k(a):=\int_0^\infty z^a He_{k+2}(z)e^{-z^2/2}dz
=m_0(a+d)\prod_{j=d,d+2,\ldots,k}(a-j).
\]
Indeed, the Gaussian moment recurrence \(m_0(b+2)=(b+1)m_0(b)\) gives a monic polynomial of degree \((k-d)/2+1\) after dividing by \(m_0(a+d)\). Full-line Hermite orthogonality supplies precisely the displayed roots, since for each listed power the integrand extends evenly and its polynomial test has degree less than \(k+2\). The recurrence and polynomial formula are valid also at the listed nonnegative integer roots. They are simple. If \(a=d+2m\) and \(k=d+2n\), then
\[
J_k'(a)=m_0(a+d)\,2^n m!(-1)^{n-m}(n-m)!.
\]
By differentiation under the integral, \(\int K_k=J_k(\alpha)\) and \(\int sK_k(s)ds=J_k'(\alpha)\).

Both \(K_k\) and, when its integral is zero, its primitive \(H_k(s)=\int_{-\infty}^s K_k(t)dt\) have exponential decay at \(-\infty\) and superexponential decay at \(+\infty\). For either such function \(U\), the assumption \(g_k(t)=C_k t^\beta+o(t^\beta)\), \(C_k=c_++(-1)^kc_-\), implies
\[
(U*g_k)(L)=L^\beta\big(C_k\int U+o(1)\big).
\]
For \(s\le L/2\), use pointwise convergence and dominance by \(C(1+|s|)^{|\beta|}|U(s)|\) after division by \(L^\beta\). For \(s>L/2\), boundedness of the extended \(g_k\) on compact sets and its polynomial growth bound the ratio by a polynomial in \(L\) times the superexponentially small tail of \(U\).

If \(J_k(\alpha)\ne0\), the exponential bound on \(K_k*g_k\) forces \(C_k=0\). If \(J_k(\alpha)=0\), integration by parts yields \(\int H_k=-J_k'(\alpha)\ne0\). Fubini and the fundamental theorem of calculus yield
\[
(H_k*g_k)(2L)-(H_k*g_k)(L)=\int_L^{2L}(K_k*g_k)(t)dt=O(e^{-\eta L}).
\]
All integrals converge absolutely by the exponential tails and polynomial bound on \(g_k\). But the asymptotic of the left side is \(L^\beta[C_k(2^\beta-1)\int H_k+o(1)]\). Thus \(C_k=0\) whenever \(\beta\ne0\) also in the vanishing-moment parity. For nonintegral \(\alpha\), neither parity has a vanishing moment, so \(c_+=c_-=0\) for all \(\beta\). For integral \(\alpha\), the parity opposite to \(\alpha\) has a nonzero moment, enforcing \(c_-=(-1)^\alpha c_+\); for \(\beta\ne0\) the matching parity enforces the remaining equation, so both coefficients are zero.

A direct pointwise expansion of \(K_k*g_k\) to order \(L^{\beta-1}\) generally requires a stronger remainder. For example let \(g_k(t)=t^\beta+t^{\beta-1/2}\sin(\xi t)\) for large \(t\). Choose sufficiently small \(\xi\ne0\) so that \(\widehat K_k(\xi)\ne0\), which is possible because its derivative at zero is \(-iJ_k'(\alpha)\ne0\). The oscillating part of the convolution has leading term \(L^{\beta-1/2}\operatorname{Im}(e^{i\xi L}\widehat K_k(\xi))\), by the same dominated convergence after factoring out the exponential. The power-only part is \(O(L^{\beta-1})\) by splitting at \(|s|=L/2\) and using \(\int K_k=0\). The oscillating term is larger along subsequences. The geometric-difference argument above works with this size of remainder.

## Local verification

`mellin_audit.py` verifies the polynomial factorization and the derivative at every root for \(1\le k\le18\) by exact SymPy arithmetic. It separately constructs the primitive through a Gaussian times a polynomial in \(e^s\), checking the derivative algebraically. `mellin_results.json` records mpmath 60- and 85-digit calibration of the primitive integral and its geometric difference for \((\alpha,k,\beta)=(4,8,1/2)\) and \((3,7,-1/2)\), with a larger oscillating remainder than a pointwise first logarithmic correction. The finite numerical log window is \([-80,4]\), the input scales are 20, 100, 500, and 2000, and the proof relies on the analytic bounds.
