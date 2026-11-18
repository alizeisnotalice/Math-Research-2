# Matched tail-area coefficients near a quadratic origin

Assume \(f\) is bounded and nonnegative, with
\[
f(u)=c u^2+b u^3+d u^4+o(u^4)\quad (u\to0),\qquad c>0,
\]
and no other zero on the positive axis. Suppose on the positive tail
\[
f(x)=1+Kx^{-3}+Lx^{-4}+Px^{-5}+o(x^{-5}).
\]
Write \(A=e^{1/(2c)}/\sqrt c\), and let \(S\sim N(1/c,1/c)\). Define
\[
A_1=-\frac{bA}2 E[S^3],\qquad
A_2=A E[-S^2/2-dS^4/2+b^2S^6/8].
\]
The exact Gaussian transform at positive \(t\) forces
\[
K=2A,\qquad
L=K^2/4+2A_1,\qquad
P=2A_2-6K+KL/2-K^3/24.
\]

Indeed the compact contribution around zero, after substitution \(s=tu\), is
\(C_t=A/t+A_1/t^2+A_2/t^3+o(t^{-3})\). Shrink that compact interval so that \(f(u)\ge(c/2)u^2\). A Gaussian bound then controls all fixed polynomial moments of the difference quotients, or one may restrict to \(|s|\le t^\eta\) with \(\eta\) sufficiently small and bound the omitted tails exponentially. Expanding \(e^{s-cs^2/2}\) times the remaining exponential supplies the displayed coefficients. Away from the compact interval and the far tail, the negative real source has the bound \(e^{-rt}\) after discarding the nonnegative variance term. Any fixed intervening positive source interval has a strictly positive minimum of \(f\), so contributes only \(O(e^{-c' t^2})\).

On the tail substitute \(x=t+w\). Its missing mass relative to the standard shifted Gaussian is
\[
D_t=E\left[1-\exp\left(-\frac{t^2}{2}(f(t+W)-1)\right)\right],
\qquad W\sim N(0,1),
\]
where the expectation is restricted to the far-tail source. Truncate \(|w|\le t^\eta\). The rest is negligible by Gaussian tails and eventual \(f(x)-1=O(x^{-3})\), taking the fixed start of the tail sufficiently large. In detail, on sources \(x\le t/2\) in the far tail, make \(|f(x)-1|\le1/16\) and use \(|x-t|\ge t/2\). On \(x\ge t/2\), use the \(O(t^{-3})\) bound on the deficit. In the retained region,
\[
\begin{aligned}
\frac{t^2}{2}(f(t+w)-1)
={}&\frac{K}{2t}+\frac{L-3Kw}{2t^2}
 +\frac{P-4Lw+6Kw^2}{2t^3}+o(t^{-3})+
 O((1+|w|^3)t^{-4}),\\
D_t={}&\frac K{2t}+\frac{L/2-K^2/8}{t^2}
+\frac{P/2+3K-KL/4+K^3/48}{t^3}+o(t^{-3}).
\end{aligned}
\]
For exactness the compact contribution and missing tail mass agree up to negligible errors, giving the relations.

The weighted tail-area integrand is
\[
B(x)=1-(1+\delta)^{-1/2}\exp[-x^2\delta/(2(1+\delta))],\qquad\delta=f(x)-1.
\]
Its first and second coefficients coincide with those of \(D_x\). Its third coefficient is \(P/2+K/2-KL/4+K^3/48\), differing from that of \(D_x\) by \(-5K/2\). In the numerical grid take \(P=0\), \(K=2A\) and \(L=K^2/4+2A_1+\Delta\). Then
\[
t^2(D_t-C_t)\to\Delta/2,\qquad
T\left(\int_T^{2T}B(x)\,dx-\int_T^{2T} C_t\,dt\right)\to\Delta/4 .
\]
The next scaled residuals tend to \(3K-KL/4+K^3/48-A_2\) and \((3/8)(K/2-KL/4+K^3/48-A_2)\), respectively. The sampled profiles are reference cores and tails. Matching coefficients cannot establish an exact global extension.

The calculation integrates the compact area by interchanging the two real integrals and completing the square in the inner normalized tilt variable. Independent binary64 adaptive quadrature checks 47-digit calculations. In either method the Gaussian-scaled source and tail variables are truncated to \([-32,32]\). In the local profiles, \(f(u)\ge(21/40)u^2\), so the omitted compact terms admit a Gaussian bound, and the omitted far-tail terms admit the bounds above. The area integrand on \(x/T\in[1,2]\) is integrated directly.
