# Independent localization audit for the dyadic area law

Assume the exact Cantelli identity \(F(t)=1\) for a bounded nonnegative continuous \(f\), and assume \(f\ge1\) for \(x\ge R>0\). Let \(\nu_t(dx)=e^{tx-t^2 f(x)/2}\gamma(dx)\), \(c_t=\nu_t(x<R)\), and \(\phi=p_1\). The previous center-area proof is valid, including its allowance of negative zeros away from the origin. The following bounded-kernel comparison gives a more direct quantitative boundary estimate.

## Nonnegative kernel comparison

For \(x\ge R\), write \(h(x)=f(x)-1\) and set
\[
K(x,t)=\phi(t-x)(1-e^{-t^2h(x)/2}),\qquad
D_t=\int_R^\infty K(x,t)dx,
\]
so \(0\le K\le\phi(t-x)\). Completing the Gaussian square in \(t\) gives the exact pointwise expression
\[
B(x)=1-f(x)^{-1/2}e^{-x^2(1-1/f(x))/2}=\int_{\mathbb R}K(x,t)dt.
\]
Separating the exact identity at the cutoff gives \(D_t=c_t-\mathbb P(N(t,1)<R)\). For any \(b\ge a\ge\max(1,R)\),
\[
\left|\int_a^b B(x)dx-\int_a^b D_tdt\right|\le(\sqrt{2\pi}+2)(D_a+D_b).\tag{A}
\]
To prove this, subtract the two integrals using Tonelli and define the signed single-boundary error by integrating \(K(x,t)({\bf1}_{t\ge y}-{\bf1}_{x\ge y})\). This boundary expression is absolutely integrable because the crossing part of \(\phi(t-x)\) is integrable. On a crossing let \(d=|x-y|\), \(s=|t-y|\). Then \(|t-x|=d+s\) and \(\phi(t-x)\le \phi(x-y)e^{-s^2/2}\). The function \(u\mapsto1-e^{-u}\) is concave and nondecreasing on \([0,\infty)\), so
\[
1-e^{-t^2h(x)/2}\le\max(1,t^2/y^2)(1-e^{-y^2 h(x)/2})
\le(1+s)^2(1-e^{-y^2h(x)/2})
\]
for \(y\ge1\). Integration in \(s\) bounds the crossing error by \(D_y\int_0^\infty(1+s)^2e^{-s^2/2}ds=(\sqrt{2\pi}+2)D_y\). The two endpoint errors give (A).

A separate inequality is convenient when endpoint decay is unproved. For \(0<H<(b-a)/2\), the bound \(K\le\phi(t-x)\) and positivity give
\[
\int_{a+H}^{b-H}D_tdt\le\int_a^b B(x)dx+2\int_H^\infty\Phi(-s)ds.\tag{B}
\]
For the parts with \(x<a\) and \(x>b\), respectively, integrate \(\phi(t-x)\) in \(x\) first and then use the distance of the inner \(t\)-interval from the endpoint.

## The exact criterion at a quadratic origin

Suppose \(f(x)=c x^2+o(x^2)\) at zero, \(c>0\). For sufficiently small fixed \(r\in(0,R)\), the origin contribution \(\nu_t([-r,r])\) is \(A/t+o(t^{-1})\) as \(t\to+\infty\), where \(A=e^{1/(2c)}/\sqrt c\). The substitution \(z=tx\) and the lower bound \(f(x)\ge c x^2/2\) give a Gaussian majorant. The negative remainder has \(\nu_t(x\le-r)\le e^{-rt}\), regardless of its zeros. Define \(E=[r,R)\), \(g(t)=\nu_t(E)\). Up to immaterial endpoints and Gaussian tail from the baseline,
\[
D_t=A/t+g(t)+o(t^{-1}).\tag{C}
\]
The dyadic law \(\int_T^{2T}B(x)dx\to A\log 2\) holds if and only if \(\int_T^{2T}g(t)dt\to0\).

For sufficiency, put \(W_t=x-tf(x)\). Differentiating twice under the Gaussian integral gives \(\mathbb E_{\nu_t}W_t^2=\mathbb E_{\nu_t}f\le M\). Also \(|g'(t)|\le \sqrt{M g(t)}\) by Cauchy–Schwarz. Thus \(\sqrt{g(t)}\) is uniformly \(\sqrt{M}/2\)-Lipschitz, with the zero-measure case identically zero. This implies that \(\int_T^{2T}g\to0\) forces \(g(T),g(2T)\to0\): a positive lower bound at an endpoint forces a fixed interior mass in a fixed-width interval. Combining (A) with (C) proves the dyadic law. With no positive zero outside the quadratic origin core, \(g(t)\) is exponentially small by compact positivity, and the boundary error in (A) is \(O(T^{-1})\).

For necessity, assume the dyadic law and apply (B) with \(a=T,b=2T,H=\sqrt T\). Integrating the origin term in (C) over the trimmed interval gives \(A\log 2+o(1)\). The dyadic law then forces \(\int_{T+\sqrt T}^{2T-\sqrt T}g\to0\). The two shifted dyadic windows with left endpoints \(0.9T\) and \(1.2T\), after trimming each by the square root of its endpoint, cover \([T,2T]\) for large \(T\). Nonnegativity of \(g\) completes the argument.

## Center measure and sufficient local density conditions

The area formula from the original note says \(\int_\mathbb R g(t)dt=\rho(E)\), where
\[
 d\rho(x)=f(x)^{-1/2}e^{x^2(1/f(x)-1)/2}dx,\qquad m=x/f(x)
\]
where \(f>0\). Consequently \(\rho(E)<\infty\) is a sufficient local condition even in the presence of positive zeros. More generally, let \(\mu=m_*(\rho|_E)\). Then \(\int_T^{2T}g(t)dt\to0\) is equivalent to \(\mu([T,2T])\to0\). To see one direction, use the area formula on \(t\in[T,2T]\). Centers \(m\in[T/2,4T]\) have mass \(o(1)\) by the stated center condition. Centers below \(T/2\) have Gaussian crossing probability \(O(e^{-r T/4})\) and mass \(O(T)\) by the packing bound, using \(\sigma^2=m/x\le m/r\). Centers above \(4T\) give a summable geometric sequence of such exponential errors. Conversely any center in \([T,2T]\) lands in \([T/2,4T]\) with probability bounded below uniformly for large \(T\), again using \(\sigma^2\le2T/r\). The \(g\)-integral over the enlarged interval is \(o(1)\) by summing three dyadic intervals. In particular an eventual density \(d\mu/ds=o(s^{-1})\) suffices.

Positive zeros can make a nonnegligible conditional contribution. For the local model of round 6 at \(a_0>0\),
\[
f(a_0\pm e^{-\ell})=\frac{a_0^2}{2(\ell-\beta\log\ell+d)},
\]
its tilted local mass is asymptotic to \(2\phi(a_0)e^d\sqrt\pi\,(a_0t/2)^{1/2-\beta}\). Thus \(\beta=3/2\) adds \([4\phi(a_0)e^d\sqrt\pi/a_0]\log 2\) to the dyadic limit if it is an isolated additional contact and the other residuals are negligible. When \(\beta>3/2\), the local \(\rho\) integral is finite, since after \(u=e^{-\ell}\) its density in \(\ell\) is comparable to \(\ell^{1/2-\beta}\). These are local consequences, and the use of these models does not establish a continuous exact solution.
