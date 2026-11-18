# Uniform center derivative from a radial entropy envelope

Let \(f\ge0\) have Gaussian mean one, and suppose \(\int e^{t x-t^2 f(x)/2}\,d\gamma(x)=1\) for all real \(t\). Put \(g_0=1/\sqrt{2\pi}\), \(z=\sqrt{2q}>0\), and let \(\rho\) be the pushforward of \(\gamma\) under \(x\mapsto(r,s)=(|x|,f(x))\). The \(r\)-marginal is half-normal and \(\int s\,d\rho=1\). The transform identities at opposite parameters give
\[
 \int e^{-q s}\cosh(zr)\,d\rho=1,\qquad
 \int_{x>0} e^{-q f(x)}\sinh(zx)\,d\gamma(x)
 =\int_{x<0} e^{-q f(x)}\sinh(z|x|)\,d\gamma(x).
\]
Define \(C=\sqrt{\pi/2}=1/(2g_0)\), \(w=e^{-qs}\), \(K=C\sinh(zr)/z\), and \(Z=\int wK\,d\rho\). Since the right-wing center resolvent is \(R(q,0)=2\int_{x>0} e^{-q f(x)}\sinh(zx)/z\,d\gamma(x)\), the second identity gives \(Z=R(q,0)/(2g_0)\). For positive \(q\), differentiate under the integral using the Gaussian tails and the bound on \(s e^{-qs}\) on compact positive \(q\)-intervals.

With \(\lambda=9/5\), \(\varphi=\log\cosh(\lambda r)-\lambda^2s/2\) satisfies \(\int e^\varphi d\rho=1\). Apply the entropy variational inequality to probability \(wK\,d\rho/Z\) relative to \(\rho\). All terms are integrable since \(K\) grows at most exponentially in \(r\), \(K\log K\) extends continuously by zero at the origin, and \(w s\) is bounded. It gives
\[
 a\int wK s\,d\rho\ge \int wK[\log\cosh(\lambda r)-\log K]d\rho+Z\log Z,
 \qquad a=\lambda^2/2-q.
\]
For the required range \(z\in[1/100,3/5]\) one has \(a>0\) and \(Z_0=1-q>0\). Use the supporting tangent \(Z\log Z\ge(1+\log Z_0)Z-Z_0\), and subtract \(a\int wK_q\,d\rho\), where \(K_q=\partial_q K\). With
\[
 H=K[\log\cosh(\lambda r)-\log K+1+\log(1-q)]-aK_q,
\]
we obtain \(-a\partial_qR(q,0)/(2g_0)\ge\int wH\,d\rho-(1-q)\).

For any real \(H\) and \(s\ge0\), the minimum of \(e^{-qs}H+(3/2)qs\) over \(s\) is
\[
 M(H)=\begin{cases}H,&H\le3/2,\\(3/2)[1+\log(H/(3/2))],&H>3/2.\end{cases}
\]
It is increasing in \(H\). Since the mean mark is one, integrating this pointwise inequality gives
\[
 -a\partial_q R(q,0)/(2g_0)\ge \int_0^\infty M(H(r,q))\,C^{-1}e^{-r^2/2}dr-1-q/2.
\]
The right side is certified uniformly larger than **0.0702584**, so \(-\partial_q R(q,0)/(2g_0)>0.04336\) on the full specified parameter interval. Equality of the opposite center wing follows from the same transform identity.

## Interval certification

The bound \(0\le K_q/K\le r^2/3\) follows termwise from the positive power series for \(\sinh(zr)/z\), since differentiating the \(n+1\) term and comparing with the \(n\) term gives the factor \(r^2/(2n+3)\). These series also prove that \(K\) and \(K_q\) are increasing in both \(r\) and \(z\). For \(0<r\le0.02\), \(K\le Cr e^{0.6r}<1\). In that region \(1+\log(1-q)\ge1-q/(1-q)>3/4\), while \(aK_q/K\le2r^2/3\), giving \(H>0\). This region contributes nonnegatively.

For \(r\ge12\), \(K\ge1\) and \(\log K\le1+1.6r\le2r\). Thus \(H\ge-K(2r+2r^2/3)\ge-r^2K\). Also \(M(H)\ge-r^2K\), since \(M=H\) for negative \(H\) and \(M\ge0\) for positive \(H\). Multiplying the half-normal density and using \(0.6r\le0.05r^2\) gives a negative-tail magnitude bounded above by
\[
 \int_{12}^\infty r^3e^{-0.45r^2}dr=65.8\,e^{-64.8}/0.405<1.171\times10^{-26}.
\]

The remaining range \(r\in[0.02,12]\) uses 1198 exact decimal spatial cells of width 0.01. There are 236 parameter cells of width 0.0025 in \(z\). On each rectangle, the monotonicity of \(K,K_q\), \(\log\cosh\), and the Gaussian density gives endpoint bounds. The program encloses \(H\) by interval arithmetic and applies the increasing function \(M\) to its lower endpoint. It multiplies that lower endpoint by the lower density endpoint when positive and by the upper density endpoint when negative. A rational alternating series for the Machin identity verifies \(C\in[1.2533,1.2534]\). Decimal basic operations use directed rounding at 48 digits, and `exp` and `ln` use correctly rounded values expanded on both sides by adjacent representable values. The weakest certified cell is \(z\in[0.5975,0.6]\), with exact lower endpoint printed in `summary.json`.

## Reproduction

MacBook Python 3.12.13 ran 15 CPU workers with `bash run.sh`, which calls `cover.py` and `interval.py`. The workers' common computation interval was 08:43:38.080112–08:43:46.929038 UTC, with 67.463 cumulative worker CPU seconds. All workers and their supervisor exited zero, confirmed absent by a separate PID query. All 41 transferred hashes were checked. `python verify.py` confirms exact coverage, positive margins and successful exit receipts. Grid construction and scalar evaluation are deterministic.
