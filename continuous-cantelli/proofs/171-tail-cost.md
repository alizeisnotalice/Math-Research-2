# Early-interval heat cost in both Gaussian tails

Let \(a<b\) be fixed finite endpoints, \(g=p_1\), and
\[
 B(x)=\int_a^b\int_0^1p_s(x-u)\,ds\,g(u)\,du .
\]
The interval is used only as an early source. Write \(X=|x|\), \(\beta=b\) on the right wing and \(\beta=-a\) on the left wing. Define \(C_\beta=\sqrt{2/\pi}\,e^{-\beta^2}\). On either wing,
\[
 \frac{B(x)}{g(x)}=\frac{C_\beta e^{\beta X}}{X^3}
 \left[1+\frac{4\beta}{X}+\frac{11\beta^2-7}{X^2}
       +\frac{26\beta^3-48\beta}{X^3}+O(X^{-4})\right].
\]
In particular, the leading exponential on the left wing is \(e^{a x}\), and on the right wing it is \(e^{b x}\).

## Endpoint integration

Reflect the left wing to the right. Write the resulting interval as \([\alpha,\beta]\). Put \(w=X(\beta-u)\), \(r=X-u=X-\beta+w/X\), and \(\rho=r^2(1/s-1)/2\). Then exactly
\[
 \frac{B(x)}{g(x)}=\frac{C_\beta e^{\beta X}}{X^3}
 \int_0^{X(\beta-\alpha)}\!\!\int_0^\infty e^{-w-\rho}
 \frac{e^{2\beta w/X-w^2/X^2}}{(1-\beta/X+w/X^2)^2}
 \left(1+\frac{2\rho}{r^2}\right)^{-3/2}d\rho\,dw .
\]
The heat endpoint contributes the factor \(2/r^2\) before extracting \(X^{-3}\); the outer spatial endpoint contributes \(1/X\). For sufficiently large \(X\), the amplitude before multiplying \(e^{-w-\rho}\) is bounded by \(4e^{w/2}\), since \(r\ge X/2\). This proves the leading term by dominated convergence. For the displayed expansion one can Taylor expand on \(w,\rho\le X^{1/4}\), where derivatives in \(1/X\) through order four are bounded by a polynomial in \(w,\rho\) times an integrable exponential after weighting. The rest is exponentially small by the preceding bound. The coefficients then follow by integrating monomials against the two unit exponential factors. `symbolic.py` independently applies the factorial moments in exact SymPy arithmetic and checks the displayed coefficients.

## Inverse time increment

Let \(d=d(x)>0\) solve
\[
 \int_1^{1+d}p_t(x)\,dt=B(x).
\]
Its solution is unique since the integrand is positive and the integral diverges as \(d\to\infty\). Put \(A=B(x)/g(x)\), \(\Lambda=X^2/2\), and \(R=\Lambda d/(1+d)\). The exact inverse equation is
\[
 \Lambda A=\int_0^R e^v(1-v/\Lambda)^{-3/2}dv .
\]
If \(\beta<0\), then \(\Lambda A\to0\) exponentially, and \(d=A[1+O(X^2A)]\). If \(\beta=0\), then \(\Lambda A=O(X^{-1})\), and Taylor expansion of the integrand in time gives \(A=d+(X^2-1)d^2/4+O(X^4d^3)\). This yields
\[
 d=\sqrt{2/\pi}\,X^{-3}-\frac{1}{2\pi}X^{-4}+O(X^{-5}) .
\]
If \(\beta>0\), the inverse equation implies \(R\le\log(1+\Lambda A)=O(X)\). Thus \((1-v/\Lambda)^{-3/2}=1+O(X^{-1})\) throughout the integration range. Consequently
\[
 R=\beta X-\log X-\beta^2-\tfrac12\log(2\pi)+O(X^{-1}),
 \qquad
 d=\frac{2\beta}{X}+\frac{-2\log X+2\beta^2-\log(2\pi)}{X^2}
      +O\left(\frac{\log X}{X^3}\right).
\]
This distinguishes the two boundary cases \(b=0\) and \(b>0\) on the right, and \(a=0\) and \(a<0\) on the left. For an interval with \(a<0<b\) both inverse times decay only as \(X^{-1}\), with respective coefficients \(2b\) and \(-2a\). The calculation gives the isolated heat cost of a chosen early source.

## High-precision stress test

A MacBook batch used 100 decimal digits for both wings of \([-1,0]\), \([0,1]\), and \([1/2,1]\) at \(X=128,4096,1048576\). The scaled outer integral was computed with two independent evaluations of the inner heat factor. One uses the finite alternating Taylor expansion of \(\int_0^\infty e^{-\rho}(1+2\rho/r^2)^{-3/2}d\rho\). Its term after degree 44 bounds the error by \(91!!/r^{90}\), at most \(1.461\times10^{-119}\) for tested \(r\ge127\). The other uses the exact identity
\[
 \int_0^1p_s(r)ds=2g(r)-2r\overline\Phi(r)\qquad(r>0),
\]
evaluated through the complementary error function. For the scaled outer integral, any omitted part above \(w=256\) is bounded by \((4/0.98)e^{-0.98\cdot256}<4.52\times10^{-109}\) in the tested domain. The largest observed discrepancy between the two quadratures was \(3.04\times10^{-84}\).

The inverse equation was solved after factoring out \(e^R\) in the growing-exponential cases and after dividing the equation by \(\Lambda A\) in the small-cost cases. A second inversion check integrated in the original time coordinate after the endpoint blowup \(t=1+d(1-w/R)\) for growing cost. The largest absolute discrepancy in log cost was \(7.49\times10^{-96}\). These discrepancies are numerical cross-checks.

Representative results at \(X=2^{20}\) follow. The normalization of \(B/g\) is division by \(C_\beta e^{\beta X}X^{-3}\).

| \(\beta\) | Normalized \(B/g\) | Delay \(d\) | Relative error of the displayed delay approximation |
|---:|---:|---:|---:|
| 0 | 0.9999999999936335371 | \(6.920544077178315\times10^{-19}\) | \(-6.32\times10^{-12}\) |
| 1 | 1.000003814700903585 | \(1.907323563614084\times10^{-6}\) | \(-5.29\times10^{-11}\) |
| \(-1/2\) | 0.9999980926475018530 | \(2.219411237824526\times10^{-227714}\) | \(1.80\times10^{-17}\) |

For the last table row the approximation to \(d\) uses the displayed heat expansion through \(X^{-2}\). For the other two rows it uses the corresponding two-term inverse. Reproduction uses `bash run.sh` with `tail_cost.py` and `symbolic.py`. The run used Python 3.12.13, mpmath 1.3.0 and SymPy 1.14.0, and performed 18 numerical cases in six MacBook CPU processes. Their computation began 10:52:05.686344 UTC and ended 10:52:14.377226 UTC, totaling 38.058 cumulative CPU seconds. All six numerical processes and the symbolic check exited zero. An independent process query found no owned workers or supervisor, and all 26 transferred hashes matched. `summary.json` keeps full outputs and process receipts.
