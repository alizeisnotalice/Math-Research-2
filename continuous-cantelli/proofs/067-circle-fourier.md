# Periodic variance against a uniform source on the circle

Let \(X\) be uniform on \(\mathbb R/(2\pi\mathbb Z)\), let \(Z\) be independent standard Gaussian, and put \(Y=X+\sqrt{v(X)}Z\pmod {2\pi}\). If \(v(x)=v_0+V(kx)\) with positive integer \(k\), its Fourier coefficient is
\[
 \mu_n=\frac1{2\pi}\int_{-\pi}^{\pi}e^{inx-n^2v(x)/2}\,dx .
\]
Partitioning into \(k\) periods gives \(\mu_n=0\) unless \(n=kj\), and
\[
 e^{(kj)^2v_0/2}\mu_{kj}
 =\frac1{2\pi}\int_{-\pi}^{\pi} e^{iju-(kj)^2V(u)/2}\,du .
\]
The selection rule is exact.

## Cusp and quadratic limits for a single sinusoidal well

For \(V(u)=c(2|\sin(u/2)|)^\alpha\) with \(c>0\), \(0<\alpha\le2\), put
\[
 L=(c(kj)^2/2)^{1/\alpha},\qquad \eta=j/L,\qquad
 T_\alpha(\eta)=\frac1\pi\int_0^\infty e^{-t^\alpha}\cos(\eta t)\,dt .
\]
Then
\[
 L e^{(kj)^2v_0/2}\mu_{kj}=T_\alpha(\eta)+O(L^{-2}).
\]
The error is uniform in \(\eta\) and in \(\alpha\) lying in any fixed compact subinterval of \((0,2]\). To see this, substitute \(t=L u\) on the positive half of the defining integral. For \(0\le u\le\pi\), \(b(u)=2\sin(u/2)/u\ge2/\pi\) and \(1-b(u)^\alpha=O(u^2)\). By the mean value theorem the difference of the two kernels is at most \(C L^{-2}t^{\alpha+2}e^{-\kappa t^\alpha}\). This has a uniform finite integral. The omitted reference tail \(t>\pi L\) is also \(O(L^{-2})\). The same argument after one Taylor correction gives
\[
 L e^{(kj)^2v_0/2}\mu_{kj}
 =T_\alpha(\eta)+\frac{\alpha}{24\pi L^{2}}\int_0^\infty
         t^{\alpha+2}e^{-t^\alpha}\cos(\eta t)\,dt+O(L^{-4}).
\]
For this remainder use \(b(u)^\alpha=1-\alpha u^2/24+O(u^4)\) and the exponential second-derivative bound between \(b(u)^\alpha\) and \(1\).

For a cusp \(\alpha<2\), at fixed \(k\) and \(j\to\infty\) we have \(\eta\to0\). Therefore the scaled value tends to \(\Gamma(1/\alpha)/(\pi\alpha)\), with errors \(O(L^{-2}+\eta^2)\). For a quadratic minimum \(\alpha=2\), the fixed-\(k\) limiting value is
\[
 \frac1{2\sqrt\pi}\exp\big[-1/(2ck^2)\big].
\]
Here the smooth phase survives because \(\eta=\sqrt{2/c}/k\). Increasing the period count drives \(\eta\to0\) uniformly over \(j\ge1\) for every fixed \(\alpha\le2\). For \(\alpha=2\) there is an exact independent formula
\[
 e^{(kj)^2v_0/2}\mu_{kj}= e^{-c(kj)^2} I_j(c(kj)^2),
\]
where \(I_j\) is the modified Bessel function.

All permitted harmonics in this family are strictly positive after removing the positive floor factor. Let \(p=\alpha/2\). If \(p=1\), expand \(e^{a\cos u}\) into products of two exponential power series. Its \(j\)-th coefficient is positive for \(a>0\). For \(0<p<1\), write \(x^p=p\Gamma(1-p)^{-1}\int_0^\infty(1-e^{-tx})t^{-1-p}dt\) by integration by parts. Truncation of this positive Lévy integral represents \(e^{-a x^p}\) as a compound Poisson Laplace transform. The truncated laws are tight, since their tail at \(R\) is bounded by \([1-e^{-aR^{-p}}]/(1-e^{-1})\). Any weak limit has transform \(e^{-a x^p}\) and no atom at zero. Applying this transform to \(x=2-2\cos u\) represents the circular heat amplitude as a positive mixture of kernels \(e^{-2t}e^{2t\cos u}\), each with positive coefficients for \(t>0\).

## Rigidity for finitely many isolated minima with one-sided powers at most two

More generally assume \(V\ge0\) has only finitely many zeros \(u_i\) and
\[
 V(u_i+\sigma s)=c_{i,\sigma} s^{\alpha_{i,\sigma}}(1+o(1))
 \quad(s\downarrow0),\qquad
 c_{i,\sigma}>0,\quad0<\alpha_{i,\sigma}\le2 .
\]
Take \(\alpha_*=\max_{i,\sigma}\alpha_{i,\sigma}\). Each side with smaller exponent is \(o(n^{-2/\alpha_*})\) after removing the floor factor. At a dominant side substitute \(s=t\,n^{-2/\alpha_*}\) with \(n=kj\). Dominated convergence follows from the corresponding positive lower power bound. At fixed \(k\) the normalized contribution is
\[
 e^{iju_i}\cdot\frac1{2\pi}\int_0^\infty e^{-c_{i,\sigma}t^{\alpha_*}/2}
       e^{i\sigma \theta t}\,dt +o(1),
 \qquad
 \theta=\begin{cases}0&\alpha_*<2,\\1/k&\alpha_*=2.\end{cases}
\]
The off-minimum compact region is exponentially small. The real part of each displayed half-line integral is strictly positive. In the quadratic case this follows by taking half the Fourier transform of its even Gaussian extension. By simultaneous approximation there are arbitrarily large integers \(j\) for which every \(e^{iju_i}\) is arbitrarily close to one. For example apply the pigeonhole principle to multiples of the vector of phases on a uniform torus partition and increase the partition size. If the resulting nonzero return times are bounded, a recurring return time gives an exact period and arbitrarily large multiples. The normalized Fourier coefficients therefore have a subsequence with strictly positive real limit. The circular output is not uniform for any such fixed field, including asymmetric one-sided minima.

## A logarithmic transition from cusp to smooth

For fixed \(k,c\), let \(\alpha_j\le2\) and \((2-\alpha_j)\log j\to\beta\in[0,\infty)\). Then
\[
 \eta_j\longrightarrow\frac{\sqrt{2/c}}k e^{-\beta/2},
 \qquad
 L_j e^{(kj)^2v_0/2}\mu_{kj}
 \longrightarrow\frac1{2\sqrt\pi}
       \exp\big[-e^{-\beta}/(2ck^2)\big].
\]
This is obtained by taking the logarithm of \(\eta_j=(2/c)^{1/\alpha_j}
k^{-2/\alpha_j}j^{1-2/\alpha_j}\) and using dominated convergence in \(T_{\alpha_j}\). It resolves why taking a fixed cusp exponent first and sending it to two need not reproduce the fixed-period quadratic limit.
