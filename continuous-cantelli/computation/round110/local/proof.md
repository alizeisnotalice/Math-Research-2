# A sparse profile for testing full Gaussian windows

Let \(G\) be standard normal and define \(I(\lambda)=\mathbb E\exp[-\lambda^2W(\lambda+G)/2]\). The unit-window integral \(I_1\) uses the same integrand restricted to \(|G|\le 1\). For \(\rho=4\), put \(R_j=4\rho^j\), \(j\ge0\), and \(w_j=e^{-aR_j^2}/8\). Let \(\psi(t)\) equal 1 for \(|t|\le 1/2\), \(2(1-|t|)\) for \(1/2\le |t|\le1\), and 0 elsewhere. The continuous profile is \(W(r)=\delta-(\delta+\gamma)\sum_j\psi((r-R_j)/w_j)\), where \(a=1/32,\delta=1/4,\gamma=1/8\). The windows are disjoint.

For fixed \(c>0\) along \(\lambda=cR_N\), put \(E(u)=au^2+(1-u)^2/2\). Then
\[
\lim_{N\to\infty}\lambda^{-2}\log I(\lambda)
=\max\left\{-\delta/2,\ \gamma/2-\inf_{\ell\in\mathbb Z}E(\rho^\ell/c)\right\}.
\]
For the lower bound, \(W\le\delta\) gives the baseline and one plateau at \(j=N+\ell\) gives the exponent \(\gamma/2-E(\rho^\ell/c)\). Its Gaussian distance can be bounded by \(|R_j-\lambda|+1/16\), with \(O(\lambda)\) effect on the logarithm for fixed \(\ell\). For the upper bound use the baseline plus the entire dip contributions. On dip \(j\), \(W\ge-\gamma\) and the Gaussian distance is at least \((|R_j-\lambda|-1/8)_+\). For \(R_j\le C\lambda\), this yields the exponent \(\gamma/2-E(R_j/\lambda)\) with log error \(O_C(\lambda)\) uniform in \(j\). There are \(O_C(\log\lambda)\) such indices. For \(R_j>C\lambda\) and \(C\ge2\), the factors \(e^{-aR_j^2}\) alone yield an upper bound \(e^{\gamma\lambda^2/2-aC^2\lambda^2}\) times a uniformly bounded multiplier. Indeed consecutive squared centers increase by a factor of 16, so the sum beyond the first such center is bounded by a geometric series. Taking \(C\) sufficiently large completes the upper bound.

The minimizer of \(E\) is \(16/17\). At \(c=5/4\) the winning ratio is \(4/5\), which gives the full rate \(9/400\). The distance from \(\lambda\) to the nearest dip tends to infinity, so \(I_1=\mathbb P(|G|\le1)e^{-\delta\lambda^2/2}\) eventually and has rate \(-1/8\). At \(c=2\) the winning ratio is \(1/2\), which gives the full rate \(-9/128\). Thus \(I-1\) has both signs arbitrarily far along the same continuous profile.

There are finite sign controls too. For \(\lambda=320\), the dip at \(R=256\) has plateau width \(e^{-2048}/8\). The heat factor is \(e^{6400}\), the Gaussian distance is at most \(64+1/16\), and the normalization exceeds \(1/3\). With \(\log 24<4\), the plateau alone gives \(\log I>2295\). Every dip is farther than one unit, so \(I_1<1\). At \(\lambda=512\), the baseline is \(e^{-32768}\). The first four dips have distance at least \(256-1/8\), total length at most 1, and heat factor at most \(e^{16384}\). For the remaining dips, use the Gaussian density bound by 1 and \(R_{4+k}^2/32 =32768\cdot16^k\ge32768+k\). Their sum is at most \(e^{-16384}/[4(1-e^{-1})]\). These bounds imply \(I(512)<6e^{-15000}<1\). The script certificate.py reproduces the arithmetic.

The theorem concerns a continuous test profile and one Gaussian functional. A separate argument would be needed to realize the profile as an exact graph in the Cantelli problem.
