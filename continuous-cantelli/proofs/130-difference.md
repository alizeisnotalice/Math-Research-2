# Round 121 MacBook CPU refinement of the first negative Bernstein test

Use the exact Gaussian-marginal, three-time-layer non-graph source of round 5. Write \(e=(3/5)e^{-16/3}\) and let \(p_v\) be the centered Gaussian density with variance \(v\). The densities of \(Q_t\) on \((0,1/2), (1/2,1), (1,3/2)\) are respectively
\[
  p_t(y),\qquad p_t(y)-e p_{t-3/8}(y-2),\qquad e p_{t-7/8}(y-2).
\]
The remainder vanishes beyond \(3/2\). Let \(W(t)=e^{-jt}(1-e^{-t})^k\), \(W_*=(j/(j+k))^j(k/(j+k))^k\) for positive \(j,k\), and \(d(y)=\int W(t)Q_t(y)dt\).

## Certified total order 949

For \((j,k)=(551,398)\), and every \(y\in[2.9360917,2.9360920]\),
\[
 -1.969032\cdot10^{-9}\le\frac{\sqrt{2\pi}\,d(y)}{W_*}\le-1.393998\cdot10^{-9}.
\]
Thus a nonzero nonnegative smooth test supported in the interior of this interval gives \(D_{551,398}(\psi)<0\). The previous round certified total order 1709. In independent 70-digit heat-time quadrature at the midpoint, \(d(y)/W_*=-6.6522293075301893\times10^{-10}\).

The certificate uses directed interval arithmetic on a composite midpoint rule. A Gaussian term is evaluated at each midpoint, with one-sided limits across the time-layer jumps. Its interval includes the uncertainty over the full spatial box. For one cell of length \(h\) the quadrature error is bounded by \(h^3\sup|\partial_t^2(Wg)|/24\), where \(g\) is the Gaussian factor with its common \(1/\sqrt{2\pi}\) removed. The second derivative is enclosed using \((Wg)''=Wg[(\log(Wg))'^2+(\log(Wg))'']\). Writing \(v=t-b\), \(d=y-m\) for a component, the two log derivatives of \(g\) are
\[
 -\frac1{2v}+\frac{d^2}{2v^2},\qquad \frac1{2v^2}-\frac{d^2}{v^3},
\]
and those of \(W\) are \(-j+k/(e^t-1)\) and \(-ke^t/(e^t-1)^2\). Bounds use extrema of \(W\) and \(g\) on each cell. On \([0,1/8]\), both positive factors are increasing, so the cell is bounded by length times the right endpoint integrand. Its scaled upper bound is less than \(7\times10^{-136}\). The rest of the mesh has \(h=2^{-12}\). The program uses 48-digit directed Decimal arithmetic, interval multiplication for integer powers, and one adjacent representable value on each side of correctly rounded `exp` and `sqrt`.

## Numerical check of the preceding diagonal

At every \(t\) the only possible negative density occurs for \(1/2\le t<75/128\). The source/free ratio in that range is
\[
 (3/5)\sqrt{\frac{t}{t-3/8}}\exp\left[-K(t)(y-16t/3)^2\right],\qquad K(t)=\frac{3/8}{2t(t-3/8)}.
\]
The peak is at most \(6/5\), and \(K(t)\ge1024/675\). Since \(\log(6/5)<1/5\), any negative point lies strictly within distance \(3/8\) of \(16t/3\). The full possible spatial support of negative density is therefore contained in \([55/24,7/2]\). In particular this argument certifies nonnegativity of every \(D_{j,k}\) spatial density outside that interval.

A local numerical cross-check enumerated every pair on \(j+k=948\), at 18,303 equally spaced spatial sites in \([55/24,7/2]\), giving 17,369,547 profiles. It integrated positive and negative Gaussian contributions separately and inspected the signed ratio \((P-N)/(P+N)\). All sampled ratios were positive. The smallest refined ratio, using order 256 Gauss–Legendre quadrature and one-dimensional minimization, was approximately \(2.55549104\times10^{-5}\) at \((551,397,y\simeq2.9315364)\). The identity \(D_{j,k}=D_{j+1,k}+D_{j,k+1}\) shows that a **rigorous** positive diagonal at order 948 would certify positivity at all smaller orders. The sampled diagonal is numerical evidence, and global minimality of the order-949 witness remains unproved.

## Four-way gate and arithmetic receipts

The MacBook CPU stage used the common 07:12:00 UTC gate. The intersection of all 15 worker intervals on the remote clock was 07:12:00.023797–07:13:31.016191 UTC, or 90.992394 s. The production batch generated 58,930,100,224 jittered spatial profiles in 8,392 batches and used 1,329.719 cumulative CPU seconds. The lowest numerical negative order found was 949, at \((551,398)\). All 15 workers were runnable at all nine samples. After the startup sample, summed worker CPU percent ranged from 1445.6 to 1475.9. The three `top` idle samples were 0%, 0% and 0.6%. The supervisor and workers all exited zero. A separate `ps` found none of their PIDs, and 42 fetched file hashes agree. Seven post-run SSH measurements bracket remote minus local clock offset at +0.008566 to +0.149368 seconds, conditional on steady offset during those measurements.

Remote software was Python 3.12.13, NumPy 2.5.4 and SciPy 1.18.1. Worker \(i\) uses NumPy generator seed \(1211201+i\), with one thread per process. The original remote command was `bash run.sh` using `refine.py` and `scan.py` in the remote project directory. To reproduce a quick run after the original gate, use `python refine.py 6 pilot.json --gate 0 --duration 0.15` with BLAS thread counts set to one. For the original 15-shard duration use `run.sh` in an unused output directory. `summarize.py` reads the fetched production receipts.

The lower-order interval result is reproduced by `python certify_midpoint.py 551 398 2.9360917 2.9360920 12 certificate-551-398.json`. The independent controls are `python check.py` and `python diagonal.py`. The local Python was 3.13.18 with mpmath 1.3.0. The sources, exact parameters, receipts and outputs are kept with this report.
