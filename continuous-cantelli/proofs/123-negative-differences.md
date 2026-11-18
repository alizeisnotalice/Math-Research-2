# Real-space Bernstein tests for exact non-graph heat sources

For a finite positive space–time source \(\mu(ds,dx)\), write
\[
 Q_t=P_t\delta_0-\int_{s\le t}P_{t-s}\delta_x\,\mu(ds,dx),\qquad
 D_{j,k}(\psi)=\int_0^\infty e^{-jt}(1-e^{-t})^k Q_t(\psi)\,dt .
\]
All the tests below act on real-space Gaussian heat kernels. The section-six binomial formula and the integral agree by finite expansion and integration, including \(j=0\) for sources with bounded stopping horizon and exact terminal cancellation. Direct integration avoids cancellation between separately computed Laplace values at high order.

## Platforms and the Laplace formula

If the source has a time atom at \(t_0>0\), then
\(Q_{t_0}-Q_{t_0-}=-\mu(\{t_0\},dx)\) in the weak sense. The event \(s\le t\) and \(P_0=I\) include this term. Single times carry no weight in the time integral defining \(D_{j,k}\). Right continuity of \(Q_t(\psi)\) follows from the contraction bound for \(P_u\psi\), its continuity in \(u\ge0\), and dominated convergence on the finite source. Thus the passage from nonnegativity for almost every time to nonnegativity at every time remains valid in the presence of platforms. For the constant graph \(f=1\), \(Q_t=p_t\,dx\) for \(0<t<1\) and \(Q_t=0\) for \(t\ge1\), including exactly \(t=1\).

For \(q>0\), set \(a=\sqrt{2q}\). Resolvent integration gives, for \(y>0\),
\[
 R(q,y)=a^{-1}\left(e^{-ay}-\int e^{-qs-a|y-x|}\,\mu(ds,dx)\right)
 =\frac{2}{a}\int_{x>y}e^{-qs}\sinh(a(x-y))\,\mu(ds,dx).
\]
The last step substitutes \(\int e^{ax-qs}\,d\mu=1\) and cancels exactly for \(x\le y\). Both steps apply to the explicit finite-layer positive sources here. For a graph with marginal \(\gamma\), the integral uses \(s=f(x)\). Expansion at \(q=0\) gives the displayed section-six formula
\[
 B(y)=-R'(0,y)=2\int_{x>y}\left((x-y)f(x)-\tfrac13(x-y)^3\right)\,d\gamma(x).
\]
In particular \(B'(y)/2=-\int_{x>y}f\,d\gamma+\int_{x>y}(x-y)^2\,d\gamma\) and \(B''(y)/2=f(y)p_1(y)-2\int_{x>y}(x-y)\,d\gamma\) for continuous \(f\). The reflected pieces match in value and first derivative at zero by \(E[Xf]=0\) and \(Ef=1\), respectively. The constant-graph numerical control in `independent.json` evaluates this Laplace resolvent by Gaussian tail functions and independent heat-time quadrature at 70 decimals.

## Gaussian-marginal exact negative witness

In the round 5 example let \(e=(3/5)e^{-16/3}\), \(q_0=p_{1/8}(\cdot-2)\), \(r_0=p_{5/8}(\cdot-2)\). The source on times \(1/2,1,3/2\) is \(e q_0\), \(p_1-e(q_0+r_0)\), \(e r_0\). It is positive, sums in space to \(p_1\), and its exponential identity equals one for all complex parameters by cancellation of the shifted Gaussian rates. In the three time slabs the densities of \(Q_t\) are
\[
 p_t(y),\qquad p_t(y)-e p_{t-3/8}(y-2),\qquad e p_{t-7/8}(y-2),
\]
followed by zero at and beyond \(3/2\). These formulas use the heat semigroup at time-layer endpoints and treat each slab integral with its own one-sided endpoint contribution.

Take \((j,k)=(1000,709)\), \(I=[2.829,2.830]\) and \(W(t)=e^{-jt}(1-e^{-t})^k\). Write \(d(y)=\int W(t)Q_t(y)\,dt\). The certified bound, rounded here outward, is
\[
 -1.871\cdot10^{-6}\ \le\ \frac{\sqrt{2\pi}\,d(y)}{\max W}\ \le\ -1.161\cdot10^{-6},\qquad y\in I .
\]
Therefore every nonzero nonnegative \(\psi\in C_c^\infty(I^\circ)\) gives \(D_{1000,709}(\psi)<0\). This preserves exact positivity of the source, exact Gaussian space marginal and the full exponential identities. The model has three distinct conditional source times and lies in the non-graph relaxation.

The bounds come from rectangular integration on 2048 dyadic cells in each time slab, taken uniformly over the full spatial interval. For each Gaussian the extrema in variance lie at an endpoint or at variance \((y-m)^2\), and it decreases in \(y>m\). The lower Bernstein weight bound comes from time endpoints, and its possible interior maximum is \((j/(j+k))^j(k/(j+k))^k\). The program propagates 48-digit Decimal enclosures with directed arithmetic for rational operations, correctly rounded elementary `exp` and `sqrt` expanded by one adjacent representable value, and integer powers by interval multiplication. The integrated quantity drops the positive common factor \(1/\sqrt{2\pi}\). At the midpoint \(y=2.8295\), independent 70-digit quadrature gives raw density \(-1.3600968723960739749\times10^{-510}\).

## Two-layer control from section nine

For \(w=11/20\) the two time layers are \(w p_{1/8}\) at \(1/2\) and \(p_1-wp_{5/8}\) at 1. The second source is positive since \(w^2/(5/8)=121/250<1\). The rate matching \(1/8-1/2=5/8-1\) gives the full exponential identity. The spatial marginal has variance \(29/40\), fourth moment \(381/160\), so it is non-Gaussian. The signed density in the middle slab is \(p_t-wp_{t-3/8}\). Its maximum source/free ratio is \(w\sqrt{t/(t-3/8)}\) at zero, and nonnegativity begins exactly at \(t=50/93\), a delay of \(7/186\) from \(1/2\).

At \(t=1/2\) exact SymPy determinant signs for the leading principal minors of the even-moment matrix are positive for \(k=0,\ldots,12\) and negative for \(k=13\). The output includes explicit exact coefficients of a monic degree-26 polynomial whose squared expectation is negative. The interval routine additionally gives a real-space witness with \((j,k)=(2500,1744)\), \(I=[0.001,0.002]\), and a scaled density in \([-4.638,-3.702]\times10^{-4}\) throughout \(I\).

## Reproduction and compute receipts

Use `python certify.py 1000 709 2.829 2.830 11 certificate-1000-709.json` and `python certify.py 2500 1744 0.001 0.002 11 certificate-two-2500-1744.json --model two`. The certificate requires only the Python standard library. Use `python independent.py` with mpmath 1.3.0 and SymPy 1.14.0 for the independent quadrature, control and rational witness. The local Python was 3.13.18.

MacBook scan software was Python 3.12.13, NumPy 2.5.4 and SciPy 1.18.1, in 15 single-thread processes for each run. To recreate the initial scan, copy `scan_initial.py` as `scan.py` with `run.sh` to an unused remote directory and run `bash run.sh`. To recreate the deep scan, use `scan.py` and `bash run_deep_final.sh`. The initial grid tested \(0\le j<320\), \(0\le k<432\), with 159,390,720 three-layer and 106,306,560 two-layer profiles. The refined grid covered \(320\le j<1800\), \(\lfloor .57j\rfloor\le k\le\lfloor .84j\rfloor\), \(y\in[2.4,3.3]\), totaling 870,595,512 profiles. The earliest negative sample among the scanned refined grid by \(j\) was \((551,398,2.9361328125)\). Its higher-order value after scaling by \(\max W\) was \(-6.6518219181\times10^{-10}\); it is a numerical candidate. The 96- and 192-point quadrature comparisons are recorded in `summary.json` and the per-worker files. No random seed was used.

The initial and refined production wrappers ran 06:46:44–06:46:46 UTC and 06:51:07–06:51:09 UTC, with 3.155 and 17.174 cumulative worker CPU seconds. All worker receipts and supervisor receipts give exit status zero. A follow-up `ps` found none of their PIDs, and file hashes agree across transfer. An intermediate calibration used the array index as a displayed label for \(k\), corrected in the refined production run, and is excluded from the summary. Its results are unnecessary for reproduction.
