# An entropy proof for the small-parameter center

Write \(g\) for the standard Gaussian density. Suppose \(f\ge0\), \(\mathbb E_\gamma f=1\), and \(\int\exp(t x-t^2f(x)/2)\gamma(dx)=1\) for every real \(t\). For exactly the center \(y=0\), put \(q=z^2/2\), \(A=g(0)\), \(\mu(dx)=xg(x)\mathbf 1_{\{x>0\}}dx/A\) and \(m=\int x\,d\mu=\sqrt{\pi/2}\). Set \(w(x)=e^{-qf(x)}\), \(Z=\int w\,d\mu\), and \(d\mu_q=w\,d\mu/Z\). The proof uses no upper bound on \(f\). The derivative passing through the defining integral is immediate for bounded \(f\), and also for \(q>0\) by the boundedness of \(f e^{-qf}\) locally in positive \(q\).

Throughout \(0<z\le1/100\), \(q\le1/20000\). Since \(1-w\le qf\), splitting at \(x=4\) and using the Gaussian mean of \(f\) gives
\[
 1-Z\le 4q/A+e^{-8}<1/1000,\qquad
 m-\int xw\,d\mu\le 16q/A+\left(4+\tfrac14\right)e^{-8}<1/250 .
\]
The tail terms use \(\int_4^\infty x e^{-x^2/2}dx=e^{-8}\) and one integration by parts for the second moment, followed by \(\int_4^\infty e^{-x^2/2}dx\le e^{-8}/4\). The rational estimates \(A^{-1}<8/3\) and \(m>5/4\) follow from \(3.14<\pi<22/7\), and \(e^{-8}<1/2500\) follows from \(e>8/3\) and \((8/3)^8>2500\). In particular, with \(C_q=\int x\,d\mu_q\) one has \(C_q>623/500\) and \(C_q\le m/Z\).

Let \(F_q=\int f\,d\mu_q\). Since the probability density of \(\mu_q\) relative to \(\gamma\) equals \(xw(x)/(AZ)\) on \(x>0\), its entropy is bounded by
\[
 D(\mu_q\|\gamma)
   =\int\log(x/(AZ))\,d\mu_q-qF_q
   \le\log\pi-2\log Z-qF_q .
\]
Here we used Jensen for \(\log x\), \(m/A=\pi\) and \(C_q\le m/Z\). The exact transform gives \(\lambda C_q-\lambda^2F_q/2\le D(\mu_q\|\gamma)\) by the entropy variational inequality. Use \(\lambda=9/5\), \(\log\pi<23/20\), and \(-\log Z\le(1-Z)/Z<1/999\); it follows that
\[
 (\tfrac{81}{50}-q)F_q>
        \tfrac95\,\tfrac{623}{500}-\tfrac{23}{20}-\tfrac2{999}>0 .
\]
For the logarithmic bound, the first seven terms in the series for \(e^{23/20}\) already exceed \(22/7\).

The right-wing center resolvent is
\(R(q,0)=2\int_0^\infty e^{-qf(x)} k_q(x)g(x)dx\) with \(k_q(x)=\sinh(\sqrt{2q}\,x)/\sqrt{2q}\). Since \(k_q(x)\ge x\) and its \(q\)-derivative is nonnegative, the preceding inequalities and the Round 122 free integral give
\[
 \frac{-\partial_qR(q,0)}{2A}
 \ge ZF_q-\tfrac23 e^q
 >\frac{999}{1000}\,
   \frac{(9/5)(623/500)-23/20-2/999}{81/50}
    -\frac{2}{3(1-1/20000)}
 >0.0059 .
\]
The free integral bound follows by differentiating
\(\int_0^\infty k_q(x)g(x)dx=A\int_0^1e^{q(1-s^2)}ds\), which is an identity by completing the square. The last sign is checked by exact rational arithmetic in **verify.py**. The reflected conclusion follows by replacing \(f(x)\) by \(f(-x)\). The endpoint \(q=0\) is included for bounded \(f\), with \(k_0(x)=x\). This center estimate is independent of any upper bound on the mark.
