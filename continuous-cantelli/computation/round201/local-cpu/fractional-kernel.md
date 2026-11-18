# A fractional derivative of the Gaussian heat kernel

Let \(p_t\) be the centered Gaussian density of variance \(t\), so that its Fourier transform is \(e^{-t\xi^2/2}\). For \(0<\kappa<1\), define both the Fourier-normalized and unnormalized Marchaud kernels at variance 1 by
\[
 K_\kappa(y)=\frac{1}{2\pi}\int_{\mathbb R}e^{iy\xi}|\xi|^{2\kappa}
                     e^{-\xi^2/2}\,d\xi,\qquad
 D_\kappa(y)=\int_0^\infty(p_1(y)-p_{1+t}(y))\,\frac{dt}{t^{1+\kappa}} .
\]
The Fourier multiplier factor between them is
\(\Gamma(1-\kappa)/(\kappa 2^\kappa)\). Indeed the scalar Marchaud integral \(\int_0^\infty (1-e^{-ut})t^{-1-\kappa}dt\) equals \(u^\kappa\Gamma(1-\kappa)/\kappa\) by integration by parts. The same argument applies via \(L^1\) domination near \(t=0\) and \(+\infty\), since \(\|p_1-p_{1+t}\|_1=O(t)\) for small \(t\).

Set \(z=y^2/2\) and \(c_\kappa=2^{\kappa-1/2}\Gamma(\kappa+1/2)/\pi\). Expand the cosine and integrate its Gaussian moments. Absolute convergence follows from \(|\xi|^{2\kappa}e^{-\xi^2/2+|y||\xi|}\in L^1\) at fixed \(y\). Kummer's transformation then gives
\[
 K_\kappa(y)
  = c_\kappa {}_1F_1(\kappa+\tfrac12;\tfrac12;-z)
  = c_\kappa e^{-z}{}_1F_1(-\kappa;\tfrac12;z).
\]
Thus the Marchaud representation has prefactor
\(\Gamma(1-\kappa)\Gamma(\kappa+1/2)/(\kappa\sqrt2\,\pi)\) multiplying \(e^{-z}{}_1F_1(-\kappa;\tfrac12;z)\).

## Unique zero and negative tail

Write \({}_1F_1(-\kappa;\tfrac12;z)=1-\kappa U_\kappa(z)\). Its coefficients after the constant are strictly negative, because \((-\kappa)_n=-\kappa(1-\kappa)_{n-1}\) for \(n\ge1\). Hence \(U_\kappa\) increases strictly from zero to infinity, giving exactly one positive zero \(r_\kappa\) in the physical \(y\)-coordinate, with \(z_\kappa=r_\kappa^2/2\). It is simple. The kernels are positive for \(|y|<r_\kappa\) and negative for \(|y|>r_\kappa\).

Differentiating the hypergeometric function and using the beta integral for \({}_1F_1(1-\kappa;\tfrac32;z)\) gives the positive integral
\[
 U_\kappa(z)=
 \frac{\sqrt{\pi}}{\Gamma(1-\kappa)\Gamma(\kappa+1/2)}
 \int_0^1 (e^{zt}-1)t^{-\kappa-1}(1-t)^{\kappa-1/2}\,dt.
\]
The singularity at zero is integrable because of \(e^{zt}-1\). This provides an exact integral equation for the root, \(\kappa U_\kappa(z_\kappa)=1\).

To obtain a controlled asymptotic, split the integral at \(t=1/2\). Below \(1/2\), use \(e^{zt}-1\le zt e^{zt}\) for \(z>0\). Above \(1/2\), set \(u=1-t\) and Taylor-expand \((1-u)^{-\kappa-1}\). Extend the resulting gamma moments to infinity at exponentially small cost. For fixed order \(N\), the resulting expansion is
\[
 U_\kappa(z)=\frac{\sqrt{\pi}}{\Gamma(1-\kappa)}e^z z^{-\kappa-1/2}
 \left[\sum_{j=0}^{N-1}
     \frac{(\kappa+1)_j(\kappa+1/2)_j}{j!\,z^j}+O(z^{-N})\right] .
\]
The remainder is uniform for \(0<\kappa\le\kappa_0<1\), by bounding the Taylor remainder on \(0\le u\le1/2\) and the endpoint integral at \(t=0\). This expansion retains the separate Gaussian term in \(1-\kappa U_\kappa\), which matters near the small-\(\kappa\) zero.

For fixed \(0<\kappa<1\), the negative leading tails are
\[
 K_\kappa(y)\sim
   -\frac{\Gamma(2\kappa+1)\sin(\pi\kappa)}{\pi}|y|^{-1-2\kappa},\qquad
 D_\kappa(y)\sim
   -\frac{2^\kappa\Gamma(\kappa+1/2)}{\sqrt\pi}|y|^{-1-2\kappa}.
\]
The first expression follows from the same beta integral and the reflection and duplication identities for Gamma. After the leading negative constant the first two relative coefficients in \(D_\kappa\) are
\[
 \frac{(2\kappa+1)(2\kappa+2)}{2y^2},\qquad
 \frac{(2\kappa+1)(2\kappa+2)(2\kappa+3)(2\kappa+4)}{8y^4}.
\]
For initial variance \(s>0\), the Marchaud kernel is \(s^{-\kappa-1/2}D_\kappa(y/\sqrt{s})\), so the root scales by \(\sqrt{s}\) and the leading negative tail constant stays the same. The function \(P_\kappa(y)={}_1F_1(-\kappa;\tfrac12;y^2/2)\) solves \((e^{-y^2/2}P_\kappa')'+2\kappa e^{-y^2/2}P_\kappa=0\) with initial data \(1,0\). Comparing the Wronskian of two solutions up to their first zero proves strict decrease of \(r_\kappa\) with \(\kappa\). The limit at \(\kappa=1\) is \(r_1=1\), since the polynomial at \(\kappa=1\) is \(1-y^2\).

## Small \(\kappa\) root

Set \(L=\log(1/\kappa)\) and \(w=\log(L/\pi)\). The zero diverges because the hypergeometric function tends to 1 locally uniformly as \(\kappa\to0\). The expansion of the positive integral first gives \(z_\kappa\sim L\). Taking logs and using the first three coefficients yields, at the root,
\[
 z_\kappa-\tfrac12\log z_\kappa
   = L-\tfrac12\log\pi-\tfrac1{2z_\kappa}-\tfrac5{8z_\kappa^2}
      +O(z_\kappa^{-3}).
\]
Terms with explicit factors \(\kappa\) are exponentially small in \(L\). Inverting the displayed equation, whose left-side derivative tends to 1, gives
\[
 z_\kappa=L+\frac{w}{2}+\frac{w/4-1/2}{L}
           +\frac{-w^2/16+3w/8-7/8}{L^2}
           +O\!\left(\frac{(1+\log L)^3}{L^3}\right).
\]
Returning to \(r_\kappa=\sqrt{2z_\kappa}\) gives \(r_\kappa^2=2L+\log L-\log\pi+o(1)\). A distinct useful form uses the negative Lambert branch:
\[
 Z_0=-\tfrac12W_{-1}(-2\pi\kappa^2),\qquad
 z_\kappa=Z_0-\frac{1}{2Z_0-1}+O(Z_0^{-2}).
\]
The \(1/L\) and \(1/L^2\) algebra is cross-checked in exact SymPy by a formal reciprocal-\(L\) expansion.

## Reproduction and scale

Run **python symbolic.py** and **python parallel_verify.py**. The latter uses 12 CPU processes and 51 decimal digits. It compares the Kummer expression to independent variance integration of the Marchaud kernel and to a direct Gaussian-cosine integral on \([0,20]\) for a subset. It checks the negative tail against two correction terms and finds small-\(\kappa\) roots in logarithmic units to avoid underflow. A separate hypergeometric bisection at 65 digits checks the roots up to \(L=160\). These computations calibrate the formulas proved above. They use high precision to test the analytic claims and do not provide interval enclosures.
