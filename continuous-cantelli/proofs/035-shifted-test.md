# Critical-origin and tail-patch overlap

Write \(\varphi,\Phi\) for the standard normal density and CDF. For \(F(s)=\int e^{sx-s^2f(x)/2}\gamma(dx)\) and bounded \(J=[s_-,s_+]\), positivity and Tonelli give
\[
 \int_J F(s)\,ds
 =
 \int_{f>0}\sqrt{\frac{2\pi}{f(x)}}e^{x^2/(2f(x))}
 \left[\Phi\left(\sqrt{f(x)}s_+-\frac{x}{\sqrt{f(x)}}\right)
       -\Phi\left(\sqrt{f(x)}s_--\frac{x}{\sqrt{f(x)}}\right)\right]\gamma(dx)
 +\int_{f=0}\int_J e^{sx}\,ds\,\gamma(dx).
\]
An exact graph therefore gives precisely \(|J|\) total mass, so any sum of contributions from disjoint source regions is bounded by \(|J|\). The displayed formula includes arbitrary zero sets via the second integral.

For a fixed relative source interval \(I_t=t+[a,b]\), a fixed observation interval \(J_t=t+[c,d]\), and a tail cap \(f\le1+\kappa_t/t^2\) on \(I_t\), the lower bound for the tail average follows from
\[
 e^{sx-s^2f(x)/2}\gamma(dx)=e^{-s^2(f(x)-1)/2}\varphi(x-s)\,dx.
\]
Let \(G=(d-c)^{-1}\int_c^d[\Phi(b-u)-\Phi(a-u)]du>0\). Suppose the positive-origin core leaves remaining mass \(p_t+o(p_t)\) uniformly at \(s=t+u\) for \(u\in[c,d]\), where \(p_t\downarrow0\). Exactness forces \(\kappa_t>0\) eventually. Once positive, it forces
\[
 \kappa_t\ge
 \frac{2\log(G/(p_t+o(p_t)))}
 {\max_{u\in[c,d]}(1+u/t)^2}.
\]
This estimate uses no global lower bound on the tail. If \(\log(1/p_t)=o(t)\), it yields a lower bound \(2\log(1/p_t)+2\log G+o(1)\).

For the saturated model with \(L=\log t\) and \(p_t\sim (1/4-j)\log L/L\), \(j<1/4\), the leading bound becomes
\[
 2\log L-2\log\log L+2\log\big(G/(1/4-j)\big)+o(1).
\]
The uniformity needed here follows from the origin expansion because the observation shift is in a fixed compact interval and \(\log(t+u)-\log t=O(1/t)\).

Under the additional hypotheses of the conditional Gaussian conclusion from round 29, let \(p_t\) be the tail-window mass itself. Weak convergence gives
\[
 \frac{1}{p_t}\int_{[a,b]}e^{-t^2(f(t+z)-1)/2}\varphi(z)\,dz
 \ \longrightarrow\ \int_{[a,b]}\varphi(z)\,dz .
\]
A weighted mean-value argument and continuity give an actual point \(z_t\in[a,b]\) with
\[
 t^2(f(t+z_t)-1)=2\log(1/p_t)+o(1)
\]
in any fixed nondegenerate bounded interval. This additional observation uses the conditional convergence, unlike the preceding upper-cap exclusion.
