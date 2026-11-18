# Strict heat-energy uniqueness on a monotone two-wing boundary

Let \(r(t)\ge0\) be continuous and nondecreasing on a time interval. Let \(\eta,\lambda\) be finite positive measures supported on the curves \((t,\pm r(t))\), whose time marginals have no atoms. Write \(\sigma=\eta-\lambda\) and \(\mu=\eta+\lambda\). The causal heat potential of a measure \(\nu\) is
\[
U_\nu(t,x)=\int_{s<t}p_{t-s}(x-y)\,d\nu(s,y).
\]
Assume the two positive potentials agree on the exterior and their common outward trace has finite integral with respect to \(\mu\). Also assume that for \(|\sigma|\)-almost every interior time on a positive-radius wing,
\[
r(t)-r(s)=o(\sqrt{t-s})\qquad(s\uparrow t).
\tag{1}
\]
Only the backward limit is needed. The following proof establishes \(\sigma=0\), including for singular time marginals.

## Unregularized energy and measurable good sets

For a fixed target \((t,\pm r(t))\), every earlier source lies between \(-r(t)\) and \(r(t)\). Moving the target further outward increases the spatial distance to each source. The positive causal potentials therefore increase to their unregularized boundary values as the outward displacement decreases to zero, by monotone convergence. In particular their boundary values agree. With \(K((t,x),(s,y))=p_{|t-s|}(x-y)\) off the equal-time diagonal, this gives
\[
\iint K\,d\sigma\,d\sigma=2\int U_\sigma(t,x)d\sigma(t,x)=0.
\]
The equal-time pairs are null by absence of time atoms. Absolute convergence follows by the same causal decomposition for \(\mu\): \(\iint K\,d\mu\,d\mu=2\int U_\mu\,d\mu<\infty\), and \(|\sigma|\le\mu\).

Let \(E_j\) be the set of times in a growing sequence of compact interior intervals where \(r(t)\ge1/j\) and \([r(t)-r(s)]^2\le t-s\) for every \(s<t\) with \(t-s<1/j\). These sets increase, and testing rational \(s\) is sufficient for their measurability by continuity. Condition (1) shows that they exhaust the charged times up to a null set. If a nontrivial initial interval has radius zero, it can be discharged separately: for each time in that interval exterior equality holds for all \(x\ne0\); integrating \(U_\sigma(t,x)=0\) over \(x\) gives \(\sigma(\{s<t\})=0\), and the two spatial wings coincide at 0 there. Thus the restriction of \(\sigma\) to that interval vanishes.

Set \(\sigma_j=\sigma|_{\{t\in E_j\}}\), and \(K_\delta=p_{|t-s|+\delta}(x-y)\) for \(\delta>0\). On the same wing and for \(0<t-s<1/j\), write \(d=|x-y|\), \(\tau=t-s\). The later time condition gives \(d^2\le\tau\), so
\[
\frac{p_{\tau+\delta}(d)}{p_\tau(d)}=
\sqrt{\frac{\tau}{\tau+\delta}}\exp\!\left(\frac{d^2\delta}{2\tau(\tau+\delta)}\right)\le e^{1/2}.
\]
For same-wing pairs with \(|t-s|\ge1/j\), \(K_\delta\) is uniformly bounded by \(\sqrt{j/(2\pi)}\). For opposite wings the spatial distance is at least \(2/j\), and \(\sup_{v>0}p_v(z)\le C/|z|\) gives another uniform bound. Thus dominated convergence is valid against \(|\sigma_j|\times|\sigma_j|\) and gives \(E_\delta(\sigma_j)\to E(\sigma_j)\) as \(\delta\downarrow0\).

## Fourier positivity and strictness

For a finite real signed measure \(\nu\) on spacetime, set \(a(\xi)=\xi^2/2\). Fourier inversion in space and the identity \(\widehat{e^{-a|t|}}(\omega)=2a/(a^2+\omega^2)\) for \(a>0\) give
\[
E_\delta(\nu)=\frac1{(2\pi)^2}\iint e^{-\delta\xi^2/2}
\frac{2a(\xi)}{a(\xi)^2+\omega^2}|\widehat\nu(\omega,\xi)|^2d\omega d\xi .
\tag{2}
\]
For \(\delta>0\), all exchanges are justified since the \(\omega\)-integral of the spectral factor is \(2\pi\) for \(\xi\ne0\), the spatial cutoff is integrable, and \(\widehat\nu\) is bounded. The slice \(\xi=0\) is null. The right side increases as \(\delta\downarrow0\). The preceding good-set limit thus shows that the ordinary unregularized energy of \(\sigma_j\) equals its nonnegative uncut Fourier energy.

Since \(\sigma_j\to\sigma\) in total variation, their Fourier transforms converge uniformly. Absolute kernel integrability for the full measure also gives \(E(\sigma_j)\to E(\sigma)=0\). Fatou in (2) then forces \(\widehat\sigma=0\) almost everywhere for \(\xi\ne0\). A finite measure's Fourier transform is continuous, so it vanishes everywhere, and uniqueness implies \(\sigma=0\). This argument does not interchange the spatially regularized kernels on the full singular measure.

## A sufficient trace integrability argument at an initial origin

Suppose the common exterior potential is \(p_t(x)\), \(\eta\) is the first-exit law of Brownian motion from the expanding interval, and the spatial marginal of \(\lambda\) has bounded density. Positive exterior mass gives
\[
2\bar\Phi(r(t)/\sqrt t)\le\lambda([0,t])\le C r(t)
\]
if \(r(0)=0\). Write \(s=r(t)/\sqrt t\). The Mills lower bound gives \(p_t(r(t))\le C(1+s^2)\), while the raw Gaussian bound is \(O(t^{-1/2}e^{-s^2/2})\). Combining these bounds across \(s^2\le2\log(1/t)\) and its complement yields \(O(1+\log(1/t))\). For \(\lambda\) the spatial bound is even \(O(1+\log(1/r(t)))\): at small \(x=r(t)>0\), \(s\) is at least the positive Gaussian tail quantile \(s_x\) solving \(\bar\Phi(s_x)=Cx\); \(s\phi(s)\) decreases for \(s>1\), and \(s_x\phi(s_x)/x\le C'(1+s_x^2)\le C''(1+\log(1/x))\). This is integrable under the bounded spatial marginal.

For \(\eta\), the same Mills estimate forces \(s^2\ge\tfrac12\log(1/t)\) for all sufficiently small \(t\), since otherwise \(e^{-s^2/2}\le C(1+s^2)\sqrt t\) is impossible. Reflection on the interval \([u,2u]\) gives \(\mathbb P(\tau\in[u,2u])\le 4\bar\Phi(r(u)/\sqrt{2u})=O(u^{1/8})\) for \(u\downarrow0\), using monotonicity of the barrier. Summing dyadic intervals and using no atom at 0 bounds \(\eta([0,t])\) by \(O(t^{1/8})\), which integrates the logarithmic trace. If \(r(0)>0\), the initial boundary trace is uniformly bounded directly. Boundedness of the trace away from time 0 follows from \(p_t\le(2\pi t)^{-1/2}\). These establish the absolute energy condition under the stated Gaussian-boundary circumstances.

The backward condition (1) must hold for the actual signed measure's charged times, rather than merely for Lebesgue almost every time, since the time marginals can be singular.
