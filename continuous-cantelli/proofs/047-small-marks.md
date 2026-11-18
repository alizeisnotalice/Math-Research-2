# Small marks at the critical positive-tail scale

Assume the setup of `040-short-windows.md`: exact entire identity, positive quadratic zero at the origin, no other positive zero, and \(f(x)\ge1\) for \(x\ge R\). Let \(M=\sup f>1\), \(A=e^{1/(2c)}/\sqrt c\), and \(D_t\sim A/t\). Write \(h(x)=f(x)-1\) on the tail, let \(\phi=p_1\), and use the marked deficit \(\omega_t\) with marginal
\[
d\mu_t(u)=\phi(u)(1-e^{-t^2h(t+u)/2}){\bf1}_{t+u\ge R}du.
\]
Assume \(\eta_t=\sup_{x\ge t/2}t^2 h(x)\to0\), which follows from the stated uniform \(x^2h(x)\to0\). Mark displacement under the Gaussian interpretation is \(w=(-V+i\sqrt V G)/t\) for \(G\sim N(0,1)\).

## An improved, finite-band critical law

For any \(0<\lambda<\sqrt{M/(M-1)}\), let \(\psi\) be the positive band-limited kernel in round 40, with Fourier support \([-4,4]\), and put \(\psi_b(y)=b^{-1}\psi(y/b)\), \(b=4/(\lambda t)\). Uniformly in \(a\in\mathbb R\),
\[
(\psi_b*\mu_t)(a)=D_t(\psi_b*\phi)(a)+O_\lambda(t^{-2}+\eta_t/t).\tag{A}
\]
The R40 marked transform is \(Q_t(s)=(t^2/(t+is)^2)\,[K_t(s)-J_t(s)]\), with
\[
K_t(s)=e^{-s^2/2}\int_{x<R}e^{is(x-tf(x))+s^2 f(x)/2}\,d\nu_t(x),
\qquad |J_t(s)|\le\gamma_t(x<R).
\]
For \(|s|\le t/2\), R40 gives the requisite Gaussian-weighted error with \(L^1(ds)\) size \(O(t^{-2})\). For \(t/2\le |s|\le\lambda t\), the real exponent in \(K_t\), apart from \(tx\), is at most
\[
-s^2/2+(s^2-t^2)f(x)/2\le -\tfrac12\min(s^2,Mt^2-(M-1)s^2)\le-c_\lambda t^2.
\]
The additional \(tx\) is at most \(tR\). The same exterior frequency band of the proposed Gaussian target and of \(J_t\) is Gaussian-small. Thus the inverse Fourier integral of the marked transform against this filter equals \(D_t(\psi_b*\phi)(a)+O_\lambda(t^{-2})\).

For the mark comparison on \(x\ge t/2\), \(V\le\eta_t\), and the first-order Taylor term in \(\psi_b(a-u-w)-\psi_b(a-u)\), averaged in \(G\), is relatively \(O_\lambda(V)\) by \(\mathbb E_G w=-V/t\). Taylor's integral remainder is also relatively \(O_\lambda(V)\): for \(j=0,1,2\) the derivative \(\psi^{(j)}\) satisfies the same exponential-imaginary and polynomial-real bounds as in R40, and
\[
\mathbb E_G\left[(|w|/b)^2(1+|w|/b)^4 e^{4|\Im w|/b}\right]\le C_\lambda V
\]
for \(0\le V\le1\). This follows by inserting \(|w|/b\le C_\lambda \sqrt V(1+|G|)\) and using a Gaussian exponential moment. Hence the contribution of the local mark replacement is at most \(C_\lambda\eta_t(\psi_b*\mu_t)(a)\).

For the remote sources \(R\le x<t/2\), the marked kernel is bounded in absolute value by \(C_\lambda t e^{\lambda\sqrt V|G|}\). The \(G\)-average is at most \(2C_\lambda t e^{\lambda^2V/2}\), so its product with the mark density yields the factor \(e^{(\lambda^2-1)V/2}\), up to polynomial factors in \(t\). On a fixed bounded source interval \([R,R_0]\), combine \(V\le(M-1)t^2\) with \(\phi(x-t)\le C e^{-t^2/2+R_0t}\); the result is exponentially small since \((\lambda^2-1)_+(M-1)<1\). On \([R_0,t/2]\), choose \(R_0\) so large that \((\lambda^2-1)_+h(x)/2<1/16\); the Gaussian source contributes at most \(C e^{-t^2/8}\). The unmarked remote comparison is even smaller. Since \(\eta_t\to0\), absorb the relative local error to first obtain \((\psi_b*\mu_t)(a)=O(t^{-1})\), then conclude (A).

The microscopic measures \(\Lambda_t(dy)=t^2\mu_t(d(y/t))\) satisfy
\[
\int g(a-y)\,d\Lambda_t(y)\longrightarrow A\phi(0)\int g(y)dy
\]
for fixed \(a\) and every Schwartz filter \(g\) whose Fourier transform is smoothly supported in \((-\lambda,\lambda)\). One can repeat the preceding calculation using \(t g(ty)\): its entire derivatives have the same exponential type and polynomially decaying real bounds after including a polynomial in the imaginary argument, by repeated integration by parts in the Fourier formula. Use the positive reference kernel in (A) to dominate its real and shifted Taylor bounds. Additional polynomial Gaussian moments on remote marks do not alter their strict exponential gaps. The positive reference kernel and (A) also bound \(\Lambda_t\) uniformly on every unit interval. Every vague subsequential limit is therefore a translation-bounded positive measure agreeing with Lebesgue density \(A\phi(0)\) on this open Fourier band.

## Why this does not give an unsmoothed subcritical law

An entire filter with spatial scale \(b=o(1/t)\) would need frequency support wider than any fixed multiple of \(t\). Under the small-mark condition, the nearby-source Taylor estimate can remain small whenever \(t b\gg\sqrt{\eta_t}\), but the remote/compact estimate above has an independent finite frequency range. This is a genuine issue for the current comparison. For example, a continuous synthetic profile can have a compact plateau \(f=M\) on \([a,b_0]\), linear loss from \(M\) at the adjacent two ends, a positive quadratic zero around zero, and tail \(h=x^{-3}\). At \(s=\lambda t\), \(\lambda>1\), the plateau's contribution to the compact transform is
\[
 e^{[(\lambda^2(M-1)-M)t^2]/2-i\lambda M t^2}
 \int_a^{b_0}e^{(1+i\lambda)tx}\phi(x)dx
 \sim\frac{\phi(b_0)}{(1+i\lambda)t}
 e^{[(\lambda^2(M-1)-M)t^2]/2+(1+i\lambda)t b_0-i\lambda M t^2}.
\]
The adjacent strips are smaller by \(O(t^{-1})\) relative to this leading size because the positive real mass near a linear edge has width \(O(t^{-2})\). Away from those strips keep a strict gap from \(M\) before cutting to the small tail. For \(\lambda>\sqrt{M/(M-1)}\) the compact term then grows exponentially in modulus. The profile here does not satisfy the exact identity; it demonstrates that the broad frequency estimate needs more than the local origin and size assumptions.

Finite-band data cannot identify an unsmoothed microscopic limit, as the two positive measures with constant density and density proportional to \(1+\epsilon\cos(\kappa y)\), \(0<\epsilon<1\) and \(\kappa>\lambda\), agree on all the above band tests. This latter pattern is consistent with small excess on a synthetic tail: for
\[
f(x)=1+2A x^{-3}(1+\epsilon\cos(\kappa x^2/2)),
\]
\(x^2h(x)\to0\), while along \(t_n=\sqrt{4\pi n/\kappa}\), \(t_n B(t_n+y/t_n)\to A(1+\epsilon\cos(\kappa y))\). Along phase-\(\pi\) sequences the sign switches. Windows of width \(t_n^{-2}\) have normalized averages tending to \(A(1\pm\epsilon)\). This synthetic profile is evidence of a resolution limit of the finite-band inference, without making any claim of an exact solution. The marked exact identity at higher frequencies could impose further cancellations, and this argument does not resolve that possibility.
