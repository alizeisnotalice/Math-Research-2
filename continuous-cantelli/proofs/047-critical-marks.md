# Microscopic marked deficit at the \(1/t\) scale

Assume the setting of the positive-tail result in round 40. Thus \(F_f=1\), \(f\ge1\) eventually on the positive tail, \(f>0\) on the rest of the positive half-line away from zero, and \(f(x)=cx^2+o(x^2)\) on both sides of zero for \(c>0\). Write \(A=e^{1/(2c)}/\sqrt c\) and \(B(x)=1-f(x)^{-1/2}e^{-x^2(1-1/f(x))/2}\) on the positive tail.

## 1. Marked equation on the critical line

Use the positive marked measure \(\omega_t\) of round 40. Set
\[
 z=t(x-t),\qquad H_t(z)=t^2(f(t+z/t)-1),\qquad
 d_t(z)=1-e^{-H_t(z)/2}.
\]
All expressions below involving \(H_t,d_t\) are restricted to \(t+z/t\ge R\), where \(f\ge1\), and the measure is zero elsewhere. Let \(\phi\) be the standard Gaussian density and define the lifted measure
\[
 d\Lambda_t(z,H)=t\frac{\phi(z/t)}{\phi(0)}
             d_t(z)\,dz\,\delta_{H_t(z)}(dH).
\]
The depth \(H\) takes values in the compactification \([0,\infty]\). For finite positive \(H\), let \(K_H\) be exponential of rate \(1/2\) conditioned to be at most \(H\). Let \(K_0=\delta_0\), \(K_\infty=\mathrm{Exp}(1/2)\). The marginal in \(z\) of \(\Lambda_t\) is locally within a vanishing relative error of \(tB(t+z/t)dz\), by the pointwise comparison in round 40 and \(\phi(z/t)/\phi(0)\to1\). The round-40 uniform strip bound gives \(\Lambda_t([z,z+1]\times[0,\infty])\le C\) at every center \(z\), for large \(t\). Thus every sequence has a vaguely convergent subsequence, with a translation-bounded nonnegative limit \(\Lambda\).

Use the convention \(\widehat\psi(k)=\int e^{-ikz}\psi(z)dz\) and set
\[
 {\cal T}_V\psi(z)=\frac1{2\pi}\int \widehat\psi(k)
                      e^{ik(z-V)+k^2V/2}dk .
\]
Every limit satisfies
\[
 \int_{\mathbb R\times[0,\infty]}
       \left(\int {\cal T}_V\psi(z)K_H(dV)\right)d\Lambda(z,H)
                 = A\int_{\mathbb R}\psi(z)dz                 \tag{1}
\]
for Schwartz \(\psi\) with Fourier support compact in \((-1,1)\).

To prove (1), take \(\kappa<1\) with Fourier support in \([-\kappa,\kappa]\). Integration by parts in \(k\) bounds the test by
\[
 |{\cal T}_V\psi(z)|\le C_N(1+V)^{2N}e^{\kappa^2 V/2}(1+|z|)^{-N}.
\]
The \(K_H\)-average is continuous in \((z,H)\), including \(H=0,\infty\), and bounded by \(C_N'(1+|z|)^{-N}\) uniformly in \(H\). Indeed, the increasing exponential moment bound is integrable under \(K_\infty\), so it bounds all the truncated averages, and the limits at zero and infinity follow directly. The uniform strip bound controls the spatial truncation error in vague convergence. For the prelimit, the left side of (1) equals
\[
 \frac{t^2}{2\pi\phi(0)}\int\widehat\psi(k)Q_t(tk)\,dk ,
\]
with \(Q_t\) from round 40. Its estimate
\[
 |Q_t(s)-D_t e^{-s^2/2}|
 \le C_\kappa (|s|t^{-2}+s^2t^{-3}+e^{-c_\kappa t})e^{-s^2/2}
                  +C_\kappa\gamma_t(x<R)
\]
works uniformly on \(|s|\le\kappa t\): replace \(3/4\) by \(1-\kappa^2>0\) in the compact bounds for the quadratic origin. After \(s=tk\), the error tends to zero. The main term tends to \(A\widehat\psi(0)\), by \(tD_t\to A\) and \(2\pi\phi(0)=\sqrt{2\pi}\).

The explicit Fourier multiplier at depth \(H\) is
\[
 M_H(k)
 =\int e^{-ikV+k^2V/2}K_H(dV)
 =\frac{1-e^{-(1+ik)^2H/2}}{(1+ik)^2(1-e^{-H/2})} ,
 \quad M_0(k)=1,\quad M_\infty(k)=(1+ik)^{-2}.          \tag{2}
\]
The endpoint \(H=\infty\) in (2) is asserted in the open band only.

## 2. Models showing what the marked equation permits

If \(\Lambda=m(dz)\delta_0(dH)\), then (1) says that the signed distribution \(m-A\,dz\) has Fourier transform vanishing on \((-1,1)\). For example, \(m=A(1+\epsilon\cos(k_0z))\,dz\) with \(0<\epsilon<1\) and \(k_0\ge1\) obeys (1) and has nonconstant critical-window masses. The lattice \(m=2\pi A\sum_{n\in\mathbb Z}\delta_{2\pi n+z_0}\) obeys it as well, by periodization of the test. In fact any positive unmarked limit has atom mass at most \(2\pi A\): test by \(\operatorname{sinc}^2((1-\delta)(z-z_0)/2)\), whose integral is \(2\pi/(1-\delta)\). To fit the Schwartz band condition, first multiply this test by a nonnegative Schwartz band-limited cutoff tending to 1 by dilation, with bandwidth tending to zero. The \(O((1+|z|)^{-2})\) bound and strip mass bound justify the limit, and then \(\delta\downarrow0\). The lattice attains the bound.

Marks can also allow oscillations at low frequencies in the ordinary marginal. For any \(0<k<1\) and sufficiently small \(\epsilon>0\), take
\[
 \Lambda= \left(\frac A2-\epsilon\operatorname{Re}[M_\infty(k)e^{-ikz}]\right)dz\,\delta_0
   +\left(\frac A2+\epsilon\cos(kz)\right)dz\,\delta_\infty .
\]
The nonzero frequency contributions cancel in (1), but they do not cancel in the marginal. More generally \(M_\infty\) is the transform of the negative of a shape-2 rate-1 Gamma variable. If \(a_\infty\ge0\) is bounded by \(A\) and \(a_0(z)=A-\int_0^\infty a_\infty(z+r)r e^{-r}dr\), then \(a_0dz\,\delta_0+a_\infty dz\,\delta_\infty\) obeys (1). For decreasing nonconstant \(a_\infty\), the marginal can exceed \(A\,dz\) and still obey (1). These are models of the limiting equation. They do not assert the existence of exact graphs.

## 3. Supplementary conditions for critical-window asymptotics

A shallow-mass condition is
\[
  t\int_{-L}^{L} d_t(z)^2\,dz\longrightarrow0
                \quad\text{for every fixed }L>0 .              \tag{3}
\]
It is equivalent to weighted lift mass with \(H>\eta\) vanishing on these windows for every fixed \(\eta>0\). To see the equivalence, for fixed \(\eta\) bound \(d_t\) away from zero on \(\{H>\eta\}\), and split the integral of \(d_t^2\) according to \(H\le\eta\) and \(H>\eta\), using the uniform marginal mass bound. Thus (3) allows deep wells on sufficiently thin sets and ensures all limiting lifts are unmarked. A simpler stronger condition is \(f(x)-1=o(x^{-2})\) on the tail. Under (3), the only remaining obstruction in (1) is the Fourier gap ambiguity above.

- If every subsequential unmarked marginal has Fourier support contained in a closed subinterval of \((-1,1)\), the gap makes it \(A\,dz\).
- The same follows if every such marginal differs from \(A\,dz\) by a signed measure having one sign outside a compact set. For this last assertion assume the sign outside the compact is positive; negate otherwise. Apply the Fourier gap to a dilated positive band-limited Schwartz test tending pointwise to 1. The negative part of the signed measure is compact, so Fatou shows its positive part is finite. Its ordinary Fourier transform is then continuous and zero on the band. For any integer \(n\), test by \((1-\cos(hz))^n\), a trigonometric polynomial with frequencies in the band when \(n|h|<1\). Divide the resulting equality for the positive and negative parts by \(h^{2n}\), and use Fatou for the positive part and dominated convergence for the compact negative part. This bounds the \(2n\)-th moment of the positive part by the corresponding compact moment. Letting \(n\to\infty\) shows both parts are compactly supported. The entire Fourier transform then vanishes identically.

The shallow assumption can be replaced by spatial independence of the depth in the limiting lift. If \(\Lambda=m(dz)\otimes\pi(dH)\) for a probability law \(\pi\), (1) multiplies the Fourier transform of \(m\) by \(\int M_H(k)d\pi(H)\) on the open band. This multiplier is smooth and has no zeros there: for \(|k|<1\) and \(H>0\), the real part of \((1+ik)^2M_H(k)\) is positive because \(|e^{-(1+ik)^2H/2}|<1\), and the same real-part inequality holds at \(H=0,\infty\). Its value at zero is 1. Division therefore gives the same Fourier gap for the ordinary marginal as in the shallow case, including for depth laws with mass on deep marks.

Under either condition, vague uniqueness gives \(t^2\int_{t+a/t}^{t+b/t}B(x)dx\to A(b-a)\) for all fixed \(a<b\), and uniformly for \(a,b\) in bounded sets by uniform convergence of the distribution functions to the continuous limit. If the scaled profiles \(z\mapsto (t+z/t)^3(f(t+z/t)-1)\) also have local equicontinuity, convergence of these interval integrals yields the pointwise limit at the center. Indeed equicontinuity and the strip bound first bound the profile values on each compact set, using \(tB\ge c\min((t+z/t)^3(f-1),t)\); hence \(tB\) differs from half the profile by \(o(1)\) there. A persistent pointwise deviation would then persist on a fixed smaller interval. Equicontinuity alone does not eliminate the Fourier gap models.

There is also a one-sided envelope criterion without assuming (3). If every limiting lift dominates \(A\,dz\,\delta_0\), or is dominated by it, the marked equation forces equality. For the domination case, the difference is a nonnegative lift with zero right side in (1). Choose a fixed sufficiently broad positive band-limited kernel as in round 40. Its \(K_H\)-averaged complex shift differs from the positive unshifted kernel by a relative error less than \(1/2\) uniformly in \(z,H\). The proof of (1) directly applies to this kernel: its compactly supported Fourier transform is twice continuously differentiable because the kernel has finite absolute second moment. Two integrations by parts suffice for a uniformly summable \((1+|z|)^{-2}\) majorant after averaging in the mark. The real part of the test is positive everywhere, forcing the difference to vanish. The dominated case is unmarked and follows by the same argument.

In particular, either
\(\liminf_{x\to\infty}x^3(f(x)-1)\ge2A\) or
\(\limsup_{x\to\infty}x^3(f(x)-1)\le2A\) suffices for critical-window asymptotics. For the lower envelope, the lift marginal density on each compact window is at least \(A-o(1)\). Cap it at \(A\) while keeping its mark, yielding a submeasure of \(\Lambda_t\). The cap tends locally to \(A\,dz\). Its mass with \(H>\eta\) tends to zero because \(\Lambda_t\) has bounded total mass on the window and its marginal density on \(H>\eta\) is at least \(c_\eta t\). The capped submeasure thus tends to \(A\,dz\,\delta_0\). For the upper envelope, marks vanish uniformly and the lift is asymptotically bounded above by this same unmarked Lebesgue measure. The critical-window law in either case does not preclude thinner exceptional values of \(f\).

## 4. A simple one-sided hypothesis for pointwise asymptotics

Set \(b(x)=x^3(f(x)-1)\). A pointwise limit follows from the one-sided microstep condition
\[
 \eta(T):=\sup_{x\ge T,\ 0\le h\le1/x}(b(x)-b(x+h))_+
                      \longrightarrow0,                       \tag{4}
\]
or the condition with the difference reversed. No bound on a classical derivative is required. Indeed, the round-40 bound gives \(t^2\int_t^{t+1/t}B\le C\). In this interval (4) gives \(b(x)\ge b(t)-\eta(t)\). Using \(tB(x)\ge c\min(b(x),t)\) there proves \(b(t)\) bounded. For the reversed condition use the left interval, applying the one-step inequality with its left endpoint as base.

With bounded \(b\), locally \(tB(x)=b(x)/2+o(1)\) uniformly. Choose \(L(t)\uparrow\infty\) so slowly that \(L(t)=o(t)\) and \(L(t)\eta(t/2)\to0\); if needed use at most \(2L(t)+2\) steps of length between \(1/(2t)\) and \(1/t\). On the right window of length \(L(t)/t\), (4) gives \(b(x)\ge b(t)-o(1)\), and on the left window \(b(x)\le b(t)+o(1)\). The round-40 asymptotic on both windows gives their \(b\)-averages tending to \(2A\), which traps \(b(t)\to2A\). This condition covers eventual monotonicity of \(b\) and is stronger than mere equicontinuity on the critical scale.

## Local verification

The companion Python 3.13.18 script uses SymPy 1.14.0 to check the exponential-mark antiderivative, the shape-2 Gamma distribution function, and the cancellation of the two exact Fourier coefficients at \(k=1/2\). In that example \(M_\infty(1/2)=12/25-16i/25\) and the remaining ordinary marginal amplitude is \(13/25+16i/25\). mpmath 1.3.0 at 80 digits independently evaluates the finite marked transform at four depths and four frequencies, including \(k=0.99\). It also checks the lattice test with successively larger symmetric truncations and calibrates the oscillatory critical-window example. These numerical checks do not establish existence of an exact graph.
