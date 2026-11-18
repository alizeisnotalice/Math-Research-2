# Eventually nonnegative tail excess and local deficit mass

Let \(0\le f\le M\) be continuous, with the exact identity
\[
 F_f(z)=\int e^{zx-z^2f(x)/2}\gamma(dx)=1 .
\]
Assume on the positive side that \(f\ge1\) for \(x\ge R>0\), that \(f>0\) on \((0,R]\), and that \(f(x)=cx^2+o(x^2)\) as \(x\to0\) on both sides, with \(c>0\). The real identity extends to the whole complex plane by boundedness and the identity theorem. Write \(\nu_t(dx)=e^{tx-t^2f(x)/2}\gamma(dx)\), \(\phi=p_1\), and
\[
 \begin{split}
 A&=e^{1/(2c)}/\sqrt c,\\
 B(x)&=1-f(x)^{-1/2}e^{-x^2(1-1/f(x))/2},\quad x\ge R,\\
 D_t&=\nu_t(x<R)-\gamma_t(x<R)
     =\int_{x\ge R} \bigl(1-e^{-t^2(f(x)-1)/2}\bigr)\phi(x-t)\,dx ,
 \end{split}
\]
where \(\gamma_t=N(t,1)\). The quadratic substitution at the origin and compact positivity give \(D_t\sim A/t\). The proof below also works with two different positive quadratic coefficients, with \(A\) equal to the sum of the corresponding half-line Gaussian integrals.

## 1. Quantitative local law

For every fixed \(K>0\) there exists \(C_K<\infty\) such that, for every sufficiently large \(t\) and every interval \(I\subset[-K,K]\),
\[
 \left|\int_{t+I}B(x)\,dx-D_t|I|\right|\le {C_K\over t^2}.        \tag{1}
\]
This bound is additive and uniform even for intervals shorter than \(1/t\).

In particular, uniformly over \(b>a\to\infty\) with \(a(b-a)\to\infty\),
\[
 \int_a^b B(x)\,dx\sim A\log(b/a).                               \tag{2}
\]
More precisely, put \(\epsilon(a)=\sup_{s\ge a}|sD_s/A-1|\to0\). Equation (1) yields a relative error at most
\(\epsilon(a)+C/[a\min(1,b-a)]\) for (2). To see this, use \(t=(a+b)/2\) for windows of length at most 1, and partition longer windows into equal lengths in \([1/2,1]\), summing the additive errors.

Every such short window contains a point \(x\) with \(x^3(f(x)-1)\to2A\). Indeed continuity realizes its average of \(B\) at a point. Since \(f\ge1\), the relation \(B(x)\sim A/x\) implies \(f(x)-1\sim2A/x^3\) by the definition of \(B\). For example, this point can always be found within distance \(L(t)/t\) of \(t\), for any \(L(t)\to\infty\) with \(L(t)=o(t)\). If \(f\) is eventually monotone, use neighboring windows on the two sides to trap \(f(t)\) and conclude the full pointwise limit \(t^3(f(t)-1)\to2A\), without assuming a tail expansion.

For any fixed \(q>0\), the Lebesgue measure of \(\{x\in t+I:f(x)-1\ge q/t^2\}\) is
\(O_{q,K}(|I|/t+t^{-2})\), since \(B\) has a fixed positive lower bound on this high set. An entire interval at that height must have length \(O_{q,K}(t^{-2})\).

## 2. Positive marked deficit and its compact-controlled transform

Let \(u=x-t\), and use \(\mu_t(du)=\phi(u)(1-e^{-t^2(f(t+u)-1)/2})1_{t+u\ge R}\,du\), of total mass \(D_t\). This is the marginal of the positive measure
\[
 d\omega_t(u,V)=\phi(u)\tfrac12 e^{-V/2}
  {\bf1}_{t+u\ge R,\ 0\le V\le t^2(f(t+u)-1)}\,du\,dV .
\]
The mark, conditional on a source with positive marginal density, is an exponential of rate \(1/2\) conditioned to be below the displayed cutoff. Let \(z=t+is\) and \(W=x-tf(x)\). Exactness and the integral of the exponential mark give
\[
 \begin{split}
 Q_t(s)&=\int e^{is(u-V/t)+s^2V/(2t^2)}d\omega_t\\
 &= {t^2\over z^2}
 \left[e^{-s^2/2}\int_{x<R} e^{isW+s^2 f(x)/2}\,d\nu_t
                   -\int_{x<R}e^{is(x-t)}d\gamma_t\right].       \tag{3}
 \end{split}
\]
Here the factor \(t^2/z^2\) comes from integrating \(\frac12 e^{-z^2V/(2t^2)}\) in \(V\), and the complementary integrals follow by subtracting the \(f=1\) identity at \(z\).

Uniformly for \(|s|\le t/2\),
\[
 |Q_t(s)-D_t e^{-s^2/2}|
 \le C\bigl(|s|t^{-2}+s^2t^{-3}+e^{-c_0t}\bigr)e^{-s^2/2}
              +C\gamma_t(x<R) .                                \tag{4}
\]
On a fixed small origin neighborhood the substitution \(y=tx\), the bounds \(c_1x^2\le f(x)\le c_2x^2\), and \(t^2-s^2\ge3t^2/4\) give uniform marked integrals \(O(t^{-2})\) for \(|W|\) and \(O(t^{-3})\) for \(f(x)\) against \(e^{s^2f(x)/2}\nu_t\). Use
\(\lvert e^{isW+s^2f/2}-1\rvert\le e^{s^2f/2}(|s||W|+s^2f/2)\).
The remaining positive compact has a positive minimum, and negative sources beyond the origin neighborhood contribute \(O(e^{-c_0t})\), in each case even after the marked weight. The term \(t^2/z^2-1\) is \(O(|s|/t)\), and the total compact mass is \(O(1/t)\). The baseline term is at most \(\gamma_t(x<R)\), explaining (4).

## 3. Positive band-limited smoothing

Let \(\psi\) be the normalized even positive kernel proportional to
\[
 \operatorname{sinc}^4 y+\operatorname{sinc}^4(y+\pi/2)
                        +\operatorname{sinc}^4(y-\pi/2),
 \qquad \operatorname{sinc}(y)=\sin(y)/y.
\]
It has Fourier support \([-4,4]\), because \(\operatorname{sinc}\) is the characteristic function of a uniform variable on \([-1,1]\). On the real line \(\psi(y)\asymp(1+|y|)^{-4}\), by the sine and cosine terms and bounded shifts of their denominators. It is entire with
\[
 |\psi(a+ib)|+|\psi'(a+ib)|\le C e^{4|b|}(1+|a|)^{-4}.
\]
For bounded \(|a|\), the bound for sinc and its derivative follows from the uniform integral representation. For large \(|a|\), use \( |\sin(a+ib)|+|\cos(a+ib)|\le C e^{|b|}\) and \( |a+ib|\ge |a|\) in the derivative formula. Thus for all real \(y\) and complex \(w\),
\[
 \begin{split}
 |\psi(y+w)|&\le C(1+|w|)^4 e^{4|\operatorname{Im}w|}\psi(y),\\
 |\psi(y+w)-\psi(y)|
   &\le C|w|(1+|w|)^4 e^{4|\operatorname{Im}w|}\psi(y).           \tag{5}
 \end{split}
\]

Fix a sufficiently large constant \(L\), set \(h=L/t\), and put \(\psi_h(y)=h^{-1}\psi(y/h)\).
Let \(G\sim N(0,1)\) be independent of the mark, and let the complex shift be
\[
 w=(-V+i\sqrt V G)/t .
\]
For \(j=0,1\), uniform conditional moments at every admissible source give
\[
 \mathbb E\!\left[(|w|/h)^j(1+|w|/h)^4
                       e^{4|\operatorname{Im}w|/h}\mid u\right]
                     \le C_j L^{-j}.                           \tag{6}
\]
To check this, integrate in \(G\) using a Gaussian exponential moment to bound the integrand by \(C L^{-j}(1+V)^5 e^{C V/L^2}\) before integrating in \(V\). The upper bound is increasing in \(V\); conditioning an exponential below any positive cutoff only decreases the expectation. Its integral over the full exponential is bounded uniformly for \(L\) large.

Set \(P_h(a)=\int\psi_h(a-u)\,d\mu_t\) and \(P_h^*(a)=\int\mathbb E_G\psi_h(a-u-w)\,d\omega_t\).
Fourier inversion of the latter is justified by compact Fourier support and the finite marked integrals. Its transformed measure is \(Q_t\) from (3). As \(4/h\le t/2\), integrating (4) over the Fourier support gives
\[
 P_h^*(a)=D_t(\psi_h*\phi)(a)+O(t^{-2})                          \tag{7}
\]
uniformly for \(a\in\mathbb R\). From the difference bound in (5) and (6) with \(j=1\), \(|P_h^*(a)-P_h(a)|\le(CL^{-1})P_h(a)\). Choose \(L\) sufficiently large once and for all. Then (7) and \(D_t=O(t^{-1})\) imply \(P_h(a)\le C/t\) on the full line. Positivity and the lower bound for \(\psi\) on \([-1,1]\) give
\[
 \mu_t([a-h,a+h])\le C h/t,\qquad a\in\mathbb R.                  \tag{8}
\]

## 4. Endpoint errors and transfer to \(B\)

For an interval \(I=[a,b]\), put \(\Psi_I={\bf1}_I*\psi_h\) with the same entire extension. The derivative is \(\psi_h(\cdot-a)-\psi_h(\cdot-b)\). Integrate this derivative along the mark displacement and use the first bound in (5) inside the segment. By (6), with \(|w|=h|w|/h\), this yields the mark replacement estimate
\[
 \left|\int\mathbb E_G\Psi_I(u+w)\,d\omega_t
                  -\int\Psi_I(u)\,d\mu_t\right|
 \le {C\over t}\bigl(P_h(a)+P_h(b)\bigr)=O(t^{-2}).               \tag{9}
\]
The right side is independent of the interval length. Since \(\psi=O((1+|y|)^{-4})\), its one-sided tail is \(O((1+|y|)^{-3})\). Using (8) on strips of width \(h\) near each boundary shows
\[
 \left|\mu_t(I)-\int \Psi_I\,d\mu_t\right|=O(h/t)=O(t^{-2}).
\]
The same bound holds when \(\mu_t\) is replaced by \(D_t\phi(u)du\), by its bounded density. Integrating (7) over \(a\in I\), then using (9) and these two leakage estimates, proves
\[
 \left|\mu_t(I)-D_t\int_I\phi(u)\,du\right|
                            \le C(1+|I|)t^{-2}.                \tag{10}
\]
On \(|u|\le K\) multiply this measure comparison by \(1/\phi(u)\). Stieltjes integration by parts on every subinterval bounds the new error by the uniform discrepancy in (10) times the supremum and total variation of \(1/\phi\) on \([-K,K]\). Writing \(g=f-1\) and \(z=t^2g/2\), one has for \(x=t+u\), \(|u|\le K\),
\[
 B(x)=1-f(x)^{-1/2}\exp[-z x^2/(t^2 f(x))],\qquad
 |B(x)-(1-e^{-z})|\le {C_K\over t}(1-e^{-z}).
\]
For the latter, \(x^2/(t^2f)\ge1/(4M)\) and its difference from 1 is \(O_K(t^{-1})+O(z/t^2)\). The mean value theorem on the exponential, together with \(1-f^{-1/2}\le g/2\), reduces to \(z^j e^{-c z}\le C(1-e^{-z})\) for \(j=1,2\). Using (10) once more for the integrated mass proves (1).

## Resolution and limitations

The additive error permits high and narrow deviations on the critical scale. As a synthetic illustration of what (1) and (2) alone control, \(B(x)=A x^{-1}(1+\frac12\sin x^2)\) obeys an \(O(t^{-2})\) additive window law but has no pointwise limit after multiplication by \(x\). For \(t=\sqrt{2\pi n}\), \(I=[0,d/t]\), its excess over \(A\log(1+d/t^2)\), multiplied by \(t^2\), tends to \(A(1-\cos(2d))/4\). This calculation concerns only the window estimates. Whether the full exact identity gives relative asymptotics at critical or shorter scales, or unrestricted pointwise asymptotics, remains open in this argument.

## Executed local checks

The companion script checks the marked exponent and endpoint derivative in SymPy 1.14.0. In NumPy 2.5.4 and SciPy 1.18.1 it uses 72-point Gauss–Legendre and 48-point Gauss–Hermite quadrature at four truncated-exponential cutoffs and ten real kernel locations. Sample maxima for the relative marked-kernel error at \(L=8,16,32,64\) were \(0.7584,0.3448,0.1623,0.07856\). The experiment only calibrates the analytically uniform \(O(1/L)\) bound. In mpmath 1.3.0 at 80 digits, the quadratic Gaussian integral tested the compact Fourier error and the explicit synthetic critical-window residual. The latter was \(0.499987500178928716\) at \(n=10{,}000\), compared with its analytic limit \(1/2\) for \(A=1,d=\pi/2\). Python 3.13.18 completed the script.
