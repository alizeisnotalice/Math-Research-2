# Two phase thresholds for critical logarithmic oscillations

Let \(f\ge0\) be continuous, positive off zero, with \(\liminf_{|x|\to\infty} f(x)>1\). Suppose it has an even core and put \(\ell=\log(1/|x|)\). Assume throughout that, for all sufficiently small nonzero \(|x|\),
\[
 \frac{x^2}{2f(x)}=D(\ell)=\ell+\tfrac12\log\ell+c+a\cos(\omega\ell^\beta+\varphi),
 \qquad \beta,\omega>0 .
\]
Here \(c,a,\varphi\) are fixed real numbers. Set \(F(z)=\mathbb E e^{zX-z^2 f(X)/2}\) for standard Gaussian \(X\), and write \(L=\log z\), \(\ell_0=L-\log(2L)\), \(\Theta(L)=\omega\ell_0^\beta+\varphi\). The two relevant observables are \(F(e^L)\) and \(A(L)=\mathbb E_G F(e^LG)\) with an independent standard Gaussian \(G\).

## Point transform

For \(x>0\) use \(y=e^Lx\), then \(y=2L+\sqrt{2L}\,w\). The exponent from the tilt is exactly
\[
 y-\frac{y^2}{4D}=D-\frac{(y-2D)^2}{4D}.
\]
The portion \(x<0\) is \(O(e^{-L})\). Compact annuli away from zero contribute exponentially little because \(f\) has a positive minimum there. On the far tails use \(f\ge1+\delta\) and complete the square against the Gaussian source, giving an exponentially small bound in \(e^{2L}\).

For \(y\le L/2\), \(y-L\le-L/2\) bounds the exponent including the change of variables. For \(L/2\le y\le 4L\), \(D=L-\frac12\log L+O(1)\) uniformly. Hence the integrand in \(w\) is dominated by \(C e^{-c_0w^2}\), and on each bounded \(w\)-interval it differs uniformly by \(o(1)\) from
\[
 \frac{e^c}{\sqrt2}\phi(w)
 \exp\!\left[a\cos\!\left(\omega
        \big(\ell_0-\log(1+w/\sqrt{2L})\big)^\beta+\varphi\right)\right].
\]
This approximation does not expand the phase. For \(y\ge4L\) still in the core, the monotone baseline of \(D\) and bounded modulation yield \(D\le L-\frac12\log L+C\), so the Gaussian penalty absorbs the whole tail. Together these estimates justify replacing the integral by that of the displayed density with \(o(1)\) error.

For \(\beta<3/2\), the phase difference from \(\Theta(L)\) vanishes uniformly on bounded intervals. At \(\beta=3/2\), it converges to \(-\kappa w\) with \(\kappa=3\omega/(2\sqrt2)\). For \(\beta>3/2\), the derivative of the phase in \(w\) has magnitude comparable to \(L^{\beta-3/2}\), and its second derivative divided by the squared first derivative is \(O(L^{1-\beta})\) on every bounded interval. Integrate each nonzero Fourier harmonic by parts against a smooth compactly supported approximation to the Gaussian weight. Then remove the compact support by the integrable bound. The Fourier series of \(\exp(a\cos\theta)\) is absolutely convergent. Consequently, to \(o(1)\), \(F(e^L)\) equals
\[
 \frac{e^c}{\sqrt2}
 \begin{cases}
 e^{a\cos\Theta(L)}, & \beta<3/2,\\
 \mathbb E e^{a\cos(\Theta(L)+\kappa G)}, & \beta=3/2,\\
 I_0(a), & \beta>3/2 .
 \end{cases}
\]
The critical multiplier for the \(n\)-th harmonic is \(e^{-n^2\kappa^2/2}\). If \(a\ne0\), the first harmonic is nonzero, and \(\Theta(L)\) is continuous, strictly increasing and unbounded for large \(L\). Thus \(F\equiv1\) excludes these profiles for \(\beta\le3/2\). For \(\beta>3/2\) this leading calculation gives the necessary normalization \(c=\frac12\log 2-\log I_0(a)\).

More generally, replacing \(\frac12\log\ell\) by \(b\log\ell\) multiplies the displayed asymptotics by \(L^{b-1/2}\), with the same error in relative units. Since the bounded modulation keeps the remaining factor bounded above and away from zero, exactness forces \(b=1/2\) in this family.

## Multiplicative Gaussian average

Integrating first in \(G\), justified by Tonelli, gives with \(\lambda=e^{-2L}\)
\[
 A(L)=\int_{\mathbb R}\phi(x)
   \sqrt{\frac{\lambda}{\lambda+f(x)}}\,
   \exp\!\left(\frac{x^2}{2(\lambda+f(x))}\right)dx .
\]
The exterior contributes \(O(e^{-L})\) under the same assumptions. On the even core use \(y=e^L|x|=2L e^s\), \(\ell=\ell_0-s\) and \(H=1+y^2/(2D)\). The integrand in \(s\) becomes
\[
 \sqrt{\frac{2\ell}{\pi H}}\,
 \exp\!\left(c+a\cos(\omega(\ell_0-s)^\beta+\varphi)
                  -D/H-\tfrac12y^2e^{-2L}\right).
\]
Locally in \(s\), this differs by \(o(1)\) from \((e^c/\sqrt2)k(-s)\) times the same periodic factor, where \(k(v)=\sqrt{2/\pi}e^{v-e^{2v}/2}\) is the density of \(\log |G|\). For positive \(s\) the integrand is bounded by \(Ce^{-s}\), since \(D\asymp\ell\) and \(H\ge y^2/(2D)\). For negative \(s\) with \(y^2\ge L\), \(\ell\asymp L\) and \(D/H\ge c_0 e^{-2s}\), so \(C e^{-s-c_0e^{-2s}}\) dominates. When \(y^2<L\), \(\ell\ge L-O(\log L)\) and \(D/H\ge c_0\ell\); integrating the bound \(C\sqrt\ell e^{-c_0\ell}\) gives an exponentially small contribution. This proves the asymptotic replacement on the full interval.

For \(\beta<1\) the phase displacement tends to zero. At \(\beta=1\) it is exactly \(-\omega s\). At \(\beta>1\) the phase derivative in \(s\) tends uniformly to infinity in absolute value on compact sets and its second derivative divided by its first derivative squared is \(O(L^{-\beta})\). The same Fourier integration-by-parts argument yields, to \(o(1)\),
\[
 A(L)=\frac{e^c}{\sqrt2}
 \begin{cases}
 e^{a\cos\Theta(L)}, & \beta<1,\\
 \mathbb E e^{a\cos(\Theta(L)+\omega\log|G|)}, & \beta=1,\\
 I_0(a), & \beta>1 .
 \end{cases}
\]
At the critical frequency the multiplier on the \(n\)-th harmonic is
\[
 m(n\omega)=\mathbb E|G|^{in\omega}
       =2^{in\omega/2}\frac{\Gamma((1+in\omega)/2)}{\sqrt\pi},
 \qquad |m(u)|=\cosh(\pi u/2)^{-1/2}>0 .
\]
Thus Gaussian averaging at large scales erases the oscillation for every \(1<\beta<3/2\), while the point transform still oscillates. This compares only the large-scale asymptotics of these two quantities. Requiring the entire averaged identity \(A(L)=1\) at all scales imposes more than its asymptotic limit.

## Numerical scope

The 12-process scan used deterministic parameters, paired composite Gaussian and adaptive quadrature, and direct reduction of the Gaussian mixture to one integral. For the point transform it integrates \(-9\le w\le 11\), and for the average it integrates \(-4\le s\le23\). All sampled integrals are in the even core geometry on these windows. The exterior estimates above account for fixed positive continuations at large scales. The scan uses \(c=\frac12\log 2-\log I_0(a)\) to normalize the common zero harmonic. It caps the sampled point local angular frequencies at 55 and evaluates the Gaussian average only when its central angular frequency is at most 105. One initial run stopped on a coarse-grid convergence assertion, and the complete run used finer frequency-dependent panels.

## Next-order obstruction above the point threshold

Suppose \(\beta>3/2\). Let \(B\) solve \(B=\ell_b+\frac12\log\ell_b+c\), where \(\ell_b=L-\log(2B)\). The equation has a unique solution near \(L\) for large \(L\) because the left side increases and the right side decreases there. In particular \(B=L+O(\log L)\). Recenter the point transform at \(y=2B+\sqrt{2B}\,w\), and let \(\Psi_L(w)=\omega(\ell_b-\log(1+w/\sqrt{2B}))^\beta+\varphi\). With \(\delta=D-B\), the exponent after removing \(B-L\) is
\[
 \delta-\frac{(\sqrt{2B}\,w-2\delta)^2}{4(B+\delta)}.
\]
On the symmetric window \(|w|\le W=\sqrt{C\log B}\),
\(\delta=a\cos\Psi_L(w)-w/\sqrt{2B}+O((1+w^2)/B)\). Taylor's formula thus gives the integrand, relative to the prefactor \(\sqrt{2B}e^{B-L}\), as
\[
 \phi(w)e^{a\cos\Psi_L(w)}
 \left[1+\frac{w}{\sqrt B}
   \left(\sqrt2\,a\cos\Psi_L(w)-\frac1{\sqrt2}\right)
      +O((1+|w|^6)/B)\right].
\]
The error integrates to \(O(1/B)\), and the Gaussian bounds from the first section make the discarded tails \(O(B^{-2})\) by enlarging \(C\). On a slightly larger window, \(|\Psi_L'|\asymp B^{\beta-3/2}\) and for each fixed \(j\ge2\), \(\Psi_L^{(j)}/\Psi_L'=O(B^{-(j-1)/2})\). Repeated integration by parts with an even smooth cutoff therefore bounds the integral of each nonzero harmonic against \(\phi(w)\) or \(w\phi(w)\) by \(C_N(|n|B^{\beta-3/2})^{-N}\). The relevant Fourier series are absolutely convergent, and \(N\) can be chosen so that \(N(\beta-3/2)>1\). The zero mode of the correction term is odd and integrates to zero. We obtain
\[
 F(e^L)=\sqrt{2B}e^{B-L}\big(I_0(a)+O(B^{-1})\big)
 =\frac{e^c I_0(a)}{\sqrt2}
       \left[1-\frac{\log L}{4L}+O(L^{-1})\right].
\]
The second equality uses the exact identity \(\sqrt{2B}e^{B-L}=(e^c/\sqrt2)\sqrt{\ell_b/B}\). Thus the stated profile is excluded even for \(\beta>3/2\).

For a profile with baseline correction \(d(\log\ell)/\ell+o((\log\ell)/\ell)\), incorporate the displayed \(d\)-term into the defining equation for \(B\). This adds the factor \(\exp(d\log\ell_b/\ell_b)\) to the exact prefactor above and leaves the first correction term and the frequency estimates unchanged to the same order. The little-oh changes the integral by \(o((\log L)/L)\) because the derivative of the point exponent with respect to \(D\) is \(y^2/(4D^2)=O(1)\) in the central window. Exactness therefore requires \(d=1/4\) in this family.

The large-scale limitation for \(A\) is distinct from the exact condition at all scales. By Laplace transform uniqueness, the full identity \(A(L)=1\) for every real \(L\) implies \(F(z)+F(-z)=2\) for every real \(z\). Indeed change variables \(v=z^2\) inside the Gaussian average and use the uniqueness of the Laplace transform on the signed density proportional to \(v^{-1/2}(F(\sqrt v)+F(-\sqrt v)-2)\). When \(f\) is globally even, \(F\) is even, and the two full identities are equivalent.
