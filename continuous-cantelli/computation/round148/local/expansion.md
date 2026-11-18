# The first two correction terms for the two-well transfer

For the two flanks in project round 144 let \(d_i(m)=x_i(m)-b\) and define the normalized sum
\[
Z(t)=t^2\sum_{i=1}^2\int_M^\infty p_{m/d_i(m)}(t-m)
 \frac{A_i+B_i\cos(2\sqrt m)}{m^2}\,dm,
\]
where \(p_v\) is the centered Gaussian density with variance \(v\). Use \(d_1=1,d_2=2,A_1=A_2=1/2,B_1=1/4,B_2=-e^{-1/4}/4\). Then
\[
Z(t)=1+\frac{13e^{-1/2}}{64\sqrt t}\sin(2\sqrt t)
+\frac{1}{t}\left[\frac34+\frac{169e^{-1/2}}{2048}\cos(2\sqrt t)\right]
+O(t^{-3/2}).
\]
This concerns the sum of the two prescribed flanks, which is the transfer tested in Metal. Other pieces used to complete a continuous graph are separate from this expansion.

To derive it, set \(\epsilon=t^{-1/2}\) and \(m=t+\sqrt t\,w\). With \(d=a-b\) fixed the rescaled kernel is
\[
K_{d,\epsilon}(w)=p_{1/d}(w)(1+\epsilon w)^{-5/2}
 \exp\!\left(\frac{dw^2}{2}\frac{\epsilon w}{1+\epsilon w}\right).
\]
The displaced oscillation phase is \(\phi_\epsilon(w)=2w/[1+\sqrt{1+\epsilon w}]\), so
\[
\phi_\epsilon(w)=w-\epsilon w^2/4+\epsilon^2w^3/8+O_w(\epsilon^3).
\]
The coefficients of the kernel relative to \(p_{1/d}\) are
\[
k_1=-5w/2+dw^3/2,\qquad
k_2=35w^2/8-7dw^4/4+d^2w^6/8 .
\]
On \(|\epsilon w|\le1/2\), the third derivatives of the rescaled kernel and phase product are bounded by a polynomial in \(|w|\) times \(e^{-dw^2/3}\) uniformly at intermediate Taylor parameters. For \(m\notin[t/2,3t/2]\) and \(m\ge M\) the original kernel gives an exponentially small contribution, using its Gaussian drift and the fixed lower cutoff. The actual flank has \(d(m)=d+O(m^{-5/2}e^{-dm/2})\). On the central region this perturbation changes the integral by an exponentially small amount, and the same outside-region bounds apply with \(d(m)\) bounded above and away from zero. Hence Taylor integration yields the stated order of remainder.

The nonoscillating expectation has first correction \(1/(dt)\), from \(\mathbb E k_1=0,\mathbb E k_2=1/d\). In fact for the constant-\(d\) kernel on the untruncated range \(m>0\) its normalized mass is exactly \(1+1/(dt)\). For a direct verification, put \(u(m)=\sqrt d\,(\sqrt m-t/\sqrt m)\) and use the substitution \(m\mapsto t^2/m\). Averaging the original and reflected integrands gives
\[
\int_0^\infty m^{-3/2}e^{-u(m)^2/2}dm
=\frac{1}{t\sqrt d}\int_{\mathbb R}e^{-u^2/2}du .
\]
Multiplying by \(e^{-dt}\) and differentiating in \(t\) gives the corresponding \(m^{-5/2}\) identity and the exact normalized mass.

For the complex oscillating coefficient, let \(v=1/d\) and \(\vartheta=2\sqrt t\). The \(\epsilon\) and \(\epsilon^2\) integrands against \(e^{iw}p_v(w)\) have polynomial factors
\[
k_1-iw^2/4,\qquad
k_2-w^4/32+i(3w^3/4-dw^5/8).
\]
Differentiating the Gaussian characteristic function gives the integrated coefficients
\[
-i e^{-v/2}(5v+v^2)/4,\qquad
e^{-v/2}(v+v^2/32-3v^3/16-v^4/32).
\]
The \(B_i\)-weighted constants cancel. The \(B_i\)-weighted first coefficient is \(-i\,13e^{-1/2}/64\) and the weighted second is \(169e^{-1/2}/2048\). Taking the real part after multiplication by \(e^{i\vartheta}\) gives the expansion.

In particular, on the phase sequences \(t=\pi^2n^2\) and \(t=\pi^2(n+1/2)^2\), the first sine vanishes exactly. The scaled next residuals \(t(Z(t)-1)\) tend to \(3/4\pm169e^{-1/2}/2048\), both strictly positive. On the two quarter-turn phases the scaled first residuals converge to opposite nonzero limits. Thus the cancellation that gives constant leading total flux leaves a genuine first oscillation and a positive second correction on each zero-sine phase.
