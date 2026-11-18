# The next saturated-origin center deficit coefficient

Use the setting of round 42 and write \(L=\log t\) and \(d=\log\sqrt2\). Let \(D_s=D(\ell_s)\) where \(\ell_s=L-\log(2D_s)\). The raw-center density \(k(t)\) and the delivered density \(m(t)\) of the positive local source satisfy
\[
 \frac{m(t)}{k(t)}
 =
 1+\frac{1}{2D_s}+o(D_s^{-2}).
\]
In particular the coefficient of \(D_s^{-2}\) vanishes. The same proof allows the explicit additional term \(r\log^2\ell/\ell^2\) in \(D\).

In the delivered integral set \(y=tu\), then use \(w=(y-2D)/(2\sqrt D)\) and evaluate \(D=D(L-\log y)\). The normalized integral with respect to \(e^{-w^2}dw/\sqrt\pi\) has amplitude
\[
 A(w)=
 \frac{\sqrt{2D}\exp(D-L-y^2e^{-2L}/2)}
 {1+D' (2+w/\sqrt D)/y}.
\]
At \(w=0\) it equals \(k(t)\) exactly. This follows either from the Jacobian \(dy/dw\) or from the raw pushforward derivative.

Write \(h=D_s^{-1/2}\), \(a=y/(2D_s)-1\), \(\Delta=D-D_s\), and \(Q=D'(\ell_s)\). On \(|w|\le\sqrt{20\log D_s}\), derivatives of the explicit profile give \(Q=1+h^2/2+O(h^4\log D_s)\), \(D''=O(h^4)\). Thus
\[
 a=wh+h^2\Delta+wh(\sqrt{1+h^2\Delta}-1),\qquad
 \Delta=-Q\log(1+a)+O(h^6(1+|w|^2)),
\]
and
\[
 \frac{A(w)}{A(0)}=
 e^\Delta\sqrt{1+h^2\Delta}\,
 \frac{1+Qh^2}{1+\frac{Qh^2}{1+a}
       (1+\frac{wh}{2\sqrt{1+h^2\Delta}})}
 +O(h^5 \log D_s(1+|w|^K))
\]
for a fixed \(K\). The Gaussian source ratio omitted in this formula differs from one by a term exponentially small in \(L\). Substitution gives \(a=wh-wh^3+O(h^5(1+|w|^K))\). The exact formal multiplication in cancellation.py yields
\[
 \frac{A(w)}{A(0)}
 =1-wh+w^2h^2+(w/2-w^3)h^3+
          (w^4-3w^2/2)h^4+o(h^4)
\]
in Gaussian weighted integral. Indeed, the displayed polynomial remainder bound is integrable against the Gaussian and is \(o(h^4)\) after integration. The complement of the \(w\) cutoff inside \(y\in[L/2,6L]\) has negligible Gaussian mass, as \(D\asymp L\), \(dy/dw=O(\sqrt L)\), and \(e^{D-L}=O(L^{-1/2})\) there. The outer \(y\) parts are exponentially small by the round-29 estimate. The odd terms have zero expectation, and \(\mathbb E w^4=(3/2)\mathbb E w^2=3/4\) for this Gaussian weight. This proves the expansion.

Now \(D_s=L-\frac12\log L-d+O(\log L/L)\). Combine with the round-42 first-order expression for \(k(t)\) to get
\[
 k(t)-m(t)=-\frac1{2L}
   -\frac{(j/2+1/8)\log L+q/2-1/2+d/4}{L^2}
   +o(L^{-2}).
\]
Consequently, under exactness and the same exterior control as round 42, with \(L=\log T\),
\[
 \Lambda([0,T])-T
 =-\frac{T}{2L}-\frac{T}{L^2}
     \left[(j/2+1/8)\log L+q/2+d/4\right]
     +o(T/L^2).
\]
The additional \(1/2\) in the bracket is obtained by integrating \(1/(2\log t)\). The earlier exterior errors \(O(1)+O(\log T)\) are negligible on this scale.
