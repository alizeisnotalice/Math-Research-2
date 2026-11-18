# Past support, record jumps, and a logarithmic regularity threshold

Suppose \(f\ge0\) is bounded and continuous on \(\mathbb R\), with
\[
 \int e^{zx-z^2 f(x)/2}\gamma(dx)=1,\qquad z\in\mathbb C.
\]
Assume both tails have \(\liminf f\ge1+T\) for \(T>0\). Set \(H=\sup f\) and
\[
 Q_t=P_t\delta_0-\int_{f(x)\le t}P_{t-f(x)}\delta_x\,\gamma(dx).
\]
We use no positivity of \(Q_t\). At each early time the past set is compact.

## 1. Early support and the sign of the observation point

For \(y\ge0\), the one-sided signed-resolvent identity is
\[
 \int_0^H e^{-qt}Q_t(y)\,dt
 =
 \frac1{\sqrt{2q}}\int_{x>y}e^{-q f(x)}
 \left(e^{\sqrt{2q}(x-y)}-e^{-\sqrt{2q}(x-y)}\right)\gamma(dx).
                                                               \tag{1}
\]
It follows from the heat resolvent and the exact real identity at \(z=\sqrt{2q}\).
The time density in (1) is integrable for each fixed \(y\), since the time integral of any Gaussian heat kernel over a bounded horizon is uniformly bounded. The zero-age spatial layer at \(t=f(y)\) is on a single time point and does not affect this formulation.

If \(0<t<T\) and \(y>\max(0,\sup\{f\le t\})\), then \(\inf_{x\in[y,R]} f(x)>t\) on every nonempty compact interval of that form. The far tail has \(f\ge1+t+\delta\) for some \(\delta>0\). Hence the transform in (1) decays exponentially at a rate strictly greater than \(t\), estimating the far integral using \(\int e^{\sqrt{2q}x}\gamma(dx)=e^q\). Signed compact-time Laplace support (as in the previous radial proof, using the order-one quadrant argument) gives \(Q_t(y)=0\). First this holds for almost every time below the decay rate at any fixed \(y\), and then pointwise when \(t<f(y)\) by local continuity away from the graph. Reflect to the left. Thus in full generality the direct conclusion is
\[
 \operatorname{supp}Q_t\subseteq
        \operatorname{conv}\bigl(\{0\}\cup\{f\le t\}\bigr).
                                                               \tag{2}
\]
The formula does not permit the observation point \(y<0\) on the right, because the free resolvent is \(e^{-\sqrt{2q}|y|}/\sqrt{2q}\). When \(f(0)=0\), (2) is exactly the past convex hull.

## 2. A jump forces exponentially thin new sublevels

Assume \(f(0)=0\), write \(B(t)=\max\{f\le t\}\), and suppose for \(t_0\in(0,T)\) that
\[
 b=B(t_0-)<x_0=B(t_0),\qquad d=x_0-b>0.
\]
For \(x>b\), no \(f(x)<t_0\) is possible by the definition of \(b\). In a neighborhood of \(x_0\) inside \(x>b\) we have \(f\ge t_0\), so \(x_0\) is a (possibly non-strict) local minimum. The strict-past measure
\[
 Q_{t_0-}=P_{t_0}\delta_0-
          \int_{f(x)<t_0}P_{t_0-f(x)}\delta_x\,\gamma(dx)
\]
is supported on \((-\infty,b]\) by taking the weak limit from earlier times in (2). The limit follows by dominated convergence on bounded continuous tests, including the change from \(\le\) to \(<\). The endpoint \(B(t)\) is right continuous: intersect a decreasing sequence of compact sublevels in one common compact region, which is available by the tail assumption.

We first bound the old positive part better than by total variation. For \(x<b\), reflect \(x\) to \(x^*=2b-x>b\). Every strict-past mean lies at \(z\le b\), so \(|x-z|\le|x^*-z|\). The positive past potential is therefore larger at \(x\) than at \(x^*\). At \(x^*\) it equals \(p_{t_0}(x^*)\) by the old support. The identity holds pointwise there because all these means are separated from \(x^*\), so the past potential has a continuous density. Thus, for almost every \(x<b\),
\[
 (Q_{t_0-})^+(x)\le
    \big(p_{t_0}(x)-p_{t_0}(2b-x)\big)_+
       \le C_{t_0}(b-x).                                     \tag{3}
\]
All strict-past ages are positive, so the old measure is absolutely continuous. The argument uses only nonnegativity of the past measure, not of its signed complement.

Fix \(0<r<d/4\) and set \(y=x_0+d-2r\). Then
\[
 y-b=2(d-r),\qquad
 y-(x_0-r)=d-r,\qquad y\ge x_0+d/2 .
\]
By right continuity, for all sufficiently small \(h>0\), uniformly over these \(r\), one has \(y>B(t_0+2h)\). The semigroup identity and late support give the pointwise equality
\[
 P_{2h}Q_{t_0-}(y)=
  \int_{t_0\le f(x)\le t_0+2h}
               p_{t_0+2h-f(x)}(y-x)\,\gamma(dx).               \tag{4}
\]
The possible zero-age layer is separated from \(y\), so causes no difficulty. By (3),
\[
 P_{2h}Q_{t_0-}(y)
 \le C\int_0^\infty u\,p_{2h}(y-b+u)\,du
 \le C' h^{3/2}\exp[-(y-b)^2/(4h)].
\]
The last step uses the factor \(e^{-(y-b)u/(2h)}\) and \(\int_0^\infty ue^{-c u/h}du=h^2/c^2\), with \(y-b\) bounded below. On
\[
 E_{r,h}=[x_0-r,x_0+r]\cap\{f\le t_0+h\}
\]
we have \(f\ge t_0\), age in \([h,2h]\), and \(0<y-x\le d-r\). Thus the corresponding kernel in (4) is at least \((4\pi h)^{-1/2}\exp[-(d-r)^2/(2h)]\). Combining,
\[
 \gamma(E_{r,h})\le C h^2\exp[-(d-r)^2/(2h)],
       \qquad 0<r<d/4,\quad 0<h<h_0 .                        \tag{5}
\]
Both constants may depend on the fixed jump, but not on \(r,h\). In particular there is no positive-length plateau at that level in a sufficiently small neighborhood of the new extreme.

Let \(m(r)=\max (f(x)-t_0)\) on either one-sided interval of radius \(r\) from \(x_0\), or on the two-sided interval. This is nonnegative there and cannot be zero on a positive-length interval by (5). Taking \(h=m(r)\) in (5) and using a positive lower Gaussian density near \(x_0\), with \(L=\log(1/r)\), gives
\[
 L\ge\frac{(d-r)^2}{2m(r)}+2\log(1/m(r))-C.                  \tag{6}
\]
In particular \(\liminf_{r\downarrow0} Lm(r)\ge d^2/2\). More sharply,
\[
 L-\frac{d^2/2}{m(r)}-2\log L\ge -C'
       \qquad(r\text{ sufficiently small}).                  \tag{7}
\]
Indeed (6) gives \(m(r)\gtrsim1/L\), so the replacement of \(d-r\) by \(d\) costs \(O(r/m(r))=o(1)\). If \(m(r)\le d^2/L\), then \(\log(1/m(r))\ge\log L-O(1)\), and otherwise the left side of (7) tends to \(+\infty\). The reflected statement holds for left jumps. This refines the \(d^2/2\) logarithmic threshold by two log-log units. The exponential coefficient arises from the optimal split: if the new level excess is \(h=\theta\delta\), the exponent difference in the mass bound has numerator \(\theta[(y-b)^2-(y-x_0)^2/(1-\theta)]/2\), maximized at \(y-x_0=d(1-\theta)/\theta\), with value \(d^2/2\).

## 3. Connecting the past and applying the full-tail budget

A sufficient assumption for eliminating all early jumps and all nonzero zeros is the following local condition: at every nonzero local minimum \(x\) with \(f(x)<T\), one of the two one-sided upper moduli
\[
 m_x^\pm(r)=\sup_{0\le u\le r}\big(f(x\pm u)-f(x)\big)_+
\]
has \(\liminf_{r\downarrow0}m_x^\pm(r)\log(1/r)=0\). In particular local Hölder continuity is more than sufficient. A flat local minimum counts as a local minimum here. At a zero \(x\ne0\), the elementary real-parameter test over that side gives
\(\liminf m_x^\pm(r)\log(1/r)\ge |x|^2/2\), since that interval lies at absolute distance \(|x|-O(r)\) from zero and optimizing the real exponent gives the square-distance estimate of round 6. There must be at least one zero: otherwise compact positivity and the two high tails make \(F(\lambda)\to0\) as \(\lambda\to+\infty\). Thus the only zero is \(0\).

Now every positive record jump would occur at a nonzero local minimum, and (6) contradicts the stated modulus there. Reflect for the negative endpoint. Both endpoints of the early sublevel hull are continuous on \(0<t<T\), and tend to zero as \(t\downarrow0\) by unique zero and compactness. For any \(x\in[0,B(t)]\), continuity of the nondecreasing endpoint shows that \(x=B(s)\) for some \(s\in[0,t]\), hence \(f(x)\le t\). The same holds for the left wing. The actual early sublevel is therefore the full interval, not merely contained in it.

Round 84 now applies, giving the global excess budget. In particular if \(f(x)\to H=\sup f>1\) on both tails, take \(T=H-1\). Since \(f\le1+T\) globally, the excess budget is impossible. This excludes all such local-Hölder exact graphs, independently of the shape of either wing at or above level \(T\). A positive uniform logarithmic slope \(L_0\) by itself instead bounds a jump by \(\sqrt{2L_0}\); without a lower bound on the possible jump length it does not force connectedness.

## Local verification

The companion Python 3.13.18 script uses SymPy 1.14.0 to verify the shift optimization \(d(1-\theta)/\theta\), its maximal coefficient \(d^2/2\), and the finite-radius exponent \((d-r)^2/2\). At 70 decimal digits, mpmath 1.3.0 computes the reflected-wall kernel model. Its old heat term divided by \(h^{3/2}e^{-(y-b)^2/(4h)}\) is \(0.0550763771294622\) at \(h=2^{-20}\), approaching the analytic value \(0.0550764367040723\). It separately verifies the matching order \(h^2e^{-A/h}\asymp r\) for the model \(h=A/[\log(1/r)-2\log\log(1/r)+C]\). These calculations check the comparison scales and do not construct a globally exact graph.
