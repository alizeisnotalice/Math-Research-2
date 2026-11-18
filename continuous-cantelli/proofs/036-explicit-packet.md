# An explicit fourth-order Gaussian-packet coefficient

In the notation of round 27 with \(n=4\), take
\[
A=3/8,\qquad g=\sqrt3/2,\qquad
w_*=-1/4+i\sqrt3/4,\qquad
D=2cd,\qquad D_*=1/(4\sqrt2).
\]
Let \(\delta=79/100\), and write
\[
G(z)=iz-z^4/4+
\frac{\delta(iz-z^4/2-w_*)^2}{2(1+\delta z^4/2)}.
\]
An explicit sufficient coefficient in the fourth-order case is
\[
  D\ge 149/1000,\qquad\text{equivalently}\qquad
  d\ge 149/(2000c).
\]
This is about \(84.3\%\) of the original critical coefficient.

The real error exponent is controlled by
\[
 \Re G(z)=
 \frac{-z^4/4+\delta[-z^4/8+1/32-(z-g/2)^2/2]}
 {1+\delta z^4/2}\qquad(z\in\mathbb R).
\]
At \(D=149/1000\), multiply \(D-z^2(A+\Re G(z))\) by its positive denominator and substitute \(z=g x\). The result is rational of degree six, positive at zero, and has zero real roots by an exact Sturm check. Thus \(-D/z^2+\Re G(z)<-A\) for real \(z\ne0\), with a uniform strict gap on the bounded input interval by continuity and divergence to \(-\infty\) at zero. Increasing \(D\) only improves this bound.

For the leading contour, \(\operatorname{Re}(z^4)\ge-8(\operatorname{Im}z)^4\ge-1/2\) throughout \(0\le\operatorname{Im}z\le1/2\). Thus the denominator has positive real part there. On \(z=g x+i/2\), clear the positive denominator \(|1+\delta z^4/2|^2\) in \(\Re G+A\). The numerator factors as \((x-1)^2\) times a strictly negative degree-six rational polynomial. This gives the single remaining quadratic maximum on the horizontal line, at \(z_+=g+i/2\), where the complex phase is stationary and \(\Re G=-A\). The same clearing on the verticals at \(\Re z=\pm12\sqrt3\) gives negative degree-eight polynomials with no root in \([0,1/2]\).

The unweighted packet real tilt is \(q=1/4\). In round 27 replace the frequency cutoff \(1/2\) by \(b=199/200\). The normalized tilted Gaussian outside \(|\tau|\le b\lambda\) has exponential action at most
\[
 bq-b^2/(2\delta)=-597/1580<-A.
\]
For those retained frequencies, \(s\ge(1-b)\xi\). On the outer real input region \(|z|\ge R\), \(R=12\sqrt3\), reference and actual variances are at least \(c|u|^4/2\), giving the action bound
\[
 -(1-b)^2 R^4/8+q^2\delta/2=-44681/80000<-A .
\]
This leaves the same finite real-input packet, exact Gaussian integration in \(\tau\), and analytic contour argument as in round 27. Compactness transfers each strict action bound to the original analytic reference \(z^4H(az)/4\) as \(a\downarrow0\), and the stationary point continues with its nondegenerate quadratic behavior. The exact graph still forces the exponentially smaller core packet. The real and analytic budgets above therefore certify the stated coefficient.

The file exact_n4.py checks the signs, root counts, saddle factor, nonvanishing denominator, and both outer action inequalities with rational algebra. The parallel floating scan is a separate finite-size diagnostic.
