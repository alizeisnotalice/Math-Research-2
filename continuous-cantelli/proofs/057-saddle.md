# Subexponential signs and gains in an asymmetric finite-order Wick group

Write \(d=m-1\), \(\nu=k/N\), and \(a=(cN)^{-1/(2d)}\). Keep \(\nu\) in a compact subset of \((0,1/(2d))\) and \(m\) fixed. For an analytic reference with
\[
 \log H(u)=\alpha u+\beta u^2+O(u^3),
\]
the two real saddle phases of a round-51 diagonal Wick group admit a more precise expansion. Set
\[
 R=\sqrt{(1-2d\nu)/\lambda},\qquad U=R/d,\qquad S=R+U,
 \qquad C=2\lambda R^2,\qquad E=2d\nu/m.
\]
In log-radius coordinates \(x,y\), the unperturbed phase Hessian is
\[
 J=\begin{pmatrix}C+E&-E\\-E&E\end{pmatrix}.
\]
The log amplitude inserted into the phase at sign \(\sigma\in\{-1,1\}\) is
\[
 a\,q_{1,\sigma}(x,y)+a^2q_2(x,y)+O(a^3),\qquad
 q_{1,\sigma}=2\nu\sigma\alpha(e^x+e^y),\quad
 q_2=2\nu\beta(e^x+e^y)^2.
\]
The stationary radius shift to first order in \(a\) is
\(-aJ^{-1}\nabla q_{1,\sigma}\). Inserting it into the phase gives
\[
 \Phi_\sigma
 =\Phi_0+2a\nu\sigma\alpha S
   +a^2\left[2\nu\beta S^2
   -2\nu^2\alpha^2\left(\frac{S^2}{C}+\frac{U^2}{E}\right)\right]+O(a^3).
\]
In particular \(\alpha\) controls which real saddle is favorable, while the
second-order stationary-phase correction is the same for both signs.

The round-51 saddle proof gives a positive sum of the two real saddle contributions on this compact range. The logarithm of the ratio of a diagonal group to \(g(0)\) times its unperturbed algebraic value is therefore
\[
 2Na\nu|\alpha|S
  +Na^2\left[2\nu\beta S^2
   -2\nu^2\alpha^2\left(\frac{S^2}{C}+\frac{U^2}{E}\right)\right]
  +O(Na^3+1).
\]
Indeed, the log of the positive sum differs from the dominant log by a uniformly bounded amount, the Gaussian prefactors are uniformly comparable to the corresponding unperturbed prefactor, and \(g(a(\pm S))=g(0)(1+O(a))\). Here \(g\) is the standard Gaussian source translated by the contact point.

The finite diagnostic uses \(H(u)=\exp(\alpha u+\beta u^2)\). After division by \(g(0)\), the final factors of the group in its single variable \(Z=W+T\) equal \(\exp(AZ+BZ^2)\) with
\[
 A=a(2k\alpha-x_0),\qquad B=a^2(2k\beta-1/2).
\]
For even coefficient order, replacing \(A\) by \(|A|\) preserves the result. With \(B\ge0\) the finite sum is positive, and we use a log recurrence or an explicit multinomial convolution. Local 2D Cauchy integration at positive stationary radii checks separately selected cases. This diagnostic and the displayed expansion do not add a hypothesis-free result for divergent Taylor jets.
