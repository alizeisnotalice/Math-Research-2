# Logarithmic quadratic contact and the second-order exterior jet

Let \(g(u)=p_1(a+u)\), \(L=\log(1/r)\), and \(K_s(z)=\int_0^s p_t(z)\,dt\). For an exact Gaussian output the positive potential \(S(y)=\int g(u)K_{v(u)}(y-u)du\) has an analytic representative, since \(S''=2(P-g)\). Suppose for small \(0<|u|<r_0\),
\[
0\le v(u)\le C\,u^2(\log(1/|u|))^{-k},\qquad k>0 .
\]
Then there is \(\eta>0\) such that
\[
\int_{-r}^r v(u)\,du
=O\!\left(r^3\exp[-\eta(\log(1/r))^{\min(k,1)}]\right).
\]
In particular a profile with a positive multiple of the nominal logarithmic mass cannot give exact Gaussian output. This includes spikes of near-nominal height whose total width in a comparable source annulus is only a negative power of \(\log(1/r)\) times \(r\).

For \(u\ne0\), the first two derivative identities are \(K_s'(z)=-\operatorname{sgn}(z)\operatorname{erfc}(|z|/\sqrt{2s})\) and \(K_s''(z)=2p_s(z)\). These and the formula for \(K_s\) itself make the three off-diagonal coefficients through degree two absolutely integrable. At degree two the near-zero integral is bounded by a constant multiple of \(\int^\infty t^{k/2}e^{-\delta t^k}dt\), with \(t=\log(1/|u|)\) and \(\delta>0\). Uniformly for \(|y|\le r\) and \(|u|\ge 2r\), the third derivative obeys
\[
|K^{(3)}_{v(u)}(y-u)|\le C_1 |u|^{-2}
      (\log(1/|u|))^{3k/2}\exp[-\delta_1(\log(1/|u|))^k]
\]
for sufficiently small \(|u|\). On the remaining exterior bounded variance and fixed separation suffice. Taylor's formula then bounds the integrated degree-two exterior remainder by
\[
Cr^3\left(1+\int^{L+O(1)} e^t t^{3k/2} e^{-\delta_1t^k}\,dt\right)
\le Cr^2\exp[-\eta L^{\min(k,1)}].
\]
In the last step split the integral at \(L/2\) and absorb powers of \(L\) into the exponential. The contributions from \(|u|\le2r\) missing in the first three exterior coefficients obey the same \(O(r^2e^{-\eta L^{\min(k,1)}})\) bound after multiplication by the associated powers of \(y\). This follows from the same Gaussian estimates, integrating in the log variable through infinity for degree two.

Subtract the exterior Taylor polynomial from the analytic \(S\), with \(|y|\le r\), to obtain a quadratic polynomial with fixed coefficients plus the displayed uniformly small remainder for the inner potential. After \(y=rz\), integrate against fixed smooth tests. The size is \(O(r^2L^{-k})\) for each test because \(\int K_v=v\) and the source has the stated variance envelope. Successively divide by \(1,r,r^2\), using tests with nonzero corresponding moments; all three fixed polynomial coefficients must vanish. On the other hand, for \(|u|\le r/2\), a Gaussian centered at \(u\) at any variance at most \(v(u)\) puts at least a fixed positive probability in \([-r,r]\), since \(v(u)/r^2\to0\) uniformly there. Thus the area upper bound for the inner potential implies the stated variance area bound, after rescaling \(r\) by two.

For \(k<1\), one can take any \(\eta<1/(2C)\). For \(k=1\), one can take any \(\eta<\min\{1,1/(2C)\}\). Replace the exterior cut \(2r\) by \(Qr\) with a large fixed \(Q\). On this exterior the Gaussian derivative can be evaluated at the variance envelope, since the normalized displacement is eventually large and the derivative magnitude then increases with variance. Its exponent coefficient is \((1-1/Q)^2/(2C)\). For \(k<1\), split the log integral at \((1-\varepsilon)L\) with small fixed \(\varepsilon>0\). For \(k=1\), integrate \(\exp[(1-\delta_1)t]\) against the polynomial factor directly. For \(k>1\), the exterior derivative is absolutely integrable and the second-order analysis gives the area bound \(O(r^4)\).

## Numerical quantities

At \(r=e^{-L}\), the outer variable is \(w=\log(|u|/r)\), the output shift is \(\eta_0 r\), and the remainder is divided by \(g(0)r^2\). The integral remainder \(K'''_{v}=2p_v'\) is evaluated with a second integral over the line segment from zero to the output shift. The positive and negative source sides have opposite-signed Taylor remainders. The result reports the sum of their magnitudes, retaining a bound even if the signed terms cancel.

For the inner source \(|u|\le2r\), compare the integrated smoothed potential in the output interval \(r[0.6,1.4]\) with the source mass \(\int_{r[0.6,1.4]}g(u)v(u)\,du\), both divided by \(g(0)r^3\). For \(\alpha\ge0\),
\[
\int_0^1\Phi(-\alpha/\sqrt{s})\,ds
=(1+\alpha^2)\Phi(-\alpha)-\alpha\varphi(\alpha),
\]
as follows by integration by parts. This formula and reflection supply the exact time-averaged window probability in the scaled area calculation.

In the dense profile with \(q=cL^{-k}\), the inner sharp mass is asymptotic to \(q(1.4^3-0.6^3)/3\). The signed area correction is asymptotic to \(q^2(1.4^3-0.6^3)\). For example, the diagonal inverse map with uniform heat-time parameter has first correction \(\tfrac14\partial_u^2(gv^2)\), since the mean relative variance is one half. In the scaled model the variance and all fixed derivatives converge after division by \(q\) to those of \(u^2\) on neighborhoods of the window endpoints. Gaussian tails outside such neighborhoods give higher-order errors. The general envelope bound on the absolute window error is \(O(L^{-3k/2})\), by restricting crossing to endpoint neighborhoods of Gaussian width \(O(L^{-k/2})\). For thin pulses, relative errors can remain significant because their sharp source masses are smaller.

The sparse numerical profiles are fixed finite collections of pulses on the log source coordinate, of five scales within the tested core. They model both window-edge leakage and near-off-axis errors. They are reference profiles for the computation, and the exact-output theorem concerns a variance satisfying its own full equation.
