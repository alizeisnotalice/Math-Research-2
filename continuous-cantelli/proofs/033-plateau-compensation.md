# Persistent high intervals force deep, narrow compensators

Let \(f:\mathbb R\to[0,M]\) be continuous and satisfy
\[
F(t)=\int e^{tx-t^2f(x)/2}\gamma_0(dx)=1\qquad(t\in\mathbb R),
\]
where \(\gamma_t=N(t,1)\). Write \(\nu_t(dx)=e^{tx-t^2f(x)/2}\gamma_0(dx)=e^{-t^2(f(x)-1)/2}\gamma_t(dx)\), and let \(W_t(x)=x-tf(x)\). The normalized identity at \(t+s\) reads \(\mathbb E_{\nu_t}\exp(sW_t-s^2f/2)=1\), so \(\mathbb E_{\nu_t}\exp(sW_t)\le e^{Ms^2/2}\). Equivalently, for independent standard Gaussian \(Y\), \(W_t(X)+\sqrt{M-f(X)}Y\) is \(N(0,M)\) when \(X\) has law \(\nu_t\).

On the positive end assume \(f(x)\ge1-o(1)\) as \(x\to+\infty\), and that \(\nu_t([-R,R])\to0\) for each finite \(R\) as \(t\to+\infty\). The first follows from uniform continuity and the fixed sublevel width estimate below. The second holds if \(f>0\) everywhere or if its only zero is at 0 and \(f(x)\ge c x^2\) there. In the quadratic case the near-zero contribution is at most \(p_1(0)\int e^{tx-ct^2x^2/2}dx=O(t^{-1})\). Compact sets separated from the zero have positive minimum.

For fixed \(\delta>0\), \(\nu_t(f>1+\delta)\le e^{-\delta t^2/2}\) by the shifted form, while \(\nu_t(f<1-\delta)\to0\) by the lower tail hypothesis, compact vanishing, and \(\nu_t(X<-R)\le e^{-Rt}\). Therefore \(f(X)\to1\) in \(\nu_t\) probability. For each fixed \(s\), the factor \(1-e^{-s^2(f-1)/2}\) tends to 0 in \(L^2(\nu_t)\) by boundedness of \(f\); using the bound on the exponential at \(2s\) and Cauchy–Schwarz in the normalized identity shows \(\mathbb E_{\nu_t}e^{sW_t}\to e^{s^2/2}\). It follows that \(W_t\Rightarrow N(0,1)\) and \(X/t=f(X)+W_t/t\to1\) in probability. These limits use only the stated \(o(1)\) lower tail control.

## Verification of the sublevel width at a moving threshold

Write \(A_h=\{x:f(x)\le h\}\) for \(0<h<1\). For \(X>0\) and \(1<\ell<2/h-1\), test \(F\) at \(s=X/h\). On \([X,\ell X]\cap A_h\), the Gaussian-weighted exponent is at least
\[
(1/h-1)X^2/2 +(x-X)[X/h-(x+X)/2]\ge(1/h-1)X^2/2.
\]
Thus \(|A_h\cap[X,\ell X]|\le\sqrt{2\pi}e^{-(1/h-1)X^2/2}\), with Lebesgue measure denoted by bars. Summing geometric bands from \(X=R\) gives
\[
|A_h\cap[R,\infty)|\le\frac{\sqrt{2\pi}e^{-B}}{1-e^{-B(\ell^2-1)}},\qquad B=(1/h-1)R^2/2.
\]
For fixed \(h\) this implies exponentially small total widths, sufficient to prove \(\liminf f\ge1\) under uniform continuity. For fixed \(d>0\), \(0<\epsilon<1\), set \(h=1-d/t\), \(R=(1-\epsilon)t\), \(\ell=1+d/t\) for \(t>d\). This is admissible, \(B\ge d(1-\epsilon)^2t/2\), and \(B(\ell^2-1)\to d^2(1-\epsilon)^2>0\). In particular,
\[
\big|\{x\ge(1-\epsilon)t:f(x)\le1-d/t\}\big|\le C_{d,\epsilon}e^{-d(1-\epsilon)^2t/2}.
\]

## Inward compensators for a high interval

Suppose some fixed \(\delta,r>0\) and a sequence \(t_n\to+\infty\) satisfy \(f\ge1+\delta\) throughout \([t_n-r,t_n]\). Write \(t=t_n\). For any \(0<a<r\), the event \(-a\le W_t\le-a/2\) has limiting \(\nu_t\) probability \(c_a=\Phi(-a/2)-\Phi(-a)>0\). Remove \(f>1+a/(4t)\), whose mass is at most \(e^{-at/8}\), and \(X<(1-\epsilon)t\), whose mass tends to zero. Remaining points cannot have \(X\ge t\), because those points satisfy \(W_t=X-t-t(f-1)\ge-a/4\). The plateau \([t-r,t]\) itself has vanishing mass, and in fact is already outside the retained \(f\)-sublevel for large enough \(t\). Thus \(\nu_t\) assigns probability at least \(c_a-o(1)\) to
\[
\{(1-\epsilon)t\le X<t-r,\ f(X)\le1-(r-a)/t\}.
\]
Here the depth bound follows from \(f(X)=(X-W_t)/t\le(t-r+a)/t\). The moving-threshold width estimate shows that this region of positive tilted mass has ordinary Lebesgue width \(O_{r,a,\epsilon}(e^{-(r-a)(1-\epsilon)^2t/2})\). In particular it has exponentially small \(\gamma_t\) mass and gives a fixed positive lower bound on total variation of the two probability laws along the plateau sequence.

Pick a point \(x_t\) in this region. Fix any \(0<d<r-a\) and \(0<c<d(1-\epsilon)^2/2\). If \(f\le1-d/t\) on \([x_t,x_t+e^{-ct}]\), that interval alone would violate the moving-threshold width estimate for all large \(t\). Hence some \(y_t\) satisfies
\[
0<y_t-x_t\le e^{-ct},\qquad f(y_t)-f(x_t)\ge(r-a-d)/t.
\]
This gives the same lower bound for the variation of \(f\) on that short interval, without requiring differentiability.

If \(\omega_f(h)=\sup_{|x-y|\le h}|f(x)-f(y)|\), these increments give \(\limsup_{h\downarrow0}\log(1/h)\omega_f(h)\ge c(r-a-d)\). Letting \(\epsilon,a\downarrow0\), \(c\uparrow d/2\), and choosing \(d=r/2\), yields the threshold \(r^2/8\). For differentiable \(f\), the mean-value theorem gives a derivative of size at least \((r-a-d)t^{-1}e^{ct}\) at some point between \(x_t\) and \(y_t\). This point is at most \(t\) for large \(t\). Letting \(\epsilon,a\downarrow0\), \(d\uparrow r\), and \(c\uparrow d/2\), yields \(\limsup_{x\to\infty}\log^+|f'(x)|/x\ge r/2\). Both assertions apply on the negative end by reflection; the high interval is always taken inward of its outward endpoint.

## Application and extension to conditioned tilted mass

If \(f\) is uniformly continuous and the end has \(\limsup f>1\), the uniform modulus supplies a fixed high interval length \(r>0\). Therefore a global modulus \(\omega_f(h)=o(1/\log(1/h))\) or uniform continuity plus differentiability with \(\log^+|f'(x)|=o(|x|)\) forces the tail limits to be 1 under compact tilted mass vanishing at both ends. If \(f\) is real analytic, its zeros ensure that vanishing by round 6. Since \(\mathbb E_{\gamma_0}f=1\), a nonconstant solution with the two tail limits 1 would attain a global maximum greater than 1, and round 23 excludes it. The result allows arbitrary uniformly continuous functions in the plateau obstruction and uses analyticity only for the last rigidity application.

For a conditional variant, take \(Q_t=\nu_t(\cdot\mid E_t)\), \(E_t=\{|X-t|\le t/2\}\), and suppose \(p_t=\nu_t(E_t)>0\), \(\log(1/p_t)=o(t)\), and \(W_t\Rightarrow N(0,1)\) under \(Q_t\). A bound \(f\ge1-o(1)\) on the positive tail then gives \(X/t\to1\) under \(Q_t\), since the upper-excess bound survives division by \(p_t\). The plateau proof gives the same oscillation and width constraints, and the positive tilted mass is now at least \((c_a-o(1))p_t\). The hypotheses on the conditional \(W_t\) limit hold in the established nonsaturated origin models of rounds 28 and 29 even if the lower deficit hypothesis there is weakened to \((1-f)_+=o(1)\), while retaining the given contact asymptotics and positivity away from the origin. Indeed, localization outside \(E_t\) uses just an eventual fixed small bound on the deficit. In the core, \(W_t=O(\log(t)/t)\), \(f=O(\log(t)/t^2)\) on the dominant region and the preceding controlled compact errors remain \(o(p_t)\). Subtracting those contributions leaves \(\mathbb E_{Q_t}e^{sW_t-s^2f/2}\to1\). This gives uniformly bounded conditional moments of \(e^{2sW_t}\) using \(f\le M\). With \(f\to1\) in \(Q_t\) probability, Cauchy–Schwarz again gives the conditional Gaussian \(W_t\) limit. The stronger pointwise bound on \((1-f)_+\) used in rounds 28 and 29 controls the additional drift of \(X-t\); this application uses only \(W_t\).

The round 29 application requires a previously proved positive equivalent for the residual (strict subcritical coefficient at either level). The double equality case remains open. For vanishing p_t the fixed positive total variation bound compares Q_t to the shifted Gaussian. For the unconditioned law the lower bound from this set is (c_a−o(1))p_t−γ_t(A).
