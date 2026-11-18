# Reduction to the trace of the accumulated heat remainder on its own graph

Let \(\gamma=N(0,1)\), with density \(g\), and let \(f:\mathbb R\to[0,H]\) be continuous with
\[
\int e^{zx-z^2 f(x)/2}\,\gamma(dx)=1\qquad(z\in\mathbb R).
\]
Write \(p_t\) for the centered Gaussian density of variance \(t\), \(P_t\) for its semigroup, and
\(Q_t=P_t\delta_0-\int_{f(x)\le t}P_{t-f(x)}\delta_x\,\gamma(dx)\), with \(P_0\) the identity. Thus \(Q_t=0\) for \(t\ge H\) and \(E_\gamma f=1\). Define the accumulated remainder and its graph trace by
\[
W(t,y)\,dy=\int_t^H Q_s(dy)\,ds,\qquad w(y)=W(f(y),y).
\]
We prove that \(w\ge0\) on all of \(\{f>1\}\) implies \(f\equiv1\). If \(f\) is nonconstant, we locate a negative trace in each connected component of the plateau set \(\{E>f\}\) that meets \(\{f<1\}\), where \(E\) is the local extinction function of [round 169](169-local-extinction.md). Both its depth and a lower bound for its mass under a reverse hitting distribution follow from any fixed early interval in the component.

## Regularity and reverse hitting

Set \(c_0=\sqrt{2/\pi}\). By Tonelli the integrated kernel formula is
\[
W(t,y)=\int_t^H p_s(y)\,ds-
   \int\!\int_{(t-f(x))_+}^{H-f(x)}p_s(y-x)\,ds\,\gamma(dx).
\tag{1}
\]
Zero-age spatial sources make no contribution to the time integral. Each of the two nonnegative terms is at most \(M=c_0\sqrt H\), so \(|W|\le M\). For any fixed \(0\le u\le v\) the function \(\int_u^v p_s(y)\,ds\) is 1-Lipschitz in \(y\), using
\(\int_0^\infty |\partial_y p_s(y)|\,ds=1\) at \(y\ne0\), and continuity at zero. Also the integral of \(p_s(y)\) over any interval of length \(h\) is at most \(c_0\sqrt h\) uniformly in \(y\). Consequently
\[
 |W(t,y)-W(s,x)|\le 2|y-x|+2c_0\sqrt{|t-s|}.
\tag{2}
\]

The occupation identity \(\int_0^H Q_s\,ds=\int_0^1 p_s\,ds\,dy\) follows by integrating the heat equation for \(Q\). The second spatial derivatives of the two integrated measures both equal \(2(\gamma-\delta_0)\), and integrability removes their affine difference. Subtract the contribution up to \(t\) to get
\[
 W(t,y)=\int_0^1p_s(y)\,ds-\int_0^t p_s(y)\,ds+
      \int\!\int_0^{(t-f(x))_+}p_s(y-x)\,ds\,\gamma(dx).
\tag{3}
\]
It follows that \(W(0,y)>0\), \(W(t,y)\ge0\) for \(t\le1\), and \(w(y)\ge0\) wherever \(f(y)\le1\). For \(t>0\), the distributional heat equation for (3) is
\[
 (\partial_t-\tfrac12\partial_{yy}) W(t,y)
        =-g(y)1_{\{f(y)>t\}} .
\tag{4}
\]
For example, the three terms of (3) give, respectively, \(\delta_0-g\), \(-\delta_0\), and \(g1_{\{f\le t\}}\). In the open set \(D=\{(t,y):t>0,\ t>f(y)\}\), \(W\) is caloric by interior heat regularity.

For a point \((t,y)\in D\), let \(Y_u=y+B_u\) and run time backwards. Define
\[
\zeta=\inf\{u\in[0,t]:t-u\le f(Y_u)\}.
\]
The set is nonempty as \(f\ge0\). Path continuity and continuity of \(f\) give \(t-\zeta=f(Y_\zeta)\), including a hit at time zero. For the latter, \(f(Y_t)>0\) would have forced a hit earlier. By Itô's formula on relatively compact subdomains of \(D\), followed by bounded convergence using (1)–(2), \(W(t-u,Y_u)\) stopped at \(\zeta\) is a martingale. Thus, writing \(K_{t,y}\) for the law of \(Y_\zeta\),
\[
 W(t,y)=\int w(x)K_{t,y}(dx),\qquad
 W^-(t,y)\le\int w^-(x)K_{t,y}(dx).
\tag{5}
\]
In particular, a nonnegative graph trace makes \(W\) nonnegative throughout \(D\). The openness in this representation is where continuity of the single-valued graph is used.

## Forced late negative trace and its size

Suppose \(f\) is nonconstant. As \(E_\gamma f=1\), there is a nondegenerate closed interval \(J=[a,b]\) strictly inside \(\{f<1\}\). Put \(\ell=b-a\) and select \(0<\tau<1\) with \(\sup_J f<\tau\). One may take \(\sigma=H\) in the argument below. More precisely, by [round 169](169-local-extinction.md), \(E>1\) at every finite point and \(E\) is constant on each connected component of \(\{E>f\}\). Denote by \(I_*\) the component containing \(J\), and put \(\sigma=E|_{I_*}>1\). Formula (1), the definition of the local extinction function, and continuity give \(W(\sigma,x)=0\) throughout \(I_*\).

With \(\phi(x)=\sin(\pi(x-a)/\ell)\) on \(J\), \(\lambda=\pi^2/(2\ell^2)\), and \(\kappa=\pi/(2\ell)\), set \(N(t)=\int_J W(t,x)\phi(x)\,dx\). Since \([\tau,\sigma]\times J\) is inside \(D\), integration by parts gives
\(N'+\lambda N=\kappa[W(t,a)+W(t,b)]\).
By (3),
\[
 N(\tau)\ge A:=\int_\tau^1\int_J p_s(x)\phi(x)\,dx\,ds>0,
 \qquad N(\sigma)=0.
\]
Integrating the ODE and using \(W\ge0\) on the end traces through time 1 as in [round 161](161-late-endpoint.md) yields
\[
 \int_1^\sigma e^{\lambda(t-\tau)}
       [W^-(t,a)+W^-(t,b)]\,dt\ge A/\kappa .
\tag{6}
\]
This forces strict negativity of \(W\) at some \((t,a)\) or \((t,b)\) with \(1<t<\sigma\). Such a point lies in \(D\), so (5) forces a strictly negative graph trace. This proves the stated reduction using only negative accumulated occupation, directly supplied by round 161.

For the quantitative form put
\[
 B=\int_1^\sigma e^{\lambda(t-\tau)}\,dt,\qquad
 \nu=\frac{1}{2B}\int_1^\sigma e^{\lambda(t-\tau)}
                   (K_{t,a}+K_{t,b})\,dt,\qquad
 m_*=\frac{A}{2\kappa B}.
\]
Here \(\nu\) is a probability measure on spatial locations of reverse graph hits. Equations (5) and (6) imply \(\int w^-\,d\nu\ge m_*>0\). In particular \(m_*\le M\) for an exact graph, and for \(0<d<m_*\),
\[
 \sup_x w^-(x)\ge m_* ,\qquad
 \nu\{w^-\ge d\}\ge \frac{m_*-d}{M-d}.
\tag{7}
\]
The latter uses \(w^-\le d+(M-d)1_{\{w^-\ge d\}}\). All points counted in (7) have \(f>1\), by (3). They also have \(f<E\), since \(W(E(x),x)=0\) by extinction and continuity. With the choice \(\sigma=E|_{I_*}\), they belong to \(I_*\): each finite endpoint of \(I_*\) has \(f=\sigma\) by round 169, and the backward path starting at either endpoint of \(J\) at time \(t<\sigma\) cannot reach a spatial endpoint of \(I_*\) before hitting the graph. Every counted hit has \(f(x)<\sigma\). These conclusions allow either or both spatial endpoints of \(I_*\) to be infinite.

There is also a location bound. A hit with \(f(x)>1\) happens in backward time \(u<t-1\le\sigma-1\). By the reflection principle,
\[
 \nu\{w^->0,\ \operatorname{dist}(x,J)\ge R\}
      \le 4\overline\Phi(R/\sqrt{\sigma-1}),
\tag{8}
\]
where \(\overline\Phi\) is the standard normal upper tail. Choose \(R\) so that (8) is less than \((m_*-d)/(M-d)\). Then there exists a graph point \(x_*\in[a-R,b+R]\) with \(w(x_*)\le-d\). Let \(\omega\) be any modulus for \(f\) on \([a-R-1,b+R+1]\), and take \(0<r\le1\) with \(2r+2c_0\sqrt{\omega(r)}\le d/2\). Equation (2) gives \(w\le-d/2\) and \(f>1\) throughout \([x_*-r,x_*+r]\). In particular the Gaussian mass of this set is at least
\[
 2r\,g\big(\max(|a-R|,|b+R|)+r\big).
\tag{9}
\]
This converts reverse hitting mass to a positive ordinary mass using the continuity modulus of the graph.

## Scope and available controls

For \(s=f(y)>1\), the remaining global target is exactly the inequality
\[
 \int\!\int_0^{(s-f(x))_+}p_u(y-x)\,du\,\gamma(dx)
      \ \ge\ \int_1^s p_u(y)\,du.
\tag{10}
\]
Equations (1)–(9) do not prove (10). They reduce global rigidity to this single trace inequality, or to a trace estimate in one plateau component that contains an early interval. Higher fixed-order resolvent positivity is unnecessary for this reduction.

The exact measurable controls in [round 154](154-deterministic.md) and [round 157](157-deterministic-orders.md) have nonnegative graph trace despite a strictly negative early heat slice. At and after time 1 their \(Q\) is the true survivor of the residual Root transport, so \(W\ge0\) there. On the fat Cantor set their stop time is below 1 and (3) gives nonnegativity. On its complement the post-1 graph trace is zero by Tanaka's formula, because no local time is accumulated at a point after that point enters the stopping set. They do not furnish a continuous-graph counterexample to the reverse-hitting argument.

The integral of the spatial derivative used in (2) was checked exactly with local SymPy 1.14.0. The reduction and the estimates use analytic identities.
