# Early Stefan edge bounds and pointwise flux

Assume that on an early time band \((0,T)\) the Cantelli remainder is already identified with killed Brownian motion from 0 in the expanding interval
\[
D_t=(-A(t),B(t)),
\]
where \(A,B>0\) are continuous and nondecreasing. Let \(q(t,x)dx=\mathbb P(B_t^{\rm BM}\in dx,\tau>t)\), and let \(\lambda_\pm\) be the first-exit time measures on the two sides. Assume the early graph identifies them as
\[
 d\lambda_+(t)=g(B(t))\,dB(t),\qquad d\lambda_-(t)=g(-A(t))\,dA(t),
\]
with \(g=p_1\). This is the setting of the connected early sublevel result. There is no immediate exit, and for each positive time the two radii are positive. We derive local regularity away from the initial corner.

## Reflection gives a linear upper edge and finite speed

Survival up to \(t\) implies that the Brownian path has stayed below the fixed endpoint \(B(t)\) at every earlier time, since the true endpoint expands. Reflection at that fixed level gives, for \(z>0\),
\[
0\le q(t,B(t)-z)\le p_t(B(t)-z)-p_t(B(t)+z).
\tag{1}
\]
The inequality holds first almost everywhere as a density bound, and later pointwise by continuity. In particular \(q(t,B(t)-z)\le C z\) uniformly for \(t\) in any compact positive time slab, using \(\sup|p_t'|=1/(\sqrt{2\pi e}\,t)\).

A path that exits on the right during \((t,t+h]\) has survived to \(t\) and must cross the fixed level \(B(t)\) during the next \(h\) units of time. Applying the Markov property and reflection again gives
\[
\lambda_+((t,t+h])\le\int_0^\infty q(t,B(t)-z)\,2\bar\Phi(z/\sqrt h)\,dz\le Ch/2,
\tag{2}
\]
since \(\int_0^\infty u\bar\Phi(u)du=1/4\). The lost probability from paths exiting the other side is immaterial for this upper bound. The graph identity for \(\lambda_+\) and the positive minimum of \(g\) on the compact range of \(B\) now make \(B\) locally Lipschitz; the same argument applies to \(A\). In fact retaining the full right side of (1), substituting \(z=\sqrt h\,u\) in (2), and using dominated convergence gives at each \(t>0\)
\[
 g(B(t))D^+B(t)\le \frac{B(t)}{t}p_t(B(t)),
\]
where \(D^+\) is the upper right Dini derivative. This part alone excludes a graph increment \(f(x+h)-f(x)=o(|h|)\) at an early wing: the inverse bound \(|B(t')-B(t)|\le L|t'-t|\) gives the opposite linear lower bound on every local time increment. The analogous statement holds in the left radial coordinate.

## Jump relation and continuous Stefan speed

After local Lipschitz regularity the exit measures have locally bounded time densities \(\rho_+=g(B)B'\), \(\rho_-=g(-A)A'\) almost everywhere. The potential representation is
\[
q(t,x)=p_t(x)-U_+(t,x)-U_-(t,x),
\quad U_+(t,x)=\int_{s<t}p_{t-s}(x-B(s))d\lambda_+(s),\quad
U_-(t,x)=\int_{s<t}p_{t-s}(x+A(s))d\lambda_-(s).
\]
It gives joint continuity on compact positive slabs: local recent ages cost \(O(\int_0^\delta\tau^{-1/2}d\tau)\) uniformly, and the older part is continuous by separation in time. The density is zero laterally and is caloric on the open time-space domain by the same kernel separation from source points near every interior point.

Fix a right-edge time \(t\). The other-wing derivative \(\partial_xU_-\) is continuous through \(B(t)\), since the other sources are spatially separated. In \(\partial_x U_+(t,B(t)+w)\), let \(\tau=t-s\), \(d=B(t)-B(s)\), and split
\[
p_\tau'(w+d)=p_\tau'(w)+[p_\tau'(w+d)-p_\tau'(w)].
\]
For recent times local Lipschitz gives \(|d|\le L\tau\), while \(\|p_\tau''\|_\infty\le C\tau^{-3/2}\). Hence the bracket is bounded by \(C\tau^{-1/2}\) and converges under the integral to \(p_\tau'(d)\). The older part is smooth. The integral of the static term \(p_\tau'(w)\) against \(\rho_+(t-\tau)\) tends to \(-\operatorname{sgn}(w)\rho_+(t)\) at left Lebesgue points of \(\rho_+\). Indeed \(|p_\tau'(w)|\), viewed as a function of \(\tau>0\), has integral 1 and is an approximate identity of scale \(w^2\); its rescaled kernel is bounded by an integrable decreasing envelope \(C(1+\tau)^{-3/2}\). The older static term tends to zero even if the older exit time measure is only finite. Thus use of the exterior derivative \(q_x=0\) gives almost everywhere
\[
\rho_+(t)=-p_t'(B(t))
 +\int_{s<t}p_{t-s}'(B(t)-B(s))d\lambda_+(s)
 +\int_{s<t}p_{t-s}'(B(t)+A(s))d\lambda_-(s).
\tag{3}
\]
The right side is continuous on every compact positive time slab. In the same-wing integral the recent part is uniformly \(O(\sqrt\delta)\) by the \(\tau^{-1/2}\) bound, and the remaining part is smooth in the spatial argument and continuous in time with finite measures. One may use a continuous cutoff in \(t-s\) for this continuity step. The opposite kernel is separated and continuous through age zero. Thus \(\rho_+\) has a continuous representative. Since \(g(B)>0\), absolute continuity of \(B\) yields \(B'(t)=\rho_+(t)/g(B(t))\) at every positive time, and likewise for \(A\). The jump relation now applies at every such \(t\). Its interior limit is the exterior derivative minus \(2\rho_+\), or
\[
-q_x(t,B(t)-)=2g(B(t))B'(t),\qquad q_x(t,-A(t)+)=2g(-A(t))A'(t).
\tag{4}
\]
These are also the linear density coefficients from the inside because \(q\) is zero at the endpoints. Formula (3) recovers the upper speed bound since both integrals use \(p_\tau'(d)\le0\) at \(d\ge0\). If the opposite wing already has positive exit mass, its term is strictly negative.

The same formula gives local \(1/2\)-Hölder slopes. For same-wing time increments of length \(h\), integrate the recent \(O(\tau^{-1/2})\) bound over ages less than \(2h\), and use \(|\partial_\tau p_\tau'(d)|+|\partial_dp_\tau'(d)|\le C\tau^{-3/2}\) when \(|d|\le2L\tau\) and \(\tau\ge h\). Since the spatial increment of \(B\) is \(O(h)\), the old-age difference costs \(O(h\int_h^{\delta}\tau^{-3/2}d\tau)=O(\sqrt h)\). Older separated time pieces are Lipschitz in the change of parameters, as are the other-wing integrals and the free term. Thus \(A,B\in C^{1,1/2}_{\rm loc}\).

## Positive edge slope by a translating strip

The lower counterpart to (1) follows by comparison, without an assumed Hopf lemma on a moving boundary. The killed density is positive at every interior point. To see this, there is some positive survival at a sufficiently small time. Because the domains expand, confining Brownian motion to that fixed early interval preserves a positive amount of probability up to any later finite time. For any interior point \((t,x)\), choose \(s<t\) so close that \(x\in D_s\). The strictly positive Dirichlet kernel of the fixed interval \(D_s\) takes the positive surviving mass to positive density at \(x\) by time \(t\). Positivity of that kernel follows, for example, by conditioning on two interior endpoints and restricting the Brownian bridge to a tube around the straight segment joining them. Continuity makes \(q\) bounded below on each compact subset of the interior.

Fix a compact positive time slab, let \(L\) be a Lipschitz bound for the right endpoint on a slightly larger slab, and let \(k>L\). Choose fixed small \(d,\delta>0\) so that for every target time \(t\) on the slab the translated strip
\[
I_s=(B(t)-k(t-s)-d,\ B(t)-k(t-s)),\qquad t-\delta\le s\le t,
\]
has its left edge inside the true interval and its right edge strictly inside before time \(t\). One can take \(d\) and \(k\delta\) each smaller than a quarter of the positive minimum of \(B\) on the enlarged slab. At the initial time the right margin is at least \((k-L)\delta\), so the union of all initial closed strips is a compact interior set. Let \(z\) be the coordinate to the right of the translating strip's left edge. The explicit function
\[
 e^{-(k^2+\pi^2/d^2)(s-t+\delta)/2}e^{-kz}\sin(\pi z/d)
\]
is caloric in the strip, positive inside and zero laterally. Multiply it by the uniform positive lower bound of \(q\) on the initial sets and compare on each strip up to time \(t\). At the right endpoint this gives \(q(t,B(t)-z)\ge c z\) for small \(z>0\), uniformly on the slab. Reflect for the left endpoint. Therefore (4) has strictly positive flux on both sides, bounded positively below on compact positive time slabs. This also excludes a superlinear density vanishing at the edge. The corresponding inverse graph is \(C^{1,1/2}_{\rm loc}\) on both early wings with slope bounded above and positively below on their compact interior subintervals.

## Executed checks

`edge_check.py` verifies the translating sine solution and the reflection moment \(1/4\) exactly with SymPy. It also calibrates the one-sided reflection rate and the finite-time mass of the jump kernel at 50 and 75 digits with mpmath. The proof above uses the exact identities and analytic estimates, without inferring any regularity at the initial vertex.
