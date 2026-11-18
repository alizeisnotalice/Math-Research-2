# Analyticity of the canonical boundary and rigidity from the high-wing germ

Use the symmetric Gaussian first-exit embedding of round 94, with boundary radius \(R(t)>0\) on positive finite time slabs and inverse graph \(F_*\). The probability of an exit on one wing has time measure \(d\lambda(s)=g(R(s))\,dR(s)\), with \(g=p_1\). The fixed-wall reflection bound and the flux argument of round 88 apply on every such slab: \(R\) is locally Lipschitz, then \(C^1\), and \(R'\) is strictly positive. In particular
\[
 R'(t)=\frac{-p_t'(R(t))+
   \int_{s<t}\big[p_{t-s}'(R(t)-R(s))+p_{t-s}'(R(t)+R(s))\big]d\lambda(s)}
     {g(R(t))}.                                               \tag{1}
\]
We first prove local analyticity of this canonical solution. The subsequent low-wing argument does not require its actual inverse to be Lipschitz.

## 1. A right-sector holomorphic continuation from (1)

Fix \(t_0>0\) and put \(b_0=R(t_0)>0\). Use complex increments \(h\) in \(0<|h|<r,\ |\arg h|<\theta\) for fixed \(0<\theta<\pi/2\), choosing the principal square root for Gaussian kernels at complex ages. Given a bounded holomorphic \(U\) on this sector, let
\[
 X(h)=b_0+h\int_0^1U(vh)\,dv.
\]
This is holomorphic and \(X'=U\), by integration along radii or by holomorphic primitives. In the right side of (1), substitute \(X(h)\) for \(R(t_0+h)\), keep the history \(s\le t_0\) fixed, and write the later integration at \(s=t_0+vh\) as \(h\,dv\) with radius \(X(vh)\) and exit-density factor \(g(X(vh))U(vh)\). This defines a map \(\mathcal T U(h)\).

Here are uniform bounds for the map and its difference on \(\|U\|_\infty\le M\). For the recent same-wing integral, \(\tau=h(1-v)\) and \(d=X(h)-X(vh)\) satisfy \(|d|\le M|\tau|\), giving
\[
 |p'_\tau(d)|\le C_M|\tau|^{-1/2},\qquad
 |p''_\tau(d)|\le C_M|\tau|^{-3/2}.
\]
For two candidate inputs the change in \(d\) is at most \(|\tau|\|U-\widetilde U\|_\infty\). Allowing also changes in \(g(X)\) and \(U\), the recent same-wing term and its Lipschitz constant thus cost \(O_M(\sqrt r)\). On the opposite wing the displacement is near \(2b_0>0\); by shrinking \(r\), the sector gives \(\Re(d^2/\tau)\ge c/|\tau|\). This dominates any fixed derivative power of \(|\tau|^{-1}\), giving \(O_M(r)\) cost for the recent opposite-wing term.

For the same-wing historical integral, set \(\tau=t_0-s+h\), so \(|\tau|\asymp t_0-s+|h|\). On a fixed small historical band, \(d_0=b_0-R(s)=O(t_0-s)\), and \(d\lambda(s)\le Cds\). The reference derivative \(|p_\tau'(d_0)|\) is bounded by \(C(t_0-s)|\tau|^{-3/2}\), yielding an \(O(1)\) historical bound independent of \(M\). The replacement by \(X(h)-R(s)\), or the difference between two such replacements, is bounded after integration by
\[
 C_M |h|\int_0^\eta (u+|h|)^{-3/2}du
          =O_M(\sqrt{|h|}),
\]
times \(1\) or \(\|U-\widetilde U\|_\infty\), respectively. This uses the kernel's second derivative bound along the short complex displacement segment. The older history has a fixed time gap and finite mass, so both its bounds and holomorphic dependence follow directly. The opposite historical wing is spatially separated from zero on the recent band and admits the same complex Gaussian off-diagonal bounds; the older band again has a time gap. Finally the free term is regular, and \(g(X(h))\) is bounded away from zero for small \(r\).

All the integrals define holomorphic functions of \(h\): their integrands are holomorphic almost everywhere in the parameter of integration and the indicated majorants hold on compact sub-sectors, including the integrable \((1-v)^{-1/2}\) at the recent endpoint. Choose \(M\) larger than twice the fixed historical and free bounds and larger than the known real speed bound near \(t_0\). Then choose \(r\) small so \(\mathcal T\) maps the closed ball to itself and has Lipschitz constant less than one. Banach's theorem gives a holomorphic fixed point. On the real interval, the same bounds give uniqueness on the ball of bounded functions, and (1) makes the real canonical speed its fixed point. The resulting \(X\) thus equals \(R(t_0+h)\) on \(0<h<r\).

The choices can be made uniform when \(t_0\) varies over any sufficiently small compact positive slab: the local speed bound and historical recent densities are uniform, the older mass is finite, and \(b_0\) stays bounded away from zero. Taking \(t_0\) slightly to the left of a given positive target time shows real analyticity at that target. Since the real speed is positive, \(F_*\) is real analytic at corresponding positive radii.

## 2. Analytic extension of the canonical interior density across its left edge

Fix the critical positive time \(T\) and write \(b=R(T)\). In the coordinate \(z=x+R(t)\) the canonical interior heat equation is
\[
 u_{zz}=2u_t+2R'(t)u_z .
\]
Its boundary data on \(z=0\), from the flux relation, are
\[
 u(t,0)=0,\qquad u_z(t,0+)=2g(R(t))R'(t).
\]
Both are analytic in \(t\) near \(T\). For clarity, a local analytic solution across the edge can be built by setting \(u_0=0,\ u_1=2g(R)R'\), and
\[
 u_{n+2}=2u_n'+2R' u_{n+1},\qquad
 \widetilde u(t,z)=\sum_{n\ge0}u_n(t)z^n/n!.
\]
To verify convergence, on a complex time disc \(D(T,\rho)\) where \(R'\) and \(u_1\) are bounded use \(\|a\|_n=\sup|a(t)|(\rho-|t-T|)^n\). Cauchy's bound with radius \((\rho-|t-T|)/(n+1)\) gives \(\|a'\|_{n+1}\le e(n+1)\|a\|_n\) for \(n\ge1\); for \(n=0\) use half the remaining radius. Thus
\[
 \|u_{n+2}\|_{n+2}\le C\rho(n+1)\|u_n\|_n+C\rho\|u_{n+1}\|_{n+1}.
\]
Induction gives \(\|u_n\|_n\le C_0 K^n n!\) for fixed large \(K\). On the half-radius time disc the series converges normally for \(|z|<\rho/(2K)\), and solves the stated equation. This is the local analytic Cauchy construction.

In physical coordinates put \(\widetilde q(t,x)=\widetilde u(t,x+R(t))\). Restrict it to \(x>-R(t)\) and extend by zero to the left. On a rectangle near \((T,-b)\) away from the right edge, this truncated function satisfies
\[
 (\partial_t-\tfrac12\partial_x^2)
        [\widetilde q\,{\bf1}_{x>-R(t)}]
 =-g(R(t))R'(t)\,\delta_{x=-R(t)} ,
\]
because its value at the interface is zero and its normal spatial slope is \(2g(R)R'\). This is exactly the source in the canonical survivor equation on that rectangle. Their difference is therefore distributionally caloric there and is zero on the open exterior. A local caloric distribution is spatially real analytic at each interior time; one can obtain this directly by multiplying by spatial and temporal cutoffs and representing with the Gaussian heat kernel. The time-cutoff forcing is separated from the evaluation time, and the space-cutoff forcing is spatially separated, so the Gaussian and its derivatives have uniform local complex spatial bounds. This also recovers a smooth representative. For every fixed nearby time, the zero exterior interval thus forces the difference to vanish throughout the connected spatial neighborhood. Consequently \(\widetilde q\) is indeed the analytic extension of the canonical interior, with uniform spatial derivative bounds on a smaller rectangle. In particular
\[
 |\widetilde q(t,x)|\le C|x+R(t)|
\]
there.

## 3. Equality of the high wing excludes the graph with no low-wing inverse regularity

Assume the continuous single-well hypotheses and unequal limits \(H_+>H_-=1+T>1\) of round 94. The canonical and actual early data through \(T\) agree. In particular their interval at \(T\) is \([-b,b]\) with no critical plateau: a positive-length \(T\)-level plateau outside the entry endpoint \(b\) would contradict the strict-past terminal survivor's exponential transform bound in round 91. Suppose the actual high wing coincides with \(F_*\) in a right neighborhood of \(b\). Then for all \(T<t<T+h_0\), after choosing \(h_0>0\) small, both right endpoints equal \(R(t)\). Denote the actual left radius by \(A(t)=\sup\{r:F_-(r)\le t\}\); the canonical one is \(R(t)\). Monotonicity and continuity of the physical graph give \(A(t)\to b\) as \(t\downarrow T\), even if there are later plateaus. No inverse Lipschitz or Hölder bound on \(A\) is imposed.

Let \(D(t,x)=q_{\rm act}(t,x)-q_{\rm can}(t,x)\), \(\Delta I(t,x)={\bf1}_{[-A(t),R(t)]}(x)-{\bf1}_{[-R(t),R(t)]}(x)\), and \(W(t,x)=\int_T^t D(s,x)\,ds\). All common sources cancel; \(D\) is the canonical minus actual recent **left** heat potential. At all times except the countable set of possible low-wing plateau levels, the only possible zero-age low source in either graph occurs at its current left endpoint. Uniform complex Gaussian bounds off that endpoint, and positive age locally at all other means, give real analyticity of \(D(t,\cdot)\) on \(x>-\min(A(t),R(t))\). The one-sided high-tail support in round 90 applies for \(t<H_+-1\), so both remainders are zero to the right of \(R(t)\). Hence by spatial analytic continuation \(D=0\) on the entire indicated half-line. This conclusion is used only for almost every \(t\).

The time primitive obeys distributionally, at almost every \(t\),
\[
 W_{xx}=2(D+g\Delta I),\qquad W_t=D.                           \tag{2}
\]
For example integrate the heat equation for the difference and the two activated graph measures, starting from their common time-\(T\) data. There is no need to differentiate either inverse. Also \(D\in L^2((T,T+h_0)\times\mathbb R)\) without any low regularity: by Minkowski, the spacetime \(L^2\) norm of a heat potential from any finite recent source measure is at most \(C h_0^{1/4}\) times its mass, since \(\int_0^{h_0}\|p_s\|_2^2ds=O(\sqrt{h_0})\). A zero-age layer at one of the countably many excluded times is immaterial. Hence \(W\) is continuous in time in \(L^2\), with \(L^2\) time derivative. Equation (2) yields \(W(t)\in H^2(\mathbb R)\) for almost every \(t\). On \(x>-\min(A,R)\), both \(D\) and \(\Delta I\) vanish, so \(W\) is affine there, and its \(L^2\) bound forces it to be zero there with zero spatial derivative. In one dimension \(H^2\) gives continuous \(W,W_x\), including at the right edge of the strip between the left endpoints.

On the strip, we now obtain the needed sign in both directions.

- If \(A(t)<R(t)\), then \(\Delta I=-1\). As \(D=U_{\rm can,new}-U_{\rm act,new}\), the upper bound \(D\le U_{\rm can,new}=O(\sqrt{t-T})\) follows from nonnegativity of the actual potential and the canonical bounded recent time density. Thus \(D-g<0\) uniformly in this strip for short times.
- If \(A(t)>R(t)\), then \(\Delta I=+1\). The actual remainder at this fixed time is real analytic on \(x>-A(t)\) in a fixed spatial neighborhood of \(-b\): its left endpoint is the only possible adjacent zero-age source, and the high sources are separated to the right. On \(x>-R(t)\) it equals the canonical interior by \(D=0\). Hence it agrees with the analytic continuation in Section 2 throughout the strip. There \(q_{\rm can}=0\) and
  \[
   |D(t,x)|=|\widetilde q(t,x)|\le C(A(t)-R(t)).
  \]
  The width tends to zero, so \(D+g>0\) for sufficiently short times.

It follows from (2) and the zero right Cauchy data that \(W\Delta I\ge0\) everywhere for almost every short time: integrate \(W_{xx}\) twice to the left across the strip. If the strip is empty the claim is immediate. The \(L^2\) time energy identity with spatial integration by parts now gives
\[
 \frac12\frac{d}{dt}\|W(t)\|_2^2
   =-\tfrac12\|W_x(t)\|_2^2-\int g W\Delta I\,dx\le0.
\]
Since \(W(T)=0\), \(W=0\), then \(D=0\) and \(\Delta I=0\) almost everywhere. In particular at some \(t>T\) the actual signed remainder is compactly supported.

This contradicts the exact future transform at that time. For positive \(z\) applied to the left wing via parameter \(-z\), restrict to the Gaussian window \(x\in[-z,-z+1]\). There \(f(x)\to H_-=1+T\), so with any \(0<\varepsilon<t-T\) the future integral is at least
\[
 c\exp\!\left[\frac{t-T-\varepsilon}{2}z^2\right]
\]
for large \(z\). This follows by completing the Gaussian square, and the interval belongs to \(\{f>t\}\) if the short time \(t\) is chosen below \(1+T\). Compact support, even for a signed finite measure, bounds its Laplace transform by \(C e^{Cz}\), a contradiction.

Thus the higher wing cannot agree with the canonical continuation on a right neighborhood of \(b\), regardless of low inverse regularity. In particular the actual high wing cannot be analytic through \(b\): the canonical inverse is analytic there and both agree on the early side. The same holds in any quasianalytic germ class containing the analytic germs and stable under subtraction. The argument uses actual single-well geometry for identifying the two activation endpoints, and does not assert rigidity of arbitrary disconnected postcritical sublevel sets.
