# Sign crossings, radial monotonicity, and early heat support

Assume \(f\) is continuous, nonnegative and bounded on \(\mathbb R\), with
\[
 F(\lambda)=\int e^{\lambda x-\lambda^2 f(x)/2}\,\gamma(dx)=1
\]
for all real parameters. Write \(g=p_1\), \(h=f-1\), and \(H=\sup f\). The original conditional variance is a constant minus \(f\), so its radial direction is reversed. For a nonconstant exact graph,
\[
 \int h\,d\gamma=\int xh\,d\gamma=0,\qquad
 \int x^2h\,d\gamma=\tfrac12\int h^2\,d\gamma>0.               \tag{1}
\]
Expansion through degree four yields (1). Jensen under the Gaussian shift gives
\[
 \int e^{\lambda x}h(x)\,\gamma(dx)>0\quad(\lambda\ne0),        \tag{2}
\]
because the shifted expectation of \(e^{-\lambda^2h/2}\) is 1 and \(h\) is nonconstant.

## 1. Radial centers and given sign patterns

Suppose \(f(x)=F_0(|x-c|)\) with \(F_0\) nonconstant and monotone. If \(c>0\), pairing \(x\) with \(-x\) at \(x>0\) gives \(|x-c|<|x+c|\), so \(\int x f\,d\gamma\) has a strict sign, opposite to the monotonicity direction. Equality would give \(F_0(r)=F_0(r+2c)\) for all \(r\ge0\), forcing constancy. Reflect for \(c<0\). Hence \(c=0\). At this center \(F_0\) cannot be nonincreasing, since the independent-copy covariance inequality gives \(\operatorname{Cov}(X^2,F_0(|X|))\le0\), contrary to (1). The nondecreasing direction is excluded in Sections 2–5 by early heat support and uniqueness of its boundary measure.

A genuine single sign crossing of \(h\) about \(c\) is impossible by multiplying by \(x-c\) in (1). The same argument excludes a single crossing between any two distinct exact \(f\)'s, using the equality of their first two displayed moments. For \(h=f-1\) with signs \((-,+,-)\) at two cuts \(c_1<c_2\), multiply by \((x-c_1)(x-c_2)\). The integral is nonpositive by the signs and is strictly positive by (1). For signs \((+,-,+)\), (2) is automatic once the first two moments vanish: subtract the secant affinely interpolating \(e^{\lambda x}\) at \(c_1,c_2\). The difference has precisely the exterior-positive and interior-negative signs and is strict away from the cuts. Thus (2) by itself cannot exclude this pattern.

More generally, for a finite pattern with a negative rightmost sign block, compare \(e^{\lambda x}\) times all positive mass below the last cut with the strictly negative contribution from a compact interval separated to its right. For large \(\lambda>0\) the latter dominates, contrary to (2). The leftmost nonzero block cannot be negative either, by reflection. This applies to genuine sign blocks, including an eventual zero suffix after the last negative block.

## 2. First-order origin contact in the remaining radial direction

Let \(f\) be even and nondecreasing on \([0,\infty)\). If it is nonconstant then \(H=\lim_{r\to\infty}f(r)>1\). The exact identity forces \(f(0)=0\): otherwise the compact positive contribution decays for \(\lambda\to+\infty\), the negative half is \(O(1/\lambda)\), and the positive tail with \(f\ge1+\delta\) contributes at most \(e^{-\delta\lambda^2/2}\). A positive-length plateau at zero is impossible by testing the real identity on a subinterval away from zero. Thus \(f(r)>0\) for \(r>0\). In fact,
\[
 f(r)\sim \frac{r^2}{2\log(1/r)}\qquad(r\downarrow0).          \tag{3}
\]

For the lower bound, fix \(0<\vartheta<1\). Restrict \(F(\lambda)\) to \([\vartheta r,r]\), use \(f(x)\le f(r)\) there and optimize at \(\lambda=\vartheta r/f(r)\). It follows that
\[
 f(r)\ge\frac{\vartheta^2 r^2}{2\log(1/\gamma([\vartheta r,r]))}.
\]
Send \(r\downarrow0\) and then \(\vartheta\uparrow1\). For the upper bound suppose there is a sequence \(r_n\downarrow0\) with
\(f(r_n)\ge(1+\delta)r_n^2/[2\log(1/r_n)]\), where \(0<\delta\le1\). Put \(\kappa=\delta/16\) and choose \(\lambda_n\to\infty\) with
\(r_n=(2-\kappa)L_n/\lambda_n\), \(L_n=\log\lambda_n\). Write \(x=L_n z/\lambda_n\), so the Jacobian supplies \(L_ne^{-L_n}\). On \(|z-2|\le\kappa/2\), monotonicity propagates the elevated value at \(r_n\) and the exponent including \(e^{-L_n}\) is at most
\[
 L_n\big[z-1-(1+\delta)(2-\kappa)^2/4\big]\le -c_\delta L_n
\]
for large \(n\), since \(\log(1/r_n)\le L_n\). On \(1/2\le z\le8\) outside that interval use the already proved lower envelope
\(f(x)\ge(1-\varepsilon)x^2/[2\log(1/x)]\), for fixed sufficiently small \(\varepsilon\). As \(\log(1/x)\le L_n\), the exponent is
\(L_n[z-1-(1-\varepsilon)z^2/4]\), uniformly negative by the gap of the square \(-(z-2)^2/4\). For \(z\ge8\) in the small physical origin window the same bound is at most \(-cL_n z^2\); for \(z<1/2\) use \(\lambda_n x\le L_n/2\) and nonnegativity of \(f\). The compact positive sources away from zero, negative half and far positive tail all vanish as before. Therefore \(F(\lambda_n)\to0\), a contradiction.

## 3. Signed heat support at early times

Let
\[
 Q_t=P_t\delta_0-\int_{f(x)\le t}P_{t-f(x)}\delta_x\,\gamma(dx),
\]
so \(Q_t=0\) for \(t\ge H\). For fixed \(y>0\), its possible signed-density version has finite time variation: both the free and past terms are integrable in time up to \(H\) by \(\int_0^H p_s(z)\,ds\le C\sqrt H\) uniformly in \(z\). The one-dimensional heat resolvent, and \(F(\sqrt{2q})=1\), give for each \(q>0\)
\[
 \int_0^H e^{-qt}Q_t(y)dt
 =\sqrt{\frac2q}\int_{x>y}
   e^{-qf(x)}\sinh\bigl(\sqrt{2q}(x-y)\bigr)\,\gamma(dx)>0.      \tag{4}
\]
Indeed the free term is \(e^{-\sqrt{2q}y}/\sqrt{2q}\), and the subtracted past resolvent at mean \(x\) is \(e^{-qf(x)-\sqrt{2q}|x-y|}/\sqrt{2q}\). Insert the exact real identity into the free term and cancel the entire part \(x\le y\).

For a radial nondecreasing \(f\), (4) decays exponentially in \(q\) at every rate less than \(\min(f(y),H-1)\). On \(y<x<R\), use \(f(x)\ge f(y)\). On \(x>R\) use \(f\ge H-\varepsilon\) and the unrestricted Gaussian integral \(\int e^{\sqrt{2q}x}\gamma(dx)=e^q\). For reference, a finite signed time measure supported in \([0,H]\) with Laplace transform \(O(e^{-aq})\) as real \(q\to+\infty\) is supported in \([a,H]\). For any \(\delta<a\), multiply by \(e^{\delta z}\). On the imaginary and positive real axes this is bounded, and its entire growth is order at most one. The Phragmén–Lindelöf argument in the two right quadrants gives a bound throughout \(\Re z\ge0\). For a smooth test supported in \(t<\delta\), Fourier inversion for the exponentially tilted translated measure, with the same two integrations by parts as in round 58, forces that test pairing to vanish. Send \(\delta\uparrow a\).

It follows that \(Q_t(y)=0\) for \(t<\min(f(y),H-1)\). First this is true for almost every \(t\) at each fixed \(y\). Away from \(t=f(y)\) the local source is separated from zero age and the remote source is spatially separated; the same bounds give continuity, so the equality holds pointwise on the open set in question. Reflect for \(y<0\). For all \(0<t<H-1\),
\[
 \operatorname{supp}Q_t\subset \{y:f(y)\le t\}=[-a(t),a(t)],
 \qquad a(t)=\sup\{r\ge0:f(r)\le t\}.                        \tag{5}
\]
The support is in the past interval, without any assertion of positivity.

There are no plateaus and every outer edge at \(0<t_0<H-1\) is steep on both sides:
\[
 \frac{t_0-f(a-u)}{u^2}\longrightarrow\infty,\qquad
 \frac{f(a+u)-t_0}{u^2}\longrightarrow\infty,\quad a=a(t_0).    \tag{6}
\]
Here is a purely temporal proof. By (5), the past output at time \(t_0\) matches \(p_{t_0}\) immediately to the right of \(a\). All its means are to the left; monotone convergence from the right gives its positive-age trace \(p_{t_0}(a)\). Fatou at fixed \(y=a\) as \(t\downarrow t_0\) from above shows
\(\limsup Q_t(a)\le0\); more precisely fixed past sources separated from \(a\) can supply arbitrarily close to that trace. If either ratio in (6) were bounded along \(u_n\downarrow0\), monotonicity would give a strip on that side of length \(u_n\) with time differences at most \(C u_n^2\). For every \(\delta\in[K u_n^2,2K u_n^2]\), \(K>2C+1\), all its sources have at time \(t_0+\delta\) Gaussian ages bounded above and below by positive multiples of \(u_n^2\). Their additional contribution at \(a\) is at least a fixed \(c>0\). These strips are disjoint from the fixed separated past sources, so \(Q_{t_0+\delta}(a)\le-c/2\) throughout the stated time window for large \(n\). On the rest of sufficiently small positive \(\delta\) its upper bound can be taken to be any prescribed \(\eta>0\). Since it vanishes below \(t_0=f(a)\), testing its shifted time-Laplace transform with \(q=(K u_n^2)^{-1}\) gives a strictly negative value for large \(n\), after multiplying by \(q\) and choosing \(\eta\) small. The times separated from \(t_0\) are exponentially negligible by finite time variation. This contradicts (4). Any plateau touching the left of \(a\) has \(t_0-f(a-u)=0\) near the endpoint, and is therefore impossible in the early band.

At a fixed early time, stiffness (6) also implies spatial continuity through each edge: near \(a\) the age \(v=t-f(x)\) obeys \(v\gg(a-x)^2\), so \(p_v(0)\le C p_v(a-x)\); this last function of \(x\) is integrable by the known exterior trace. Dominated convergence now gives the interior trace zero for \(Q_t\). This fixed-time argument does not provide joint continuity in \((t,y)\).

## 4. Exclusion with regularity of the inverse boundary

Suppose in addition that \(a(t)\) is locally Hölder with exponent strictly greater than \(1/2\) on \(0<t<H-1\). This includes, for example, \(f\) of class \(C^1\) with strictly positive radial derivative wherever \(0<f(r)<H-1\). Then the centered radially nondecreasing graph cannot be exact.

Indeed, on every compact subslab of the early time band, the past heat term is jointly continuous all the way to the moving edges. The time distribution of the boundary-source measure on both wings assigns an interval of length \(\eta\) mass \(O(\eta^\beta)\) for some \(\beta>1/2\) on the compact subslab. The contribution from ages at most \(\eta\) to the density is bounded uniformly by \(C\eta^{\beta-1/2}\) by dyadic summation against \(p_s(z)\le C s^{-1/2}\). The remaining positive ages are continuous by dominated convergence. Hence \(Q\) is continuous and zero at the lateral boundary of the expanding interval in (5). Inside it is caloric, since graph times have already passed.

Fix a small initial time \(t_i>0\). Add free heat propagation of the negative part of \(Q_{t_i}\) to \(Q\) in the moving interval, making the initial and lateral data nonnegative. The forward maximum principle on each bounded continuous moving domain gives nonnegative values. But \(\|Q_{t_i}^-\|_1\le\gamma\{f\le t_i\}\to0\), since \(f>0\) away from zero. Sending \(t_i\downarrow0\) proves \(Q_t\ge0\) throughout the early band.

This contradicts its moments near \(H-1\). By (5) and positivity,
\[
 \int y^2Q_t(dy)\le a(t)^2\gamma\{f>t\}.
\]
The exact future transform instead makes the second moment
\[
 \int_{f>t}\bigl[x^2-(f(x)-t)\bigr]\gamma(dx).
\]
For \(a>0\), \(E[X^2-a^2\mid |X|>a]=1+a g(a)/\bar\Phi(a)-a^2>1\), using the strict Mills bound. As \(t\uparrow H-1>0\), \(a(t)\) remains in a compact interval bounded away from zero because \(f(0)=0\) and its tail limit is \(H\). Thus the excess over 1 has a positive uniform margin whereas \(f-t\le H-t\to1\). The second moment is strictly too large. Pointwise steepness does not directly give the joint continuity used in this elementary proof. Section 5 bypasses that regularity requirement through finite boundary energy.

## 5. Removal of the inverse-boundary regularity assumption

In fact the centered radially nondecreasing direction is excluded under continuity alone. The key is uniqueness of the early boundary measure using its finite heat energy.

Fix any \(0<T<H-1\). The plateau exclusion after (6) gives a continuous strictly increasing early inverse \(a(t)\) on \(0<t\le T\), with \(a(0)=0\). By (3),
\[
 a(t)^2/t \sim \log(1/t).
\]
Let \(\mu\) be the graph measure \(g(x)dx\,\delta_{f(x)}(dt)\) restricted to \(f(x)\le T\).
For a fresh Brownian motion \(B\) from zero, let \(\tau\) be its first exit from \((-a(s),a(s))\) at positive times up to \(T\), with \(\tau=\infty\) if no exit occurs, and let \(\eta\) be the law of \((\tau,B_\tau)\) restricted to \(\tau\le T\).
The expanding boundary is continuous, so for a positive exit time \(|B_\tau|=a(\tau)\). The reflection bound on each dyadic interval \([2^{-n-1},2^{-n}]\), using the smaller radius \(a(2^{-n-1})\) and the larger Brownian horizon \(2^{-n}\), gives
\[
 {\mathbb P}(\tau\le t)\le C t^{c}\quad (t\downarrow0)
\]
for some \(c>0\) (for example \(c=1/8\) once \(a(s)^2/s\ge\frac12\log(1/s)\)). In particular there is no immediate exit. Neither measure has a time atom: for \(\eta\) at any fixed \(t>0\) use the zero probability of \(|B_t|=a(t)\), and for \(\mu\) use strict radial increase. Also \(\mu([0,t]\times\mathbb R)\le C a(t)\le Ct^{1/3}\) near zero.

For every \(0<t\le T\), the heat potentials from earlier \(\mu\) and earlier \(\eta\) both equal \(p_t(y)\) on \(|y|>a(t)\). For \(\mu\) use (5); for \(\eta\) use strong Markov at \(\tau\), since a Brownian path observed outside must have exited. There is no atom at time \(t\), and these exterior potentials have continuous densities there. Since all prior source means lie in \([-a(t),a(t)]\), let \(y\downarrow a(t)\) or \(y\uparrow-a(t)\) and use monotone convergence of each positive Gaussian contribution. This gives both prior potentials at each boundary point exactly equal to \(p_t(a(t))\).

Write \(\kappa=\mu+\eta\), \(\sigma=\mu-\eta\), and use the kernel
\[
 {\cal K}((s,x),(t,y))=p_{|t-s|}(x-y)\quad(s\ne t).
\]
The diagonal has \(\kappa\otimes\kappa\) mass zero. Positivity and the two traces give
\[
 \iint {\cal K}\,d\kappa\,d\kappa
       =4\int p_t(a(t))\,\kappa(dt,\mathbb R)<\infty .
\]
Indeed \(p_t(a(t))=t^{-o(1)}\), and both measures have a polynomial near-zero time mass bound. Thus the energies for the signed differences are absolutely integrable. The past trace of the difference is zero at both graph points of every positive \(t\), so by time symmetry
\[
 \iint {\cal K}\,d\sigma\,d\sigma=0.                           \tag{7}
\]

We justify strict positive definiteness for this signed energy without assuming any uniform steepness. Choose increasing subsets \(E_j\) of the graph, requiring radius \(|x|\ge1/j\) and
\[
 |f(x')-f(x)|\ge |x'-x|^2
 \quad\text{for every same-wing }x'\text{ with } |x'-x|\le1/j .
\]
They are Borel by testing on a countable dense set and then continuity, and exhaust the positive-time graph by (6) and reflection. Let \(\sigma_j=\sigma|_{E_j}\).
For same-wing pairs in \(E_j\) within distance \(1/j\), writing \(d=|x-y|\), \(\ell=|t-s|>0\), gives \(d^2\le\ell\) and
\[
 p_{\ell+\varepsilon}(d)/p_\ell(d)
   =\sqrt{\ell/(\ell+\varepsilon)}
     \exp\!\left(\frac{\varepsilon d^2}{2\ell(\ell+\varepsilon)}\right)
   \le e^{1/2}.
\]
For remaining pairs \(d\ge1/j\), the mollified kernel is bounded by \(Cj\). Thus fixed-\(j\) spatially mollified energy converges to the bare energy as \(\varepsilon\downarrow0\), dominated by the already finite positive \(\kappa\) energy and total mass. Opposite-wing pairs are also covered by the distance bound. Time diagonals remain null.

At spatial frequency \(\xi\), let \(\widehat{\sigma_{j,\xi}}(dt)=\int e^{i\xi x}\sigma_j(dt,dx)\) and
\[
 E_j(\xi)=\iint e^{-\frac12\xi^2|t-s|}
                  \widehat{\sigma_{j,\xi}}(dt)\,
                   \overline{\widehat{\sigma_{j,\xi}}(ds)} .
\]
This is nonnegative. For \(\xi\ne0\) the temporal kernel has strictly positive Fourier transform \(2\lambda/(\lambda^2+\omega^2)\), \(\lambda=\xi^2/2\); for \(\xi=0\) the energy is the squared mass. Spatial Fourier inversion expresses the mollified energy as
\((2\pi)^{-1}\int e^{-\varepsilon\xi^2/2}E_j(\xi)\,d\xi\).
By monotone convergence as \(\varepsilon\downarrow0\), and the physical-space domination just proved, this yields
\((2\pi)^{-1}\int E_j(\xi)d\xi=\iint{\cal K}\,d\sigma_j\,d\sigma_j\).
Since \(\sigma_j\to\sigma\) in total variation, the frequency energies converge pointwise to \(E(\xi)\) for \(\sigma\). By dominated convergence in physical space and then Fatou in Fourier space, (7) forces \(\int E(\xi)d\xi=0\). The function \(E\) is continuous, hence zero everywhere. Its strict temporal positivity for each \(\xi\ne0\) yields \(\int e^{i\xi x}\sigma(dt,dx)=0\) as a measure in time. Uniqueness of the joint Fourier transform, also using continuity at \(\xi=0\), implies \(\sigma=0\), hence \(\mu=\eta\).

We obtain \(Q_t\ge0\) for all \(0<t\le T\), since it is the actual surviving Brownian mass. As \(T<H-1\) was arbitrary, the compact-support moment contradiction in Section 4 now applies to the continuous case without additional inverse regularity. Together with Section 1 this excludes every nonconstant continuous graph monotone in \(|x-c|\), in either monotonicity direction and for any center.

## 6. Two distinct wings of a continuous single-well graph

Assume more generally that \(f\) has a single valley: it is nonincreasing as \(x\) moves from \(-\infty\) to the minimum interval and nondecreasing from that interval to \(+\infty\). Write \(H_-,H_+\) for its two tail limits. In a nonconstant exact graph both limits exceed 1. For if \(H_-\le1\), the whole left wing has \(h\le0\), with strictly negative \(h\) at an interval to the left of all positive mass, and (2) fails for \(\lambda\to-\infty\). The argument is identical on the right. A positive compact minimum would make the transform vanish at large real parameters. A zero interval is impossible, and if the unique zero were \(c\ne0\), take large real parameters of sign opposite to \(c\); compact mass near the zero then also vanishes, a contradiction. Thus the unique minimum is \(f(0)=0\). Section 2 applies to each wing separately, with the opposite half suppressed by the real parameter. Both wings have the same leading origin contact (3).

Let \(H_*=\min(H_-,H_+)>1\) and \(T_*=H_*-1\). Reflect (4) at negative \(y\), and use the corresponding wing in each direction. The same proof gives the early support \([-A(t),B(t)]\) for \(0<t<T_*\), where \(-A(t),B(t)\) are the endpoints of \(\{f\le t\}\), as well as absence of plateaus and both-sided steepness on each moving edge by its own temporal positivity. The first-exit argument of Section 5 now takes the interval \((-A(t),B(t))\). For the dyadic bound near zero use the smaller of \(A(t),B(t)\); both are asymptotic to \(\sqrt{t\log(1/t)}\). The prior potentials at their respective edges equal \(p_t(-A(t))\) and \(p_t(B(t))\), each at most \(t^{-o(1)}\). The positive-definite argument uses the local steep property on each wing and the same separated opposite-wing bounds. Hence \(Q_t\ge0\) and has the asserted compact support throughout the early band.

Use the quadratic test \(P(y)=(y+A(t))(y-B(t))\), nonpositive on that support. Its exact moment is
\[
 \int P(y)\,Q_t(dy)=\int_{f>t} \bigl(P(x)-(f(x)-t)\bigr)\gamma(dx)\le0.
\]
For \(m(z)=g(z)/\bar\Phi(z)\),
\[
 E[P(X)\mid X>B]=1+A[m(B)-B],\qquad
 E[P(X)\mid X<-A]=1+B[m(A)-A].
\]
These follow by the elementary first two Gaussian tail moments, and each excess over 1 is positive. Let \(A,B\) now stand for the limits of the endpoints as \(t\uparrow T_*\) (the entry endpoints if a plateau occurs at \(T_*\)), and let
\[
 p_+=\bar\Phi(B),\quad p_-=\bar\Phi(A),\qquad
 d_+=H_+-E[f\mid X>B],\quad d_-=H_--E[f\mid X<-A].
\]
All radii are positive finite, and both \(d_\pm>0\), since a continuous wing immediately past the endpoint still takes values below its limiting \(H_\pm>T_*\). Dominated convergence gives the necessary inequality
\[
 p_+(H_+-H_*)+p_-(H_--H_*)
 \ \ge\ p_+A[m(B)-B]+p_-B[m(A)-A]+p_+ d_++p_-d_- .          \tag{8}
\]
In particular equal tail limits are impossible for any nonconstant continuous single-well graph, without any symmetry. With unequal limits, (8) quantifies the required high-tail excess. It does not exclude all sufficiently unequal tails.

The first-moment compatibility at these endpoints is also
\[
 -A(p_++p_-)\le g(B)-g(A)\le B(p_++p_-),
\]
because the exact future first moment has no variance term and the nonnegative early measure is supported between the endpoints.

## Local verification

The companion Python 3.13.18 script uses SymPy 1.14.0 and Fraction to check the temporal Fourier factor \(2\lambda/(\lambda^2+\omega^2)\), the Gaussian mollification ratio exponent, the exponent gap in the monotone-origin argument, and the near-zero dyadic energy margin. At 80 decimal digits, mpmath 1.3.0 checks Gaussian tail-moment margins for both symmetric and unequal radii. For \(A=0.7,B=2\) the conditional quadratic means are \(1.26125087297598861\) on the right and \(2.18099867891633309\) on the left.

For scope, a continuous centered radial ramp supported at radii above 3 was centered and rescaled to satisfy the first four moment relations exactly by its covariance and variance. Its minimum is approximately \(0.9652953883\). The secant argument proves its exponential sign for every nonzero real parameter, but its full transform at \(\lambda=2\) is approximately \(0.9045882280\). Thus the crossing and low-moment constraints do not by themselves prove the radial result. The energy and support argument above does.
