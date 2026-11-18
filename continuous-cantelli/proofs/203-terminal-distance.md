# Recent mark distance and late caloric traces

Use the signed heat potential \(Q_t=P_t\delta_0-\int_{f(x)\le t}P_{t-f(x)}\delta_x\,\gamma(dx)\) of a bounded continuous normalized exact graph. Fix a spatial point \(y\) with \(f(y)<\sigma\). By continuity it has a source-free spacetime neighborhood at the terminal time. Assume in addition every left time derivative of \(q(\cdot,y)\), including order zero, vanishes at \(\sigma\). This is guaranteed in the setting of R164/R166 for \(y\) in its low component. At other points it is an assumption. For \(0<h<\sigma/10\) set
\[
d_h(y)=\operatorname{dist}\big(y,\{x:\sigma-8h<f(x)<\sigma\}\big).
\]
If \(d_h(y)\ge\sqrt h\), then for every \(0<\eta<1/4\),
\[
 |q(\sigma-h,y)|\le C_\eta h^{-1/2}
 \exp\big(-(1/4-\eta)d_h(y)^2/h\big).
\]
If the set is empty the trace at \(\sigma-h\) is zero. The hypothesis \(f(y)<\sigma\) alone does not imply the stated terminal jet condition.

**Proof.** Choose \(\tau<1\) sufficiently close to one that \(7\tau>e\) and \(\tau/[2(1+\tau)]>1/4-\eta\). For \(u>0\) use the Cauchy circle \(|z-u|=\tau u\) in the time variable of the Gaussian. Since \(\operatorname{Re}(1/z)\ge1/((1+\tau)u)\),
\[
 |\partial_u^n p_u(y-x)|\le C_\tau n!(\tau u)^{-n}u^{-1/2}
 e^{-|y-x|^2/(2(1+\tau)u)}.
\]
Recent sources are at distance at least \(d=d_h(y)>0\); their time kernels extend smoothly by zero at nonpositive ages. In the left terminal interval, graph marks at time \(\sigma\) or later contribute zero, including all left jets. Sources marked at or before \(\sigma-8h\), and the starting kernel \(p_t\), have ages at least \(7h\) for \(\sigma-h\le t\le\sigma\). Thus, after division by \(n!\) and multiplication by \(h^n\), their derivative bound is \(C_\tau h^{-1/2}(7\tau)^{-n}\). For the recent sources put \(D=d^2/(2(1+\tau))\). The corresponding bound is
\[
 C_\tau(h/\tau)^n \sup_{u>0}u^{-n-1/2}e^{-D/u}
 = C_\tau(h/\tau)^n((n+1/2)/(eD))^{n+1/2}.
\]
These are uniform bounds, permitting differentiation under the Gaussian source probability measure. Take \(n=\lfloor \tau D/h\rfloor\). When \(\tau D/h\ge1\), \((n+1/2)/(\tau D/h)\le1+1/(2n)\), so the last bound is \(O_\tau(h^{-1/2}e^{-n})\). The same holds for the old sources. Taylor's remainder at \(\sigma\) proves the estimate. If \(\tau D/h<1\), use the order-zero estimate with \(d\ge\sqrt h\), increasing the constant. In the empty case use the old bound alone with \(n\to\infty\).

For a varying observation point \(y_h\), the same uniform estimate applies if the jet hypothesis holds at each point and \(f(y_h)<\sigma\). In particular \(d_h(y_h)\ge c h^\alpha,\ 0\le\alpha<1/2\), gives an attenuation factor \(\exp(-c' h^{-(1-2\alpha)})\). For a fixed cut and \(d_h\ge c h^{-\rho}\), \(\rho>0\), the factor is \(\exp(-c'h^{-(1+2\rho)})\). The terminal normalized Laplace integral
\[
 \int_s^\sigma e^{-k^2(\sigma-t)/2}|q(t,y)|\,dt
\]
is then \(O(\exp(-c'' k^{(1+2\rho)/(1+\rho)}))\), by minimizing \(k^2h/2+c'h^{-(1+2\rho)}\) and retaining half the exponential for integration. One needs the trace integrable before the terminal regime, as holds on the R166 cylinder. If \(d_h\to\infty\) at a fixed cut, take \(d_h\ge M\) sufficiently near the terminal time for each fixed \(M\). Minimizing \(k^2h/2+c'M^2/h\) then gives \(o(e^{-Ak})\) for every fixed \(A>0\). None of these upper bounds alone contradicts common terminal activity.

**Spatial phase bound at a source-free time slab.** Given \(1<s<\sigma\), on an R166 component set \(F=\int_s^\sigma e^{k^2(t-s)/2}q(t,\cdot)\,dt\) and \(V=(F,F'/k)^T\). For real \(k>0\) and \(a,y\in C\),
\[
 \|V(y)-T(k(y-a))V(a)\|_2\le 2\|Q_s\|_{\mathrm{TV}([\min(a,y),\max(a,y)])}/k\le4/k,
 \qquad T(z)=\begin{pmatrix}\cos z&\sin z\\-\sin z&\cos z\end{pmatrix}.
\]
Indeed \(F''+k^2F=-2q(s,\cdot)\); the forced convolution uses a unit vector. The bound is uniform in interval length. If \(A=\|V(a)\|\) and \(4/(kA)=\epsilon<1/3\), a whole nodal slot about any zero of the first component of the rotating vector contains one simple zero of \(F\) within angular distance \(\arcsin\epsilon\). More precisely one can use the interval between consecutive rotating extrema around that zero: away from the indicated core the first component has the sign of the rotating wave, and inside the core the derivative has constant strict sign since \(\sqrt{1-\epsilon^2}>\epsilon\). Only complete such slots inside the component are covered.
