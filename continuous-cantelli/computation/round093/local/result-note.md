# Fixed-interval Brownian cap and the optimal quadratic budget

Use the hypotheses of Round 84. Let \(I=(-A,B)\) be the limiting early interval at time \(T\), and let \(\tau_I\) be the exit time of standard Brownian motion started at zero from this fixed interval. The moving survival event up to \(t<T\) is contained in \(\{\tau_I>t\}\), up to events of probability zero at the boundary. Hence the measures \(Q_t\) of Round 84 satisfy
\[
 0\le Q_t(dy)\le p_I(t,0,y)\,dy,
\]
where \(p_I\) is the Dirichlet heat kernel with generator \(\frac12\partial_{yy}\). Their limit \(Q\) as \(t\uparrow T\) satisfies the same domination at \(T\). One way to take the limit is to use continuity of paths, bounded convergence of the decreasing survival indicators, and continuity of the fixed-interval heat kernel. The future-source moment identities of Round 84 give
\[
 M=Q(I)=\bar\Phi(A)+\bar\Phi(B),\qquad
 N=\int_I y\,Q(dy)=\phi(B)-\phi(A).
\]

Write \(p(y)=p_I(T,0,y)\), \(C=\int_I p(y)\,dy\), and \(w(y)=(y+A)(B-y)\). The Round 84 excess budget has a nonnegative additional term \(\int_I w\,dQ\). With \(E=\int_{\mathbb R\setminus I}(f(x)-1-T)\gamma(dx)\) and \(E_{84}=A(\phi(B)-B\bar\Phi(B))+B(\phi(A)-A\bar\Phi(A))\), the exact moment identity reads \(E=E_{84}+\int_I w\,dQ\). This follows by expanding \(w\) and applying the mass, first and second moments computed in Round 84. The boundary carries no Gaussian or Brownian survival mass.

If \(C<M\), the cap is infeasible. Otherwise let \(F(u)=\int_{-A}^u p(y)\,dy\). The smallest first moment achievable at mass \(M\) is that of the leftmost cap interval of mass \(M\). The largest is that of the rightmost cap interval. These facts follow by subtracting the threshold from the linear objective and checking the signs on the two sides. For \(C=M\), the only candidate is the entire density \(p\). For \(C>M\), let \(q\in[0,M]\) and choose
\[
 l=F^{-1}(q),\qquad u=F^{-1}(C-M+q).
\]
The first moment of \(Q_*=p\,\mathbf1_{(-A,l)\cup(u,B)}\,dy\) is continuous in \(q\) and strictly decreasing because its interior derivative is \(l-u<0\). It spans precisely the feasible range. If it equals \(N\), \(Q_*\) minimizes the quadratic cost over the relaxed feasible class. Indeed \(w(y)-(y-l)(u-y)\) is affine and integrates to the same value under all candidates. The remaining quadratic is negative on the endpoint strips, where \(Q_*\) saturates the cap, and positive on the central interval, where \(Q_*\) vanishes. This proves
\[
 E\ \ge\ E_{84}+\Delta(A,B,T),\qquad
 \Delta(A,B,T)=\int_{(-A,l)\cup(u,B)} w(y)p(y)\,dy>0.
\]
Strict positivity follows since \(M>0\), the cap is absolutely continuous, and \(w>0\) in the interval. The fixed-interval cap and two moments alone cannot improve \(\Delta\). Moving survival and exact graph structure can impose more restrictions.

## Numerical reproduction

Run **python bounded_heat.py --duration 30** followed by **python validate.py**. Each grid point is computed twice using integrated Gaussian images and eigenfunction primitives. The separate 55-digit check solves the two moment equations simultaneously at representative points. These finite computations calibrate the proved inequality and the relaxed feasibility conditions. The numerical classification uses ordinary floating arithmetic and the stated root tolerance. The grid is deterministic, stratified via irrational increments in the logarithms of \(A,B\in[0.24,3.6]\) and \(T/(A+B)^2\in[0.002,0.36]\). The duration controls how many fully completed grid points each worker retains.

## Large symmetric boxes

For \(A=B=R\to\infty\) with \(0<T<1\) fixed, the improvement factor satisfies
\[
 \frac{\Delta(R,R,T)}{E_{84}(R,R)}\sim (1-T)R^2.
\]
Indeed, the cap is feasible for all sufficiently large \(R\) since its mass tends to one and the required first moment is zero. The symmetric threshold \(c_R\in(0,R)\) satisfies \(\int_{c_R}^R p(y)dy=\bar\Phi(R)\). For \(|y|<R\), reflection at each boundary bounds the probability density of paths hitting either boundary by the sum of their separately reflected densities. Thus
\[
 0\le p_T(y)-p(y)\le p_T(2R-y)+p_T(2R+y),
\]
where \(p_T\) is the free Gaussian density at time \(T\). On \(0\le y\le sR\) for fixed \(s<1\), the ratio of the right side to \(p_T(y)\) is at most \(2e^{-2(1-s)R^2/T}\). If \(0<v<\sqrt T\), choose \(v<s<\sqrt T\). Then \(\int_{vR}^{sR}p(y)dy/\bar\Phi(R)\to\infty\) by the Gaussian exponential rates. For \(\sqrt T<v<1\), the free upper bound gives \(\int_{vR}^R p(y)dy/\bar\Phi(R)\to0\). Hence \(c_R/R\to\sqrt T\), and under the normalized cap restricted to the optimal strips \(|y|/R\to\sqrt T\) in probability. The quadratic cost divided by \(R^2M\) tends to \(1-T\) by boundedness. Finally \(E_{84}(R,R)/M\to1\). Substitution \(y=R+u/R\) writes the ratio as the quotient of integrals of \(ue^{-u-u^2/(2R^2)}\) and \(e^{-u-u^2/(2R^2)}\), both converging to one. This asymptotic analyzes only the relaxed cap and its moment constraints.

At \(T=1\) the gain has a finite limit. With \(u=R(R-y)\), the right-edge cap satisfies
\[
 \frac{p_I(T,0,R-u/R)}{p_T(R)}
 =2e^{-u^2/(2TR^2)}\sum_{n\in\mathbb Z}
 e^{-[(4n+1)^2-1]R^2/(2T)}\sinh\big((4n+1)u/T\big).
\]
This follows by pairing the standard images at opposite distances from the right edge. For fixed \(T\), bounded \(u\ge0\), the \(n=0\) term dominates uniformly, with exponentially small relative error for \(u>0\) extending continuously to \(u=0\). For example, on \(0<u\le H\), \(|\sinh(ku/T)/\sinh(u/T)|\le |k|e^{|k|H/T}\), and the Gaussian prefactors bound the sum of every \(n\ne0\) term exponentially in \(R^2\).

At \(T=1\), the right-edge strip mass in units \(p_1(R)/R\) tends to \(\int_0^D 2\sinh u\,du\) when its rescaled width tends to \(D\). The required normalized strip mass tends to 1. Monotonicity and convergence on each bounded interval therefore force \(D=\operatorname{arcosh}(3/2)\). The normalized quadratic cost tends to
\[
 \frac{\Delta}{E_{84}}\ \longrightarrow\
 \int_0^D 4u\sinh u\,du=6\operatorname{arcosh}(3/2)-2\sqrt5 .
\]
For fixed \(T>1\), the required normalized strip mass instead equals
\[
 \frac{R\bar\Phi(R)}{p_T(R)}
 \ \sim\ \sqrt T\,e^{-(1-1/T)R^2/2}\to0 .
\]
The rescaled width \(D_R=R(R-c_R)\) tends to zero. The paired image formula then gives \(p_I(T,0,R-u/R)/p_T(R)\sim 2u/T\) uniformly for \(0<u\le D_R\). The strip mass implies \(D_R\sim T^{3/4}e^{-(1-1/T)R^2/4}\). Over either strip the quadratic cost is \(2u-u^2/R^2\), so its mean is asymptotic to \(4D_R/3\) and
\[
 \frac{\Delta}{E_{84}}\ \sim\
 \tfrac43 T^{3/4}\exp[-(1-1/T)R^2/4].
\]
The numerical scripts **symmetric_asymptotic.py** and **symmetric_edge.py** check the three regimes without ordinary floating-point underflow.

## A strictly certified new exclusion

Since \(f\ge0\) and \(\mathbb Ef=1\), the actual excess also satisfies \(E\le 1-(1+T)M\). In the symmetric box \(A=B=T=1\), the cap permits the required mass and zero first moment. The original budget passes the comparison against this normalization upper bound. The new budget nevertheless contradicts it. To certify the claim, write \(J=\int_{-1}^1 y^2p(y)dy\). Then \(\int (1-y^2)dQ\ge M-J\), so \(E\ge 2\phi(1)-J\). The first spectral term of \(J\) is \(e^{-\pi^2/8}(4/\pi-32/\pi^3)\), and the remaining terms have absolute value bounded by \((2/3)\sum_{n\ge1}e^{-(2n+1)^2\pi^2/8}\). The same tail controls the survival mass with prefactor 2. The script **exact_example.py** compares these bounds using rational arithmetic with \(3.14<\pi<22/7\), alternating integral bounds for the Gaussian tail, and geometric remainders for positive exponential series. The resulting strict exclusion is a necessary-condition obstruction at this prescribed geometry and time.
