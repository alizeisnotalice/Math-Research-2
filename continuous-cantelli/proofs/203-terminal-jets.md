# Universal exceptional point for all finite jets and for reflection

Let \(C\) be a component of \(\{f<1\}\) of a nonconstant bounded continuous normalized exact graph. Let \(q\) and \(\sigma\in(1,H]\) be the local density and common extinction time in proofs/164-common-extinction.md and proofs/166-temporal-parity.md. On \((1,\sigma]\times C\), \(q_t=q_{xx}/2\), each fixed-time slice is spatially real analytic, and all derivatives are continuous onto the terminal face on compact subintervals. The terminal slice vanishes, and every nonempty open subinterval has nonzero slices arbitrarily close to the terminal time. All space and time jets on the terminal face vanish by the equation. These are the hypotheses used below.

**Theorem.** There is a set \(E\subset C\) of cardinality at most one such that, for every \(x\in C\setminus E\) and every nonzero real polynomial \(P\), the function \(t\mapsto P(\partial_x)q(t,x)\) takes both signs arbitrarily close to \(\sigma\) from the left. The coefficients of \(P\) may depend on \(x\). For every \(m\in C\setminus E\), every \(r>0\) with \(m\pm r\in C\), and every such \(P\), each of the two temporal traces
\[
P(\partial_x)q(t,m+r)\ \pm\ P(\partial_x)q(t,m-r)
\]
also takes both signs arbitrarily close to \(\sigma\).

**Polynomial endpoint rows.** Suppose distinct \(a<b\) support eventual weak signs for the polynomials \(P_a,P_b\). Take \(1<s<\sigma\) after both sign thresholds. Write \(D=\partial_x,\ L=b-a,\ \lambda=k^2/2\) and set
\[
F(x)=\int_s^\sigma e^{\lambda(t-s)}q(t,x)\,dt,\qquad V(x)=(F(x),F'(x)/k)^T.
\]
For complex \(k\ne0\), with \(T(z)=\left(\begin{smallmatrix}\cos z&\sin z\\-\sin z&\cos z\end{smallmatrix}\right)\), integration of the equation gives
\[
F''+k^2 F=-2q(s,x),\qquad V(b)=T(kL)V(a)+W(k),
\]
where \(W=-(2/k)\int_a^b(\sin(k(b-u)),\cos(k(b-u)))^Tq(s,u)\,du\). It is \(O(1/|k|)\) uniformly in \(|\operatorname{Im}k|\le1,\ |k|\ge1\), as is the forcing up to every intermediate point. Write \(P_i(D)=A_i(D^2)+D B_i(D^2)\) and
\[
u_i(k)=(A_i(-k^2),k B_i(-k^2)),\qquad H_i(k)=\int_s^\sigma e^{\lambda(t-s)}P_i(D)q(t,i)\,dt\quad(i=a,b).
\]
Then \(H_i=u_i V(i)+J_i(k)\) with \(J_i\) a polynomial. Indeed \(D^{2j}q=2^j\partial_t^jq\), the same applies to \(Dq\), and all upper endpoint terms vanish on integration by parts. Lower endpoint terms involve only finitely many jets at time \(s\).

Write \(d_i=\deg P_i\). By opposite parity of the two components,
\(k^{-d_i}u_i(k)=c_i e_i+O(1/k)\) with \(c_i\ne0\) real and \(e_i=(1,0)\) for even \(d_i\), \(e_i=(0,1)\) for odd \(d_i\). In particular uniformly for complex \(k\) in the high positive strip,
\[
k^{-d_a-d_b}\Delta(k)=c_a c_b\det\begin{pmatrix}e_a\\ e_b T(kL)\end{pmatrix}+O(1/|k|),\qquad \Delta=\det\begin{pmatrix}u_a\\u_b T(kL)\end{pmatrix}.
\]
The leading term is a nonzero constant times \(\sin(kL+\phi)\), \(\phi\in\{0,\pi/2\}\). At real frequencies \(k_n=(n\pi-\phi)/L+o(1)\) there is a zero of \(\Delta\) for every sufficiently large \(n\), by bracketing at distances \(n^{-1/2}\). There \(u_bT(k_nL)=\gamma_n u_a\), with \(|\gamma_n|\asymp k_n^{d_b-d_a}\) and \(\operatorname{sgn}\gamma_n\) alternating since the leading vectors align and their relative sign flips with each \(n\).

Propagation yields \(H_b(k_n)-\gamma_n H_a(k_n)=O((1+k_n)^N)\) for fixed \(N\). Choose a parity making the left side have absolute value \(|H_b|+|\gamma_n|\,|H_a|\), assigning either sign to a zero trace. Each nonzero physical trace at an interior time of \((s,\sigma)\) then gives a lower bound \(c e^{k_n^2\epsilon/2}\) on its absolute transform by continuity on a compact interval after \(s+\epsilon\). This is impossible even with the polynomial factor \(|\gamma_n|\). Both physical traces thus vanish throughout \((s,\sigma)\).

**Terminal uniqueness for two polynomial endpoint conditions.** Now \(H_a=H_b=0\) for complex \(k\), so the rows \(u_a\) and \(u_b T\) determine \(V(a)\) with polynomially bounded right sides \(-J_a,-J_b-u_bW\) in the strip. Use successive rectangles with horizontal sides at \(\operatorname{Im} k=\pm1\) and vertical sides at real parts \(((n+1/2)\pi-\phi)/L\). On their boundaries \(|\sin(kL+\phi)|\ge\min(1,\sinh L)>0\). Thus \(|\Delta|\ge c|k|^{d_a+d_b}\) there for large \(n\). Cramer's rule bounds both components by \(C(1+|k|)^M\) on these boundaries. Each component is holomorphic on each rectangle by its defining time integral. The maximum principle gives the same polynomial bound on the positive real axis. Propagating to any \(a<y<b\) bounds \(F(y)\) polynomially as \(k\to+\infty\). By the compact signed Laplace support lemma of proofs/083-unimodal.md §3, applied with time reversal as in R164, a nonzero trace after \(s\), with upper support endpoint \(d>s\), has limsup exponential rate \(d-s\) in \(\lambda\). Thus \(q(t,y)=0\) throughout the terminal slab, contradicting common terminal activity. Taking \(E\) to be the set of points where some nonzero polynomial has an eventually weak sign proves the first assertion.

**Reflection.** Apply the R166 sine identity to \(u=P(D)q\) on \([m-r,m+r]\). On each parity the exponential transform of the corresponding endpoint sum or difference is \(O(1/k_n)\) for fixed initial time \(s\). If this combination has a weak sign after \(s\), it vanishes there. Use the identity for every starting time \(s'\in(s,\sigma)\). For the difference, all even-index sine coefficients of \(u(s',\cdot)\) vanish, so the slice is even about \(m\). For the sum, all odd-index coefficients vanish, making it odd. This follows from completeness and continues to the whole reflected overlap by spatial analyticity. In the even case \(D P(D)q(t,m)=0\), and in the odd case \(P(D)q(t,m)=0\). Either case forces \(m\in E\). This proves the second assertion.

The same exceptional point therefore covers all Robin traces, nonadjacent finite derivative combinations, and reflection midpoints of these fields. At an exceptional point the argument asserts existence of a weakly signed polynomial trace. It establishes no fixed parity common to different polynomials. Moving cusps and infinite endpoints require other estimates.
