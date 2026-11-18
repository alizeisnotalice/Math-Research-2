# Reverse Brownian graph hits: strict ordered minors and a nontrivial kernel

Let \(f\) be continuous and nonnegative and let \(I\) be an open interval with \(f<\sigma\) in \(I\) and \(f\ge \sigma\) at each finite endpoint. From \(y\in I\) run \(Y_u=y+B_u\) at backward absolute time \(\sigma-u\). Stop at the first \(u\in[0,\sigma]\) with \(\sigma-u\le f(Y_u)\), and let \(K_{\sigma,y}\) be the law of the stopped position. The path cannot leave \(I\) before stopping, since a finite endpoint has graph time at least \(\sigma\). Continuity gives equality with the graph time at the stop and implies that no exit is possible at an endpoint of \(I\) after positive backward time while still alive.

## Ordered minors

For \(y_1<\cdots<y_m\) in \(I\), take bounded nonnegative Borel weights \(h_1,\ldots,h_m\) carried by disjoint ordered Borel subsets of \(I\). In particular every point in the carrier of \(h_j\) is strictly less than any point in the carrier of \(h_{j+1}\). Then
\[
 \det\!\left[\int h_j\,dK_{\sigma,y_i}\right]_{i,j=1}^m\ge0.
\tag{1}
\]
If each \(h_j\) is bounded below by a positive constant on some open interval in its carrier, the determinant is strictly positive.

To prove this, use \(m\) independent Brownian paths and freeze each coordinate at its graph hit. The determinant expansion is the expectation of
\(\sum_{\pi}\operatorname{sgn}(\pi)\prod_i h_{\pi(i)}(Y_i^*)\),
where \(Y_i^*\) is the terminal position of path \(i\). On configurations with positive product, all terminal positions are distinct. At the first meeting of two frozen-or-live paths, neither can already have frozen. Indeed, at the frozen position \(x\) its graph time equals its earlier absolute time, so a live path reaching it at a strictly later reverse time would already have hit the graph. If two paths freeze at the collision they have identical terminal positions and the product is zero. Thus on contributing collision configurations the first meeting is a meeting of live paths strictly before both of their stops.

At the earliest meeting exchange the future Brownian increments of its two paths, choosing the lexicographically first ordered pair if there is a tie. The vector strong Markov property at that common time makes this a measure-preserving involution. The stopping rule depends only on the absolute time and the spatial path, so the two stopped futures, including their frozen parts, are exchanged. The pair selection and meeting time are unchanged because the initial paths are unchanged. Compose the corresponding two positions of \(\pi\) with this transposition. This preserves the product and reverses its sign. The cancellation leaves only nonmeeting paths. Their terminal positions are ordered, so the only surviving permutation is the identity. This proves (1) and expresses the determinant as the positive-weight probability integral of nonmeeting paths.

For strictness, choose a point in each of the open intervals on which the corresponding \(h_j\) is bounded below. The ordered linear interpolation from the starts to these targets stays inside a fixed compact subinterval of \(I\); all gaps remain uniformly positive. By continuity the supremum of \(f\) on a small compact tube about these interpolations is strictly less than \(\sigma\). Make the interpolation time small enough that the absolute time stays above that supremum until all the target intervals are reached. Thereafter let each prescribed path stay fixed at its target. Independent Brownian motion places positive probability in sufficiently small uniform tubes about the whole collection of paths over \([0,\sigma]\). On that event they reach the designated target intervals while alive, then stay in their disjoint intervals until their hits. This event has positive nonmeeting weight.

## Consequence for a Cantelli extinction plateau

Use \(W(t,x)=\int_t^H Q_s(x)\,ds\) and \(w(x)=W(f(x),x)\). Let \(I\) be a component of \(\{E>f\}\), where the local extinction function \(E\) is constant with value \(\sigma\) on \(I\). At each finite endpoint \(f=\sigma\). The late slice vanishes, so by the reverse hitting representation
\[
 W(\sigma,y)=\int_I w(x)K_{\sigma,y}(dx)=0
 \qquad (y\in I).
\tag{2}
\]
If continuous \(w|_I\) is nonzero, it has sequences of arbitrary finite length with alternating nonzero signs. Otherwise its set of nonzero values splits into finitely many disjoint ordered sign blocks, discarding the zeros. Each nonempty block contains an open interval on which \(|w|\) is bounded below. Use \(h_j=|w|1_{\text{block }j}\) in (1) and any equally numerous distinct ordered starting points. Equation (2) makes the signs of the blocks into a nonzero null vector of a strictly positive determinant, which is impossible.

If \(I\) meets \(\{f<1\}\), then it contains a point with \(w>0\), so the alternating sign conclusion always holds there. A finite accumulation point of alternating signs has zero trace by continuity and has \(f>1\): in a nonconstant exact graph the trace is strictly positive wherever \(f\le1\), by the positive occupation formula, even at \(f=1\). In particular a bounded \(I\) forces a finite accumulation point. The argument on an unbounded \(I\) allows the alternations to escape to infinity.

## Strict ordered minors do not imply injectivity even for this Brownian kernel

Fix \(0\le b<1\), \(T>0\), and choose the continuous Lipschitz clock
\[
 f(x)=b+\min\bigl(T,(|x|-1)_+\bigr),\qquad
 \sigma=b+T,\qquad I=(-1-T,1+T).
\]
We construct a continuous nonzero \(w\) on the real line with \(K_{\sigma,y}w=0\) for all \(y\in I\), which is strictly positive on \((-1,1)\). The clock in this example serves only as a barrier for \(K\).

Let \(t\) measure time above \(b\), \(r=1+t\), \(s=t/r\), \(z=x/r\), and \(S=T/(1+T)\). On the forward expanding interval \(|x|<1+t\), set
\[
 U(t,x)=r^{-1/2}e^{-x^2/(2r)} V(s,z).
\tag{3}
\]
The chain rule gives
\((\partial_t-\tfrac12\partial_{xx})U
 =r^{-5/2}e^{-x^2/(2r)}(\partial_s-\tfrac12\partial_{zz})V\).
Thus it suffices to take a solution \(V\) of the heat equation on \((-1,1)\), with \(V(0,z)=\phi_1(z)=\cos(\pi z/2)\), left boundary zero, and right boundary a smooth compactly supported function \(g(s)\) inside \((0,S)\). Such a solution is given by subtracting the lift \(g(s)(z+1)/2\) and using the homogeneous Dirichlet heat semigroup with the resulting smooth forcing. The initial and boundary data are compatible.

Let \(\phi_n(z)=\sin(n\pi(z+1)/2)\) and \(\lambda_n=n^2\pi^2/8\). Integration against \(\phi_n\) gives
\[
 \big\langle V(S),\phi_n\big\rangle
  =e^{-\lambda_n S}
     \left(1_{\{n=1\}}+c_n\int_0^S e^{\lambda_n s}g(s)\,ds\right),
 \qquad c_n=-\tfrac12\phi_n'(1)\ne0.
\tag{4}
\]
The eigenfunctions are normalized because the interval has length 2. To make every expression vanish, it remains to find a smooth compactly supported \(g\) whose exponential moments vanish at every \(\lambda_n,\ n\ge2\), and whose \(\lambda_1\) moment is nonzero. The latter can then be scaled to \(-1/c_1\).

Choose \(d_n>0\), \(n\ge2\), with \(\sum d_n=D<S/2\) and \(\sum e^{-\lambda_n d_n}<\infty\). For the tail take \(d_n=C\log(n+1)/n^2\) with fixed sufficiently large \(C\); by starting sufficiently late its sum is arbitrarily small. Take the finitely many earlier \(d_n\) arbitrarily small and positive. Put \(a_n=e^{-\lambda_n d_n}\) and
\[
 \mu_n=\frac{\delta_{d_n}-a_n\delta_0}{1+a_n}.
\]
The translate of \(\mu_n\) by \(-d_n\) has total variation 1 and differs from \(\delta_0\) by total variation at most \(2a_n\). The convolution products of the translates therefore converge in total variation by telescoping. Translate their limit by \(D\) to obtain a finite signed measure \(\mu\) supported on \([0,D]\); the original finite convolutions converge weakly to \(\mu\) on this compact interval. At each \(\lambda_n,\ n\ge2\), its Laplace transform is zero due to the corresponding factor. At \(\lambda_1\) the transform is strictly positive. Indeed it is
\[
 e^{-\lambda_1 D}
  \prod_{n\ge2}
 \frac{1-e^{-(\lambda_n-\lambda_1)d_n}}{1+e^{-\lambda_n d_n}}>0.
\]
The differences from one of the tail factors are summable, and all factors are strictly positive. Convolve \(\mu\) with a nonnegative nonzero smooth bump supported in \((0,S-D)\); its compactly supported smooth density \(g_0\) lies inside \((0,S)\) and has the same moment vanishing and nonvanishing properties under the Laplace transform. Reflect \(s\mapsto S-s\) and scale to obtain the requisite \(g\) in (4). Completeness of the sine basis now gives \(V(S,\cdot)=0\) everywhere.

For \(|x|\le1\), define \(w(x)=U(0,x)\). For \(1<|x|<1+T\), define \(w(x)=U(|x|-1,x)\), evaluating the lateral boundary value. Set \(w=0\) for \(|x|\ge1+T\). Continuity follows from the compatibility at time 0 and the compact support of the lateral control. On the interior of the early core the trace is \(e^{-x^2/2}\cos(\pi x/2)>0\). Apply the same bounded reverse-time Brownian representation to the caloric function (3) on the expanding open interval. The first graph hit supplies exactly its initial or lateral boundary value \(w\); the initial core clock is \(b\), a harmless shift in absolute time. Since \(U(T,\cdot)=0\) everywhere in \(I\), this proves \(K_{\sigma,y}w=0\) for every starting point. The kernel therefore has strictly positive ordered minors in the sense of (1) and nevertheless has a nontrivial continuous null function. The Cantelli trace obeys the additional global source equation, which is absent from this example.

The local SymPy 1.14.0 check verified the three identities for the coefficient, spatial coordinate and time coordinate in the transformation (3). The collision cancellation and the heat control construction are analytic.
