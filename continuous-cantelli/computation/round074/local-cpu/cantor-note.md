# Cantor distance powers and resonant Fourier coefficients

Use a circle of circumference one. Let \(b\ge3\) be an integer and \(D\subsetneq\{0,\ldots,b-1\}\) contain \(0,b-1\), with \(N=|D|\ge2\). Let \(K\) be the deleted-digit Cantor set and \(\mu_K\) its equal-weight digit measure. For \(c,\alpha>0\) take \(v(x)=c\,\mathrm{dist}(x,K)^\alpha\), where the Euclidean and circular distances agree on \([0,1]\). Uniform stationarity of \(x\mapsto x+\sqrt{v(x)}Z\bmod 1\) requires every nonzero coefficient
\[
 I(n)=\int_0^1 e^{2\pi inx-2\pi^2n^2v(x)}\,dx
\]
to vanish. The corresponding positive damping mass is \(M(n)=\int_0^1 e^{-2\pi^2n^2v(x)}dx\).

## Complementary gaps

A first-level gap of length \(w=(d_2-d_1)/b>0\) and midpoint \(m_g=(d_1+d_2)/(2b)\), written using its two boundary digit-grid indices, produces copies at all depths \(\ell\) with length \(w b^{-\ell}\). A gap of width \(l\) and midpoint \(m\) contributes to \(I(n)\), after subtracting the undamped Fourier mode, by
\[
e^{2\pi i nm}\,l
 \int_0^1 \cos\big(\pi n l(1-t)\big)
       \operatorname{expm1}\left[-2\pi^2n^2c(l/2)^\alpha t^\alpha\right]dt.
\]
At depth \(\ell\), the sum over phases from the retained prefixes is \(\prod_{j=1}^{\ell}\sum_{d\in D}e^{2\pi i nd/b^j}\). Similar subtraction computes \(M(n)\) by adding the undamped integral one. The same absolute bound applies to either tail of levels \(\ell\ge L\). In fact \(1-e^{-z}\le z\) gives
\[
\frac{2\pi^2n^2c}{2^\alpha(\alpha+1)}
 \left(\sum_{\mathrm{first\ gaps}} w^{\alpha+1}\right)
 \frac{(N/b^{\alpha+1})^{L}}{1-N/b^{\alpha+1}} .
\]
The Cantor set itself has zero Lebesgue measure. For the numerical integral, when its variance-damping parameter is large, apply \(u=2\pi^2n^2c(l/2)^\alpha t^\alpha\) to the damped part and integrate the undamped cosine exactly. Truncating this changed variable at \(100\) discards at most \(e^{-100}\) times the gap length. The finite-level integrations are floating-point approximations, with an independent SciPy check. The bound above is an analytic bound for the omitted infinite levels.

## Subquadratic powers

Set \(s=\log N/\log b\) and \(h=(2\pi^2 c n^2)^{-1/\alpha}\). The damping mass is comparable to \(h^{1-s}\). For the lower bound, take all retained intervals at the first level whose length is at most \(h\). Their total length has this order and the distance to the Cantor set is bounded by \(h\) throughout. For the upper bound, split at that level. Each coarse gap contributes \(O(h)\), and the remaining intervals have the same total length bound. Similarly, the damping probability on \([0,1]\), obtained by normalizing \(e^{-h^{-\alpha}\mathrm{dist}(x,K)^\alpha}dx\), tends weakly to \(\mu_K\) as \(h\to0\). Indeed, for each fixed prefix depth the retained intervals have equal damping mass by self-similarity, and the normalized mass of the finitely many coarse gaps is \(O(h^s)\). The shrinking diameters of the retained intervals determine the weak limit.

Fix a positive integer \(n_0\) and take \(n=n_0b^m\). Inside every retained depth-\(m\) interval, the phase due to its prefix is exactly one. In its local coordinates \(t\in[0,1]\), the residual Fourier phase is \(e^{2\pi i n_0t}\) and the damping parameter is \(2\pi^2c n_0^2 b^{m(2-\alpha)}\). The gaps before depth \(m\) contribute at most \(O(N^m h)\) in absolute value. Relative to \(M(n)\), this is \(O((b^m h)^s)\), which tends to zero when \(0<\alpha<2\). The effective damping parameter inside the intervals then tends to infinity, so
\[
 \frac{I(n_0b^m)}{M(n_0b^m)}
 \longrightarrow \int e^{2\pi i n_0t}\,d\mu_K(t).
\]
The measure \(\mu_K\) is singular. If all of its nonconstant Fourier coefficients vanished, identification by trigonometric polynomials would give the uniform measure, a contradiction. Thus some fixed \(n_0\) yields an eventually nonzero coefficient along the resonant sequence. The proof does not require the first harmonic to be nonzero.

For \(\alpha\ge2\), at any \(a\in K\) there is the bound \(v(a+u)\le c|u|^\alpha=O(u^2)\), so the previously proved quadratic-upper-contact obstruction for uniform circular stationarity applies. Combining the statements excludes this deleted-digit distance-power family for every \(\alpha>0\).

At \(\alpha=2\) there is also an exact resonant limit readily compared with computations. For \(n=n_0b^m\), multiply \(I(n)\) by \((b/N)^m\). The retained depth-\(m\) intervals give \(I(n_0)\). The coarse level \(\ell=m-j\) contributes
\[
N^{-j}\sum_{\mathrm{first\ gaps}}2\int_0^{w b^j/2}
      e^{-2\pi^2cn_0^2 t^2}\cos(2\pi n_0 t)\,dt .
\]
Each coarse endpoint phase equals one because the endpoints are digit-grid points, and the coarse integral is the sum of the two symmetric half-gaps. The sum over \(j\ge1\) is absolutely convergent.
