# Suppressed basic harmonics on symmetric deleted-digit masks

Let the base factor as \(b=Bq\), where \(q\ge2\). Take a proper digit set \(D=R+B\{0,\ldots,q-1\}\) with \(0,b-1\in D\) and \(N=|D|>1\). Its Cantor set \(K\), viewed on the unit circle, is invariant under translation by \(1/q\). This translation permutes the first digit and leaves the rest of the expansion unchanged. Consequently distance to \(K\) has period \(1/q\). For \(q\nmid n_0\), the basic Fourier harmonic with index \(n_0\) vanishes for the positive heat damping at every fixed strength, as does the \(n_0\)-th coefficient of the digit measure.

At the resonant index \(n=n_0b^m\), the deep cells are coherent across the retained prefixes. Within each depth-\(m\) cell the residual harmonic is \(n_0\), with effective damping \(2\pi^2c n_0^2 b^{m(2-\alpha)}\). In the case \(q\nmid n_0\) this entire deep-cell contribution vanishes exactly by the periodicity just noted. The global index \(n=n_0b^m\) is still divisible by \(q\) when \(m\ge1\), so its coefficient on the full circle is not forced to vanish by the global periodicity.

Let \(G\) be the number of first-level gaps. The remaining gaps have levels \(\ell=0,\ldots,m-1\) and there are exactly \(GN^\ell\) at level \(\ell\). Their endpoints are at integer multiples of \(b^{-m}\), so both endpoint phases are one. Put \(h=(2\pi^2c n^2)^{-1/\alpha}\). For a remaining gap of length \(l\), its contribution is
\[
2h\int_0^{l/(2h)}e^{-t^\alpha}\cos(2\pi n h t)\,dt .
\]
If \(0<\alpha<2\), then every such gap has \(l\ge b^{-m}\), \(b^{-m}/h\to\infty\), and \(nh\to0\). Dominated convergence holds uniformly over these gap lengths. Thus all of them supply the same leading positive constant per gap, and the suppressed basic harmonic obeys
\[
I(n_0b^m)\sim
  2h\,\Gamma(1+1/\alpha)\,G\,\frac{N^m-1}{N-1}
  \qquad (0<\alpha<2,\ q\nmid n_0).
\]
Cancellation of the limiting Cantor harmonic therefore leaves a smaller positive contribution from the coarse gaps. By comparison, a nonzero surviving digit harmonic gives a coefficient of the larger damping-mass scale \(h^{1-\log_b N}\).

The grid uses bases \(6,8,9,10,12\), with \(q=2,2,3,2,4\) respectively, and compares \(n_0=1\) with \(n_0=q\). For each of these digit sets the \(q\)-th digit-measure Fourier coefficient is nonzero. The first digit factor reduces to a sum of roots on the smaller base \(B\). For \(B=3,4,5\) this is nonzero for the retained sets in the script. Higher digit factors lie in a strict open semicircle, and their deviations from one are summable.

The primary numerical method uses the product of digit phases and sums the complementary gaps by level, with an analytic bound on all omitted levels. A separate high-precision calculation of the displayed coarse-gap formula bypasses all deep levels and phase products when \(n_0=1\).
