# A boundary-layer transition for a steep past power trace

On the old side of the round-73 outer interface, consider the local power profile \(b_-(u)=c u^p\), \(0<u<r\), with \(c,r>0\), \(1<p<2\), and Gaussian source density \(g(a-u)\). Let
\[
 \Delta_p(h)=\int_0^r g(a-u)
          [p_{c u^p}(u)-p_{c u^p}(u-h)]\,du .
\]
If the remainder of the past source is a fixed positive distance from \(a\), it changes this difference by \(O(h)\) through the off-diagonal heat kernel estimate.

For fixed \(p\), let \(\rho=1/p-1/2\). The round-73 decomposition into loss of the second Gaussian near zero and its linear cross term gives
\[
 \Delta_p(h)\sim
 \frac{g(a)}{\sqrt{2\pi c}}\,
 \frac{(2c)^{-\rho}\Gamma(1-\rho)}{p\rho}\,h^{2\rho}.
\]
The limit in the positive leading term follows by \(u=h^{2/p} z\); its domination follows from \(z^{-p/2}\min(1,z^{-p})\), integrable at both ends for \(1<p<2\). The negative cross term is smaller as estimated in round 73.

For the transition \(p\uparrow2\), put \(k=1-p/2\). Use \(y=u^k\) to get the exact identity
\[
 \frac{k\sqrt{2\pi c}}{g(a)}\Delta_p(h)
 =
 \int_0^{r^k}\frac{g(a-y^{1/k})}{g(a)}
 \left[e^{-y^2/(2c)}-
     e^{-y^2(1-(h^k/y)^{1/k})^2/(2c)}\right]dy .
\]
If \(h^k\to H\in[0,1]\) while \(h\downarrow0\), the right side tends to
\[
 \int_0^H e^{-y^2/(2c)}\,dy .
\]
For \(0<y<1\) the first factor tends to one. The exponential difference tends to \(e^{-y^2/(2c)}1_{\{y<H\}}\). Uniform boundedness of \(g(a-y^{1/k})/g(a)\) inside the fixed core gives dominated convergence.

There is a strict next-order term when \(h^k=H\in(0,1)\) is held exactly fixed. Set \(A=H^2/(2c)\). Then
\[
 \frac{k\sqrt{2\pi c}}{g(a)}\Delta_p(h)
 =\int_0^H e^{-y^2/(2c)}dy + kH C(A)+o(k),
\]
where
\[
 C(A)=
 \int_0^1\frac{e^{-A}-e^{-A(1-x)^2}}x\,dx
 -\int_1^\infty\frac{e^{-A(1-x)^2}}x\,dx <0 .
\]
To prove this, first include the source ratio in the leading integral on \((0,H)\). That modification is \(O(h)=o(k)\) by local Lipschitz continuity of \(g\). Subtract this leading integral and set \(y=H e^{kz}\). The divided difference then converges by dominated convergence to
\[
 H\int_\mathbb R \left[e^{-A}1_{\{z>0\}}-
            e^{-A(1-e^{-z})^2}\right]dz.
\]
For \(z\ge0\) the unscaled bracket at \(y=He^{kz}\le r^k\) is bounded by \(C e^{-z}\); multiplying by \(e^{kz}\) preserves integrability uniformly for small \(k\). For \(z\le-1\) the exponent \(A e^{2kz}(1-e^{-z})^2\) grows at least as \(c_0 e^{3|z|/2}\) for \(k\le1/4\), giving a uniform integrable bound. The upper limit \((\log r^k-\log H)/k\) tends to infinity. Substitute \(x=e^{-z}\) to get the displayed \(C(A)\). The first integral is nonpositive pointwise and the second strictly positive before subtraction.

This is a real local trace asymptotic for the indicated power profile, including the nonuniformity of the limit \(h\downarrow0\) as \(p\uparrow2\). Round 73 still supplies the outer-interface implication for exact graphs under its hypotheses.
