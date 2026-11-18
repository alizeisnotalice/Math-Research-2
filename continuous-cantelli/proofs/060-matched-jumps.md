# Higher jump after the leading variance coefficient matches

Use signed output coordinate \(y\) around \(a\) and standard Gaussian source \(g\). Let the two one-sided positive variance germs have finite contact order greater than two. Suppose their signed Taylor coefficients agree up to degree \(m-1\), and at degree \(m\) differ by \(A\), with right minus left convention. For smooth one-sided germs the order-\((m-2)\) density jump is
\[
\frac{m!}{2}g(a)A.
\]

The common leakage coefficients follow from the fixed-order Taylor expansion away from the moving diagonal, with bound \(v(u)^{-C}\exp[-u^2/(16v(u))]\) for each fixed order. On the moving diagonal the conditional inverse expansion gives the finite series of operators \(g+\sum_{\ell\ge1} (2^\ell\ell!)^{-1}\partial_y^{2\ell}(g v^\ell)\). If the initial variance order is \(n>2\), the first possible signed difference in the \(\ell\)-th term has degree at least
\[
m+n(\ell-1)-2\ell=(m-2)+(n-2)(\ell-1).
\]
Thus the lowest difference in output jets comes entirely from \((g v)''/2\). In particular, when one-sided Peano density coefficients match to every order, their signed smooth variance jets match to every order.

The new models take even \(n=4,6,8,10,12,14\), with the same leading coefficient \(c>0\) on both sides and
\[
v_\sigma(z)=c z^n(1+b_\sigma z+d_\sigma z^2),\qquad z=|y|,\quad \sigma=\operatorname{sgn} y .
\]
If \(b_++b_-\ne0\), the first jump is \(\tfrac12 (n+1)!g(a)c(b_++b_-)\) at derivative order \(n-1\). If \(b_-=-b_+\) and \(d_+\ne d_-\), it is \(\tfrac12 (n+2)!g(a)c(d_+-d_-)\) at derivative order \(n\).

Quantitatively, the normalized diagonal correction obtained by subtracting \(g(a+\sigma z)\) and dividing by \(g(a)z^{n-2}\) has coefficients
\[
\begin{split}
C_0&=\tfrac12 cn(n-1),\\
C_{1,\sigma}&=\tfrac12 cn(n+1)(b_\sigma-\sigma a),\\
C_{2,\sigma}&=\tfrac12 c(n+1)(n+2)
 \big[d_\sigma-\sigma a b_\sigma+(a^2-1)/2\big]
 +210c^2\mathbf1_{n=4}.
\end{split}
\]
The \(210c^2\) term comes from the second heat coefficient. It cancels between the two sides in this grid, since \(c\) is common.

The 87-digit diagonal integral uses the conditional Gaussian scaling \(u=\sigma(z+\sqrt{v_\sigma(z)}t)\) and subtracts the constant Gaussian directly using \(\operatorname{expm1}\). A binary64 adaptive integral computes the same symmetric sum by even pairing, expanding the mean and divided antisymmetric exponent in stable formulas. This second route validates the unamplified diagonal correction. The leading lower coefficients are removed in the high-precision route before extracting the higher-order jump. The computation compares the diagonal contributions, with common leakage handled by the expansion above.
