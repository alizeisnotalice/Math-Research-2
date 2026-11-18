# Directional Wasserstein derivative for a bounded measurable profile

Let \(g\) be the standard normal density, let \(X,Z\) be independent standard normals, and let \(h\) be any bounded measurable real function. Take any constant \(C\ge\operatorname{ess\,sup} h\). For \(\eta\downarrow0\) put \(a=C-h\), \(Y_\eta=X+\sqrt{\eta a(X)}Z\), and let the comparison law be \(N(0,1+C\eta)\), for \(\eta\) small enough that this variance is positive. Then in extended nonnegative reals
\[
 \lim_{\eta\downarrow0}\eta^{-1}W_1(\mathcal L(Y_\eta),N(0,1+C\eta))
 =\tfrac12 |D(gh)|(\mathbb R),
\]
where the right side is infinite if \(gh\notin BV\). Mean-zero and continuity assumptions on \(h\) are unnecessary.

For the lower bound, test against a 1-Lipschitz primitive \(\psi\) of \(\phi\in C_c^\infty(\mathbb R)\) with \(\|\phi\|_\infty\le1\). The departure of \(\mathbb E\psi(Y_\eta)\) from \(\int g\psi\), divided by \(\eta\), tends to \(\int a g\phi'/2\), by Taylor and boundedness of \(a\). The corresponding expression for \(N(0,1+C\eta)\) tends to \(C\int g\phi'/2\) by differentiating its explicit density. Taking the absolute value and supremum over \(\phi\) proves the full lower bound, including divergence for \(gh\notin BV\).

If \(gh\in BV\), then \(a\) is locally BV and \(\int g\,d|Da|<\infty\) by the product rule, since \(ga=Cg-gh\in BV\) and \(a\) is bounded. Write \(b=\|a\|_\infty\). On \(\mathbb R\times(0,b)\) form the finite signed measure
\[
 M(dx,dr)=D_x[g\,1_{\{a>r\}}](dx)\,dr
 =g'(x)1_{\{a(x)>r\}}\,dx\,dr+g(x)D1_{\{a>r\}}(dx)\,dr .
\]
The weighted local coarea formula identifies the \(x\)-projection of the total variation of the second term as \(g|Da|\). The total variation of \(M\) has absolutely continuous \(r\)-projection, by integrating the finite layer variations. The signed \(x\)-projection of \(M\) is \(D(ga)\) by Fubini. Let \(S\) support \(D^sa\) and be Lebesgue-null. The restriction \(M_s=M|_{S\times(0,b)}\) has variation at most \(\int g\,d|D^sa|\) and its signed \(x\)-projection is \(gD^sa\), giving equality. The total variation of \(M_0=M-M_s\) has absolutely continuous \(x\)-projection.

Let \(K_\eta M(y)=\int p_{\eta r}(y-x)M(dx,dr)\). This is a contraction into \(L^1(dy)\). If the total variation has absolutely continuous \(x\)-projection, \(K_\eta M\) converges in \(L^1\) to the density of its \(x\)-projection: first discard \(0<r<\delta\) of vanishing total variation, then freeze \(r\) on a finite grid in \([\delta,b]\). The uniform \(L^1\) cost of replacing the heat kernel is \(O(\mathrm{mesh}/\delta)\|M\|\), by scaling and the \(L^1\) variance derivative of the Gaussian. Each frozen spatial marginal has \(L^1\) density and ordinary convolution converges strongly. This argument justifies measurable bandwidths even when the physical variance vanishes on a plateau.

Write \(G\) for the standard Gaussian CDF. With \(F_\eta\) the CDF of \(Y_\eta\), the heat equation and integration in age, tested distributionally at zero age, give
\[
 2(F_\eta-G)/\eta=K_\eta M .
\]
For example, tested against \(\phi\), the CDF pairing equals minus the pairing of the difference of laws with its primitive; the time-integrated heat equation yields \(-\int_0^b\int g1_{\{a>r\}}(p_{\eta r}*\phi') dx\,dr\), as claimed after integration by parts. Both sides are \(L^1\) since both laws have finite first moment. For \(G_\eta\) the comparison Gaussian CDF, \(2(G_\eta-G)/\eta\to Cg'\) in \(L^1\) by the explicit variance derivative, also when \(C\le0\). Thus the absolutely continuous contribution to \(2(F_\eta-G_\eta)/\eta\) converges strongly to the density of \(D^a(ga)-Cg'\,dx=-D^a(gh)\). The singular contribution has norm at most \(\int g\,d|D^sa|=|D^s(gh)|(\mathbb R)\). Combining these estimates with the CDF formula for one-dimensional Wasserstein distance gives the matching upper bound.

The \(hg'\) and \(gh'\) terms can cancel in their absolutely continuous part. Estimating variation before taking the strong limit would lose the sharp upper bound.
