# A real-fold benchmark for a divergent quartic Taylor jet

Let \(S\) have the exponential law of mean one and define
\[
 H(u)=\mathbb E(1+\gamma u^2S)^{-1},\qquad
 v(u)=d u^4H(u)^2,\qquad \gamma,d>0.
\]
This \(H\) is a positive smooth function of real \(u\). Indeed all real derivatives of \((1+\gamma u^2S)^{-1}\) are bounded by constants times powers of \(\sqrt{\gamma S}\) of the appropriate orders, so differentiation under the integral applies. The coefficient of \(u^{2k}\) in its formal Taylor jet is \((-1)^k\gamma^k k!\). The coefficient in \(H^2\) is \((-1)^k\gamma^k\sum_{j=0}^k j!(k-j)!\), of magnitude at least \(\gamma^k k!\). The quartic variance therefore has a divergent Taylor series.

For \(b>0\), consider the conditional output map \(w=u+b u^2H(u)\). At a negative input set \(t=-u\), \(x=bt\) and \(r=\gamma/b^2\). Its stationary point equation becomes
\[
 2x\,\mathbb E(1+rx^2S)^{-2}=1.
\]
If \(0<r\le1/8\), the solution is unique in \([1/2,3/5]\). The left endpoint is strictly below one. At the right endpoint, \((1+q)^{-2}\ge1-2q\) gives a value at least \(273/250>1\). The derivative of the left-hand side in \(x\) is \(2\mathbb E[(1-3rx^2S)/(1+rx^2S)^3]\), which is at least \(73/50>0\) in this bracket by \((1-3q)/(1+q)^3\ge1-6q\). The exact rational checks use nonnegative polynomial differences after clearing denominators.

If \(\kappa(b)=-1/w_{\mathrm{fold}}\), this specifies
\[
 \frac{\kappa(b)}{4b}=\frac{1}{4[x-x^2\mathbb E(1+rx^2S)^{-1}]} .
\]
The denominator is strictly larger than one for \(r>0\). Indeed \(x-x^2H(x/b)\) is maximized in the bracket and it already exceeds \(1/4\) at \(x=1/2\). The same formula with \(|b|\) applies to the absolute reciprocal fold location for negative \(b\).

Finite-order expansions of the integrands at \(r=0\) have controlled remainders from finite exponential moments. Combining them with the stationary point equation and its derivative bound gives
\[
 x=\tfrac12+\tfrac14r+\tfrac3{16}r^2+O(r^3),\qquad
 \log\frac{\kappa(b)}{4b}=-\tfrac14r-\tfrac3{32}r^2+O(r^3).
\]
The finite expansions are valid for this smooth nonanalytic \(H\).

For \(B\sim N(0,d)\), use the probability obtained by weighting its law by \(|B|^N/\mathbb E|B|^N\). On this law \(B^2/d\) is chi-square with \(N+1\) degrees of freedom, hence \(|B|/\sqrt{dN}\to1\) in probability. With any fixed cutoff \(b_0\) such that \(\gamma/b_0^2\le1/8\), the expansion and the bound \(\kappa/(4|B|)\le1\) yield
\[
 \frac{\mathbb E[|\kappa(B)|^N 1_{\{|B|\ge b_0\}}]}{4^N\mathbb E|B|^N}
   \ \longrightarrow\ \exp[-\gamma/(4d)].
\]
This is a statement about real conditional fold locations. It provides no Taylor coefficient formula for the full output density of a divergent jet.
