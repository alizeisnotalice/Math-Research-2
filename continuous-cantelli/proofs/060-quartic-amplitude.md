# Expanding both asymmetric real-saddle coefficients to fourth order

Use the round-51 diagonal Wick group in scaled variables, with \(m\) fixed and \(k/N\) staying in a compact interior interval. For \(H(u)=\exp(\alpha u+\beta u^2+\gamma u^4)\), the log-radius phase near the saddle of either sign is
\[
 F_a(X)=F_0(X)+a A(X)+a^2 B(X)+a^4 C(X),
 \qquad
 A=2\nu\sigma\alpha S,\quad B=2\nu\beta S^2,\quad C=2\nu\gamma S^4.
\]
Here \(X=(x,y)\), \(S=e^x+e^y\), \(\nu=k/N\), \(a=(cN)^{-1/(2m-2)}\), and \(\sigma\in\{-1,1\}\). The translated Gaussian source can remain an amplitude, as in round 51.

All derivatives below are evaluated at the reference saddle \(X_0\). Set \(J=F_0''\), \(v=-J^{-1}A'\) and \(L=B'+A''v+\tfrac12 F_0'''[v,v,\cdot]\). The shifted radius is \(X_0+av-a^2J^{-1}L+O(a^3)\). The coefficients of its stationary value are
\[
\begin{aligned}
 [a]F_a&=A,\\
 [a^2]F_a&=B-\tfrac12 A'^T J^{-1}A',\\
 [a^3]F_a&=B'[v]+\tfrac12 A''[v,v]+\tfrac16 F_0'''[v,v,v],\\
 [a^4]F_a&=C+\tfrac12 B''[v,v]+\tfrac16 A'''[v,v,v]
                +\tfrac1{24}F_0^{(4)}[v,v,v,v]-\tfrac12 L^T J^{-1}L .
\end{aligned}
\]
Taylor expansion gives these formulas directly. The unspecified cubic radius shift cancels from the fourth coefficient because \(Jv+A'=0\). The remaining quadratic shift contribution is the minimum of \(J[w,w]/2+L[w]\), equal to \(-L^T J^{-1}L/2\). Under \(\sigma\mapsto-\sigma\), \(A\) and \(v\) change sign and \(L\) is unchanged. Thus the odd coefficients change sign and the even coefficients agree through fourth order. The quartic log-amplitude enters the fourth coefficient exactly as \(2\nu\gamma S^4\).

The finite coefficient test retains the Gaussian factor explicitly after division by \(g(0)\). In \(Z=W+T\), its resulting amplitude is \(e^{A_1 Z+B_2 Z^2+C_4 Z^4}\) with
\[
 A_1=a(2k\alpha-x_0),\quad B_2=a^2(2k\beta-1/2),\quad C_4=2ka^4\gamma .
\]
Only its even coefficients contribute to the even-\(N\) group, so replacing \(A_1\) by \(|A_1|\) is exact. When \(B_2,C_4\ge0\) the log recurrence for coefficients \(b_j\) uses
\[
 j b_j=|A_1|b_{j-1}+2B_2b_{j-2}+4C_4b_{j-4}.
\]
An independent check convolves the quadratic exponential with \(e^{C_4 Z^4}\). Selected cases also use shifted-radius two-dimensional Cauchy quadrature, including exact linear cancellation \(A_1=0\).
