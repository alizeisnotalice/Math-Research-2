# Independent audit of the single-valley endpoint obstruction

This audit takes as input the early positivity and support in R084 and the early symmetry in R088. Let \(0<t<T\) and \(I_t=[-b_t,b_t]=\{f\le t\}\). Write \(g\) for the standard Gaussian density and \(\bar\Phi\) for its survival function. The future transform is
\[
 \int e^{zy}Q_t(dy)=\int_{|x|>b_t}e^{zx-(f(x)-t)z^2/2}\gamma(dx).
\]
Because the normalized graph is bounded, Gaussian domination on compact parameter intervals justifies differentiation at zero. Thus the mass is \(2\bar\Phi(b_t)\) and the second moment is \(\int_{|x|>b_t}(x^2-f(x)+t)\gamma(dx)\). Positivity and support give
\[
 \int_{|x|>b_t}(f-1-t)d\gamma\ge
 \int_{|x|>b_t}(x^2-1-b_t^2)d\gamma
 =2b_t[g(b_t)-b_t\bar\Phi(b_t)].
\]
The last step is Gaussian integration by parts. Its bracket equals \(\int_{b_t}^\infty(x-b_t)\gamma(dx)>0\). The factor of two is from both tails. The asymmetric R084 formula follows from testing \(-(y+A_t)(y-B_t)\), using the first moment \(\int_{\{f>t\}}x\gamma(dx)\) from the same transform. The gap becomes
\[
 A_t[g(B_t)-B_t\bar\Phi(B_t)]+B_t[g(A_t)-A_t\bar\Phi(A_t)].
\]

Let \(L=1+T\). The early increasing radius has a positive limit \(b\), and it is finite because both tail limits exceed \(T\). Dominated convergence in the last inequality gives the displayed endpoint budget of the attachment. It only uses points strictly beyond \(b\), so an interval where \(f=T\) starting at \(b\) remains in the future outer tail budget. Endpoint atoms are absent from the Gaussian source.

The third coefficient of the defining transform is \(E[X^3]/6-E[Xf]/2\), so \(E[Xf]=0\). Early symmetry cancels the interior. With \(R\sim\mathcal L(X\mid X>b)\), \(F_\pm(r)=f(\pm r)\), and \(m=ER>0\), the equality becomes \(E[RF_+]=E[RF_-]\). If \(F_-\le L\) and \(F_+\) is nondecreasing, then the independent-copy covariance identity gives
\[
 mEF_+\le E[RF_+]=E[RF_-]\le mL.
\]
Also \(EF_-\le L\). Therefore the endpoint excess is nonpositive, contrary to the strictly positive gap. The proof also works with the weaker conditional assumption \(\operatorname{Cov}(R,F_+)\ge0\). No covariance sign on the low wing is needed.

For the quantitative version at fixed \(t<T\), replace \(b\) by \(b_t\), let \(Y\) have law \(Q_t/(2\bar\Phi(b_t))\), and put \(C_+=\operatorname{Cov}(R,F_+)\). The first calculation gives
\[
 E(F_-+F_+)=2(1+t+b_tm-EY^2),
\]
and the cubic equality gives \(E[RF_-]/m=EF_++C_+/m\). Combining them yields the attachment's exact formula with the positive sign of \(C_+/m\):
\[
 E[(1+R/m)(F_--1-t)]
 =2b_t(m-b_t)+2(b_t^2-EY^2)+C_+/m.
\]
The second term is nonnegative. These identities verify both coefficients and the direction of the contradiction.

## Calculation

Run `/Users/zzh/.venvs/main/bin/python check_endpoint.py`. Python 3.13.18 with SymPy 1.14.0 verifies the second derivative, cubic coefficient, and exact weighted residual with zero symbolic error. mpmath 1.3.0 at 65 decimals compares the two-tail polynomial integral with its claimed decomposition at two symmetric and two asymmetric positive endpoints. The largest absolute residual is about \(1.2\times10^{-65}\). Exact integration by parts proves positivity and the identity for all endpoints.
