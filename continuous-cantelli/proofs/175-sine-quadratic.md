# Two disconnected peaks pass the quadratic budgets and violate the sine budget

Use \(s=5/4\) and \(I_L=[-L/2,L/2]\). Put \(u=x/L+1/2\) on \(I_L\). On the first support \(0\le u\le3/32\), the excess \(f-s\) rises linearly from 0 to \(3L^2/4\) over the first eighth, stays flat through seven eighths, and then declines linearly to 0. On the second support \(27/32\le u\le1\), it rises linearly from 0 to \(L^2/2\) over the first eighth, stays flat through three quarters, and declines linearly to 0. On the rest of \(I_L\) take \(f=s\), and extend continuously outside with \(f\le s\). The open strict superlevel has two components with gap \(3L/4\). The Gaussian weights below are in the \(u\) coordinate, so common factor \(L\) in the physical integrals is divided out.

The low-order future identity determines \(D_0,D_1,D_2\), the zeroth through second candidate moments of \(p_s\,dx-Q_s\) after division by \(L\). With \(F=\{u:f(x)>s\}\) and \(h=(f-s)/L^2\), these are the \(p_s\)-integrals of \(1,u,u^2\) minus the future \(g\)-integrals of \(1,u,u^2-h\), respectively. Write \(D_{\sin}\) for the first Dirichlet sine budget's right side minus left side, divided by \(L\). Equivalently,
\[
D_{\sin}=\int_0^1p_s(L(u-\tfrac12))\sin(\pi u)\,du
-\int_Fg(L(u-\tfrac12))\sin(\pi u)e^{\pi^2h(u)/2}\,du.
\]

If \(D_0>0\), \(D_0D_2-D_1^2>0\), and \(D_1-D_2>0\), then these moments give a strictly positive integral for every nonzero polynomial of degree at most two that is nonnegative on \([0,1]\). To prove it, normalize moment mass to 1 and put \(m=D_1/D_0\), \(q=D_2/D_0\). The inequalities give \(m^2<q<m\), hence \(0<m<1\). Choose endpoints \(a,b\) strictly inside \((0,1)\) with \(a<m<b\) and \((m-a)(b-m)>q-m^2\). The two-point law on these endpoints with mean \(m\) has variance \((m-a)(b-m)\). A strict mixture with \(\delta_m\) matches all three moments and charges three distinct interior points. A nonnegative quadratic can vanish at all three only when it is zero.

At \(L=1\), directed interval quadrature gives the following outward-rounded enclosures:
\[
\begin{aligned}
D_0&\in[0.2544886,0.2546264],\\
D_0D_2-D_1^2&\in[0.0140336,0.0144254],\\
D_1-D_2&\in[0.0066766,0.0082073],\\
D_{\sin}&\in[-0.0761445,-0.0738414].
\end{aligned}
\]
The physical \(x^2\) slack alone is \(D_2-D_1+D_0/4\), enclosed in \([0.0554149,0.0569800]\). The same script treats the density factors as intervals in \(L\in[0,1/4]\) in the normalized integrals. Uniformly for \(0<L\le1/4\) it gives \(D_0>0.25634\), \(D_0D_2-D_1^2>0.015047\), \(D_1-D_2>0.0023287\), and \(D_{\sin}<-0.0959443\). The two physical maximum excess heights scale as \(L^2\), so this is a family with no absolute positive threshold on peak height. Since the signs at every stated nonzero width are strict, continuity also permits small perturbations of support fractions, heights, separation and level.

## Continuous mean-one completion

For \(0<L\le1\), let \(B_L=1\) on \(I_L\), decrease linearly to zero at \(|x|=3L/4\), and vanish beyond. Let \(D=1\) on \([-2,2]\), drop linearly to zero at \(|x|=3\), and vanish beyond. If \(H=f-s\) is the two-peak excess on \(I_L\), extend \(H=0\) outside \(I_L\). For \(\beta=1/64\), define
\[
f(x)=1+(s-1)B_L(x)+H(x)
 +(1-B_L(x))[\beta/(1+x^2)-\alpha D(x)],
\]
where \(\alpha=[(s-1)E_\gamma B_L+E_\gamma H+\beta E_\gamma((1-B_L)/(1+X^2))]/E_\gamma[(1-B_L)D]\). The numerator is positive, and its bound at \(L=1\) and the lower bound \(\gamma(3/4\le |X|\le2)\) on the denominator give \(\alpha<0.471\), as enclosed by the script. Thus this is a continuous nonnegative mean-one graph with \(f\le s\) off \(I_L\) and \(\sup f=s+3L^2/4\le2\). Its strict future is exactly the stated two components, so the round-170 sine budget would apply under exact normalization. The negative sine slack excludes that normalization at once.
