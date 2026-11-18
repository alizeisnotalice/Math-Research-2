# Audit of the quartic jet construction and its finite-transform extension

Write $g=p_1$ for the Gaussian source. On $|u|<R$ let $v(u)=u^4w(u)$, where $w$ is smooth, $w(0)=c>0$, and $c/2\le w\le3c/2$. Put $S_R(y)=\int_{-R}^R g(u)H_{v(u)}(y-u)\,du$ with $H_s=\tfrac12\int_0^s p_t\,dt$.

## The potential is smooth without exactness

Select a smooth spatial cutoff $\zeta$ supported in $(-R,R)$ and equal to one near 0. The Fourier transform of the $\zeta$-localized Gaussian output is
$$
J(\xi)=\int\zeta(u)g(u)e^{i\xi u-\xi^2u^4w(u)/2}\,du.
$$
For $\Lambda=|\xi|\ge1$, set $z=\sqrt\Lambda\,u$. The rescaled amplitude is $\zeta(z/\sqrt\Lambda)g(z/\sqrt\Lambda)\exp[-z^4w(z/\sqrt\Lambda)/2]$. Every fixed number of its $z$ derivatives has a uniform bound by a polynomial in $|z|$ times $e^{-c'z^4}$. The rescaled oscillation has frequency $\sqrt\Lambda$. After $N$ integrations by parts, $J(\xi)=O_N(\Lambda^{-(N+1)/2})$. Fourier inversion therefore gives a smooth density representative for the localized output.

The spatial complement inside $[-R,R]$ is separated from 0 and has analytic Gaussian density near 0. In distributions $H_s''=p_s-\delta_0$ for $s>0$, so $S_R''=D_R-g$ near 0, where $D_R$ is the smooth core output representative. The potential is continuous because $|H_s|\le\sqrt{s/(2\pi)}$, and therefore is itself smooth near 0. The raw integrand evaluated at $y=0$ may omit a finite diagonal contribution. Distributional differentiation identifies the correct smooth representative.

The near-diagonal substitution $u=y+y^2z$ gives the factor $y^4$ times the $z$ integral of $g(u)H_{v(u)/y^4}(z)$. On $|u-y|\le |y|/2$, a fixed number $j$ of $u$ derivatives of the amplitude costs $O_j(|y|^{-j})$ times a fixed polynomial-Gaussian bound. This gives the real expansion to arbitrary finite order. The even moment $\int z^{2k}H_s(z)\,dz=(2k)!s^{k+1}/(2^{k+1}(k+1)!)$ reproduces the coefficients from round 38. Thus these are classical derivatives for the non-exact construction. In particular
$$
[y^n]S_R=I_n(w)+B_n,\quad
B_n=\sum_{j=1}^{\lfloor(n-2)/2\rfloor}
\frac{(n+2j-2)!}{2^j j!\,n!}[u^{n-2j-2}](g w^j).
$$
A fixed complex disc $|y|\le |u|/4$ gives a bound $C(4/|u|)^n e^{-c'/u^2}$ for the non-diagonal coefficient integrand. This is integrable uniformly over the parameter family. Differentiating in $w$ replaces $H_{u^4w}$ by $u^4p_{u^4w}/2$, with an analogous majorant, proving uniform Lipschitz bounds on $I_n$ for each fixed $n$.

## Triangular realization

For $n=k+4$, the coefficient $B_n$ is $g(0)a_k/2$ plus a polynomial in the lower jet coefficients, where $a_k=[u^k]w$. Prescribe the target $[y^n]S_R$ for all $n\ge5$. Choose the parameter bounds $M_k$ successively to dominate the coordinate update on the previous box, using the uniform bound on $I_n$. Then use
$$
w_z(u)=c+\sum_{k\ge1}M_kz_k u^k\chi(u/\epsilon_k),\qquad |z_k|\le1 .
$$
Choose $\epsilon_k$ so the amplitude is at most $c2^{-k-2}$ and the $C^j$ norm is at most $2^{-k}$ for $j\le\lfloor k/2\rfloor$. This gives uniform convergence in every finite $C^j$ norm and continuity in the product topology. Finite Brouwer followed by diagonal compactness is valid. Alternatively, the amplitude bounds make the map $z\mapsto w_z$ uniformly Lipschitz from $\ell^\infty$ to $C^0$. Choosing $M_k$ additionally to dominate twice the coordinate update's Lipschitz constant, before setting $\epsilon_k$, gives a contraction on the closed $\ell^\infty$ cube.

With all target coefficients zero, the resulting $S_R$ does not equal its degree-four polynomial on either side near 0 if $\sup_{|u|\le R}v(u)<1$. If it did, on the positive side $D_R=g+P''$. Both sides extend real analytically throughout $y>0$ because only the source point $u=0$ has vanishing variance, away from every nonzero $y$. Decay at infinity forces $P''=0$, whereas $D_R(y)/g(y)\to0$ when the bounded core variances have supremum below 1. The negative side is identical. Thus the pending remainder is genuinely flat and nonzero on either side.

## Exterior finite-transform controls

Fix $b>0$ and distinct nonzero real parameters $t_1,\ldots,t_m$. Use a smooth transition to the constant exterior variance $b$ inside $R<|x|<2R$. Set all Borel support radii below $R/2$, so the transition uses the fixed quartic $c x^4$ where required. Choose exterior bump functions with disjoint supports fixed independently of $R$. To obtain an invertible Jacobian for the mean variance and the $m$ output moment generating values, append $t_0=0$ and choose the bump centers $x_j$ so that $[e^{t_i x_j}]_{i,j}$ is invertible. Such choices exist on either open exterior ray because the exponentials with distinct exponents are linearly independent there. The finite bump functions can be shrunk and normalized with respect to the Gaussian source so their response matrix is arbitrarily close to this evaluation matrix, with row multipliers $(t_i^2/2)e^{t_i^2b/2}$ for $i>0$.

The exterior response map in the bump parameters is independent of $w$. The discrepancy caused by the core and its transition is $O(R)$ uniformly in the Borel family, because the variance and tilt parameters are bounded. Choose $R$ small once so all discrepancies lie in one inverse-function neighborhood for the fixed exterior response. The resulting exterior coefficients are continuous and Lipschitz in $w$ in the uniform norm. They remain in a fixed small parameter ball and the exterior variances remain above $b/2$ on the bump supports.

The full output density $D$ has coefficient of order $n$ equal to
$$
(n+1)(n+2)[y^{n+2}]S_R+[y^n]g+[y^n]D_{\mathrm{ext}}.
$$
The exterior term is analytic near 0, with coefficients bounded uniformly over the exterior parameter ball and Lipschitz in its parameters. For $n\ge3$, use the preceding triangular update with $k=n-2$. The prescribed output coefficients can therefore be matched to the Gaussian target $p_{1+b}$ simultaneously with the variance constraint and all finite-transform constraints.

Finally, $B_2=0$ gives for the smooth density representative
$$
D(0)=g(0)+\int_{\mathbb R}g(u)p_{v(u)}(-u)\,du>p_{1+b}(0).
$$
At nonzero $y$ the full density is real analytic, as follows by complex Gaussian domination away from $u=0$ and the positive variance floor off a neighborhood of 0. If $D-p_{1+b}$ agreed with its degree-two Taylor polynomial near 0 on either side, real analytic continuation on that halfline would extend the polynomial there. The density difference is integrable, forcing the polynomial to vanish, contrary to the value at 0. The matching construction rules out sufficiency of finitely many output transform evaluations together with all the specified higher Taylor coefficients.
