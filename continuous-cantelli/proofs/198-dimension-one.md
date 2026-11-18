# Dimension-one gauges and a logarithmic time-difference budget

Assume $X+\sqrt{v(X)}G$ is exactly Gaussian with variance $H>0$, where $X,G$ are independent standard normals and $v\ge0$ is bounded and continuous. Let $E$ be compact, contained in $\{v=0\}$, and let $\mu$ be a probability supported on $E$. Write $r=\operatorname{dist}(x,E)$, $h_b=\mu*p_b$ and $L_b=\log(1/\sqrt b)$.

## A gauge version allowing irregular growth

Suppose $\ell$ is continuous, nondecreasing and unbounded on $[L_0,\infty)$, and $\ell(L+1)\le D\ell(L)$ for a fixed $D<\infty$. Suppose that for small $u>0$
$$
c\,u\ell(\log(1/u))\le\mu(B(a,u))
\le C\,u\ell(\log(1/u))\qquad(a\in E).
$$
For each fixed $K$, some $\theta>0$ has the following property. If $r\le K\sqrt v$ throughout a collar, then
$$
\int_{0<r<\delta}\ell\!\left(\log(1/\sqrt{v(x)})\right)^{1+\theta}\,dx<\infty .
$$
One may shrink the collar to ensure $v<e^{-2L_0}$. The same bound holds on the good part $r\le K\sqrt v$ if the bad portion has finite integral of $\ell(\log(1/r))^{1+\theta}$.

Gaussian annuli and the upper gauge bound give
$$
|\partial_b^j h_b(x)|\le C_j b^{-j}\ell(L_b)e^{-r^2/(8b)}
\qquad(j=0,1)
$$
for small $b$. In detail, bound the integral against $\mu$ of a wider Gaussian by summing balls of radius $m\sqrt b$, using monotonicity $\ell(L_b-\log m)\le\ell(L_b)$. Arbitrary ball centers use the bound at a point of $E$ in the ball with twice its radius. The extracted exponential uses the minimal source distance. When $r\le K\sqrt b$, a ball of radius $\sqrt b$ about a nearest point instead gives $h_b(x)\ge c_K\ell(L_b)$. Fixed ages away from zero are bounded. The unit-ratio condition implies $\ell(L+y)\le D^{1+\lceil y\rceil}\ell(L)$ for $y\ge0$, so the extracted Gaussian also gives $\sup_b|\partial_b^j h_b(x)|\le C_j r^{-2j}\ell(\log(1/r))$ for small $r>0$.

Push the positive Stieltjes measure $d(\ell(u)^\theta)$, $u\ge L_0$, to time $t=e^{-2u}$, and call it $\Omega_\theta$. Take $\theta$ sufficiently small so $D^\theta<e^2$. Then $\int t\,\Omega_\theta(dt)<\infty$. Put
$$
F_\theta(b,x)=\int (h_b-h_{b+t})(x)\,\Omega_\theta(dt).
$$
For $t<b$, the absolute negative cost is bounded by a constant times
$$
\frac{\ell(L_b)}{b}\int_0^b t\,\Omega_\theta(dt)
\le
\ell(L_b)^{1+\theta}
\frac{D^\theta-1}{1-e^{-2}D^\theta}.
$$
To get this bound, split $u\ge L_b$ into unit intervals and bound each Stieltjes mass by $(D^\theta-1)D^{\theta j}\ell(L_b)^\theta$. For $t\ge b$, keep the entire $h_b$ term. It is at least $c_K\ell(L_b)[\ell(L_b)^\theta-\ell(L_0)^\theta]$ on the good set. The competing $h_{b+t}$ term is bounded by a constant times
$$
\int_{L_0}^{L_b}\ell(u)\,d(\ell(u)^\theta)
\le \frac{\theta}{1+\theta}\ell(L_b)^{1+\theta}.
$$
The Stieltjes chain rule holds because $\ell$ is continuous and monotone. Both negative coefficients tend to zero with $\theta$, so for a fixed sufficiently small $\theta>0$ and all small $b$,
$$
F_\theta(b,x)\ge c'_K\ell(L_b)^{1+\theta}\qquad(r\le K\sqrt b).
$$
Splitting instead at $t=r^2$ and using the uniform off-diagonal derivative and size estimates gives $|F_\theta(b,x)|\le C\ell(\log(1/r))^{1+\theta}$ for small $r>0$, uniformly in $b>0$. Outside a fixed collar the same type of bound is uniform.

Smooth the exact output by variance $s>0$. Tonelli for the underlying nonnegative heat mixtures and absolute Fubini for their difference yield
$$
\int g(x)F_\theta(v(x)+s,x)\,dx
=
\int\mu(da)\int[p_{H+s}(a)-p_{H+s+t}(a)]\,\Omega_\theta(dt).
$$
Absolute Fubini uses the lower age $s$ and $\int t\,\Omega_\theta(dt)<\infty$. The right side is bounded uniformly in $s$ because $H>0$. A maximal separated cover from the lower gauge gives $|E|=0$. The good part lower bound holds uniformly under $v\mapsto v+s$ after shrinking the collar, while the bad and the fixed exterior have an integrable negative majorant by the distance bound. Fatou and the positive minimum of $g$ on the compact collar prove the theorem.

The ball lower bound need only hold at a nearest point on the relevant $\sqrt{v+s}$ scales for each good source, uniformly for the smoothing used. The upper gauge and distance bound control the exceptional part. This gives a localized condition without two-sided regularity of the entire zero set.

## Consequences and a limitation

For a two-sided gauge, separated coverings give $|\{r<u\}|\asymp1/\ell(\log(1/u))$. The lower tube estimate uses the upper measure bound to force at least $c/[u\ell(\log(1/u))]$ covering balls, then their separated centers. For $\ell(L)=L^\beta$, $\beta>0$, the conclusion is a logarithmic gain from $\beta$ to $\beta(1+\theta)$ in the negative-variance moment. Thus a quadratic lower envelope $v\ge cr^2$ cannot coexist with $v\le\exp[-c'(\log(1/r))^p]$ throughout the collar when $p\ge1/(1+\theta)$. The proposed norm would dominate an integral of $(\log(1/r))^{p\beta(1+\theta)}$, which diverges by the tube estimate and layer cake. This extends the information obtained from direct Fatou on $h_v$, which gives only exponent $\beta$. Other examples covered are $\ell(L)=(\log L)^\beta$, $\ell(L)=\exp(L^\beta)$ for $0<\beta<1$, and non-smooth gauges with bounded unit ratio.

The usual fixed-order fractional kernel fails already for a null dimension-one Cantor set obeying a logarithmic gauge. Start with $[0,1]$, and at level $n$ use endpoint-retained equal children with lengths $r_n=2^{-n}(n+1)^{-\beta}$, $\beta>0$, assigning each cylinder mass $2^{-n}$. Balls of radius comparable to $r_n$ meet a bounded number of level-$n$ cylinders and contain one cylinder at a comparable finer level. Thus this is exactly a logarithmic-gauge example.

At the left endpoint put $x_n=-K r_n$, $b_n=r_n^2$ for $K>0$. The rescaled measures $\nu_n(B)=2^n\mu(r_n B)$ converge locally weakly to Lebesgue measure on $[0,\infty)$. To see this, fix an ancestor depth offset $m$: the leftmost level-$(n-m)$ cylinder has rescaled length tending to $2^m$ and normalized mass $2^m$. On that cylinder, each fixed number of subsequent splits tends to the ordinary dyadic geometry, with equal descendant masses. Then send the number of splits to infinity and choose $2^m$ to cover any fixed test window. The uniform gauge gives a uniform bound on the $\nu_n$ mass of each unit interval, which controls the Gaussian tails. Hence for fixed $u\ge0$,
$$
\frac{h_{b_n(1+u)}(x_n)}{(n+1)^\beta}
\to \bar\Phi(K/\sqrt{1+u}).
$$
The uniform Gaussian size and time derivative bounds dominate the normalized difference by $C\min(u,1)$. For every $0<\kappa<1$, dominated convergence therefore gives
$$
\frac{b_n^\kappa}{(n+1)^\beta}
\int_0^\infty[h_{b_n}-h_{b_n+t}](x_n)t^{-1-\kappa}\,dt
\to \int_0^\infty
[\bar\Phi(K)-\bar\Phi(K/\sqrt{1+u})]u^{-1-\kappa}\,du<0.
$$
Thus dimension one, even with an accurate local gauge, does not support the positive fixed-fractional-order step of round 193. The Stieltjes clock uses the unbounded gauge in a different manner. For a purely dimensional assumption with no usable common gauge or lower heat bound, this argument supplies no uniform positive budget.
