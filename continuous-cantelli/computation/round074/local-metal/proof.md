# Cantor-gap Fourier identity and ternary subquadratic resonance

Identify the circle with $[0,1)$. The two outer similarities $x\mapsto\lambda x$ and $x\mapsto1-\lambda+\lambda x$, $0<\lambda<1/2$, define $C_\lambda$. Write $\delta(x)=\operatorname{dist}_{\mathrm{circle}}(x,C_\lambda)$. With uniform initial point $x$, a centered Gaussian increment of variance $d\delta(x)^p$ gives the Fourier coefficient
$$
F_n=\int_0^1 e^{2\pi i n x}e^{-2\pi^2 n^2 d\,\delta(x)^p}\,dx,\qquad n\ge1.
$$

## Exact compression by level

The $2^{j-1}$ gaps removed at level $j$ have length $l_j=(1-2\lambda)\lambda^{j-1}$. Their centers have form
$$
\frac{\lambda^{j-1}}2+(1-\lambda)\sum_{k=0}^{j-2}\epsilon_k\lambda^k,\qquad\epsilon_k\in\{0,1\}.
$$
Their center phase sum is
$$
(-1)^n2^{j-1}P_j,\qquad P_j=\prod_{k=0}^{j-2}\cos(\pi n(1-\lambda)\lambda^k).
$$
The phases factor bit by bit. The common phase is $e^{\pi i n}$, since the center expression with each digit replaced by $1/2$ is $1/2$. On each gap, distance is measured to its nearer endpoint. Subtract the coefficient at zero noise, which is zero for $n\ge1$, then integrate both half-gaps in the normalized distance $t\in[0,1]$. For $s_j=2\pi^2 n^2d(l_j/2)^p$, $\omega_j=\pi n l_j$, and $m_j=2^{j-1}l_j$, one obtains
$$
F_n=(-1)^n\sum_{j\ge1}m_j P_j
 \int_0^1\cos(\omega_j(1-t))(e^{-s_jt^p}-1)\,dt.
$$
The convergence follows from $\sum m_j=1$. This also proves the imaginary part vanishes. The scan calls a level coarse when $n l_j\ge1$ or $s_j\ge1$, and fine otherwise.

The unintegrated series tail after $J$ is bounded in absolute value by
$$
\frac{m_1s_1(2\lambda^{p+1})^J}{(p+1)(1-2\lambda^{p+1})}.
$$
Use $|P_j|\le1$, $|1-e^{-u}|\le u$, and the geometric progression $m_js_j$. In a small level integral, replacing the integrand integral by $-s_j/(p+1)$ costs at most
$$
\frac{s_j^2}{2(2p+1)}
+\frac{s_j\omega_j^2}{(p+1)(p+2)(p+3)}.
$$
The first term comes from $0\le e^{-u}-1+u\le u^2/2$. The second follows from $1-\cos z\le z^2/2$ and integrating $t^p(1-t)^2$. The declared local usage threshold is applied to this error times $|m_jP_j|$.

For the uniform cap at $J=96$, the tail bound decreases in $p$, since $l_1/2,\lambda<1$ and its denominator increases. It suffices to set $p=1/2$. The numerator's $\lambda$-dependent factor $(1-2\lambda)^{3/2}\lambda^{144}$ increases for $\lambda\le0.461$, since its log derivative has the sign of $144-291\lambda$; the denominator decreases. Use $\lambda=0.461$ conservatively, including binary32 endpoints. The exact Fraction check bounds this tail below $10^{-12}$. The per-level small-integral approximation contributes at most $96(2\times10^{-14})$ in exact arithmetic. Discarding the heat term past exponent 42, replacing expm1 by its degree-six polynomial for arguments in $[-0.06,0]$, and using degree-four sinc for $|\omega|<0.01$, together keep the analytically omitted total below $4\times10^{-12}$. This does not bound binary32 or quadrature error. Independent numerical checks address those effects.

## Ternary subquadratic resonance

Let $\lambda=1/3$, $p<2$, and $n=3^k+r$ with $r$ a fixed integer and $k\to\infty$. Put $s=(2\pi^2 n^2d)^{-1/p}$, $w_n(x)=e^{-(\delta(x)/s)^p}$, and $M_n=\int_0^1w_n(x)\,dx$. Write $\alpha=\log_3 2$ and let $\mu_C$ denote the uniform Cantor probability. Then
$$
\frac{F_{3^k+r}}{M_{3^k+r}}\longrightarrow\widehat\mu_C(1)\widehat\mu_C(r).
$$

To prove it, take $L=\lfloor\gamma k\rfloor$ with $1<\gamma<2/p$. The union of the surviving depth-$L$ intervals has total heat mass split equally among its $2^L$ intervals. Within each interval the distance to $C$ is determined by that interval's scaled Cantor copy, because both interval endpoints already belong to $C$. Outside this union, the $2^L-1$ deleted gaps have combined heat mass $O(2^L s)$ by integrating the positive kernel out from both endpoints. Meanwhile, $M_n\gtrsim s^{1-\alpha}$: choose depth $J$ with $3^{-J}\le s<3^{-J+1}$; its surviving intervals have mass measure $(2/3)^J$ and the integrand is at least $e^{-1}$ there. Thus the outside heat mass fraction is $O(2^Ls^\alpha)=o(1)$.

The phase varies by at most $2\pi n3^{-L}=o(1)$ on each surviving interval. The average of its value at the left endpoints approximates both $F_n/M_n$ and $\widehat\mu_C(n)$, since $\mu_C$ gives equal mass to the same intervals. Represent a Cantor point as $B_k+3^{-k}X'$, with $3^k B_k$ an integer and $X'$ an independent Cantor point. Then
$$
\widehat\mu_C(3^k+r)=\mathbb E e^{2\pi i r B_k}\cdot\widehat\mu_C(1+r3^{-k})
    \longrightarrow\widehat\mu_C(r)\widehat\mu_C(1).
$$
The factor $\widehat\mu_C(1)$ is positive, since the digit product is $(-1)\cos(2\pi/3)\prod_{j\ge2}\cos(2\pi/3^j)>0$. The infinite positive product converges to a positive value by summability of $1-\cos(2\pi/3^j)$. This gives final positive signs for $r=0,\pm1$. No subquadratic resonance conclusion for general $\lambda$ is asserted.

## Quadratic and superquadratic obstruction on the circle

In the quadratic case $p=2$, choose any removed gap and one endpoint $a$. For $a+u$ inside that gap near the endpoint the variance is exactly $du^2$. In the periodically extended Gaussian output at $y=a$ after adding a variance $\varepsilon>0$, this half-gap supplies a lower bound proportional to
$$
\int_{\sqrt{\varepsilon}}^{l/2}\frac{1}{u}\,du
$$
since on that interval $\varepsilon+du^2\le(1+d)u^2$ and $\exp[-u^2/(2(\varepsilon+du^2))]\ge e^{-1/(2d)}$. This diverges. Uniform output would retain constant density after the extra heat smoothing, a contradiction.

When $p>2$, any $a\in C_\lambda$ has $v(a+u)\le d|u|^p=o(u^2)$. After adding $\varepsilon$, the arc $|u|\le K\sqrt\varepsilon$ contributes density tending to $\mathbb P(|Z|\le K)$, where the mass convention is uniform density one on the unit circle. This follows by changing variables to $u/\sqrt\varepsilon$ in the centered Gaussian kernel term. A disjoint interval inside a gap supplies a fixed positive density at $a$, uniformly as $\varepsilon\downarrow0$. Take $K$ large and then $\varepsilon$ small; the smoothed output density at $a$ exceeds one. Thus the standard ternary profile cannot yield a uniform circle output for any $p>0$, combining this argument, the quadratic argument, and the preceding subquadratic resonance. For other ratios the conclusion here covers $p\ge2$.
