# Fixed-core residual and translated tilt

Let $g=p_1$, $F(t)=1$, and $\nu_t(dx)=e^{tx-t^2 f(x)/2}g(x)\,dx$. Fix $K=[-r,r]$, write $m(t)=\nu_t(K)$ and $p(t)=1-m(t)>0$, and put $E_t=\{|X-t|\le t/2\}$ for $t>0$. Let $W_t=X-tf(X)$. For all real $s$, exactly

$$
E_{\nu_t}[1_K e^{sW_t-s^2 f/2}]=m(t+s),\qquad
E_{\nu_t}[1_{K^c}e^{sW_t-s^2 f/2}]=p(t+s).
$$

This identity avoids an additive estimate on the change in the core mass. If the unfavorable side has denominator tending to infinity, its mass is $g(0)/t+o(1/t)$, by the substitution $y=t|x|$ and domination by $g(0)e^{-y}$. The bilateral core therefore matters at the $1/t$ scale.

## Relative translation

Suppose $f$ is positive outside the interior of $K$, each tail has lower limit at least $1$, and $\sup_{E_t}t(1-f)_+\to0$. If

$$
-\log p(t)=o(t),\qquad p(t+s)/p(t)\to1\quad(s\in\mathbb R\text{ fixed}),
$$

then the law of $X-t$ conditioned on $E_t$ tends to $N(0,1)$. If $f$ is uniformly continuous, this implies $f(t)\to1$.

On a fixed compact annulus outside $K$ the kernel is $O(e^{-c t^2})$ since $f$ has a positive minimum there. On both tails outside $E_t$, the bound $f\ge1-1/16$ and the tilted Gaussian representation give the same estimate. These claims hold with $t+s$ replacing $t$ while retaining $E_t$, since $|X-t|>t/2$ then implies $|X-(t+s)|\ge t/3$. The shifted identity yields

$$
E_{\nu_t}[1_{E_t}e^{sW_t-s^2f/2}]=p(t+s)+O_s(e^{-c' t^2}).
$$

In particular $\nu_t(E_t)\sim p(t)$. The bound $\nu_t(E_t\cap\{t(f-1)>\delta\})\le e^{-\delta t/2}$ remains negligible relative to $p(t)$, so $t(f-1)\to0$ conditionally. To remove the factor in the exponential transform, use

$$
E_{\nu_t}[1_{E_t}e^{2sW_t}]
\le e^{2s^2\|f\|_\infty}E_{\nu_t}[1_{E_t}e^{2sW_t-2s^2f}]
=O_s(p(t)).
$$

After conditioning, Cauchy–Schwarz applies with a uniformly bounded second moment, while $f\to1$ conditionally. The moment generating functions of $W_t$ converge to the standard Gaussian moment generating function. Conditional convergence of $X-t$ follows from $t(f-1)\to0$. Uniform continuity rules out fixed-width intervals on which $f$ remains above $1+\delta$, since their unnormalized probability is exponentially small relative to $p(t)$.

## Curvature bound for the residual ratio

Assume $f(\pm u)=u^2/(2D_\pm(\ell))$ near zero, $\ell=\log(1/u)$, and $D_+(\ell)\le C\ell$ for all large $\ell$. Exactness itself forces this upper envelope whenever $D_+(\ell+h)/D_+(\ell)\to1$ for each fixed additive shift $h$. In fact $\limsup D_+(\ell)/\ell\le1$. Otherwise choose $\ell_n\to\infty$ such that $D_n=D_+(\ell_n)\ge(1+\eta)\ell_n$, put $u_n=e^{-\ell_n}$ and $t_n=2D_n/u_n$. At $x=v u_n$ with $v$ in a short fixed interval about $1$, the exponent divided by $D_n$ tends pointwise to $2v-v^2$. The subset where it is at least $1-\varepsilon$ has asymptotically full length in the interval, by bounded convergence of its indicators. For $\varepsilon$ sufficiently small, the mass lower bound $c\exp[(1-\varepsilon)D_n-\ell_n]$ diverges, contradicting exactness.

Put $L=\log t$. Differentiating the compact integral in the parameter gives

$$
m''(t)=E_{\nu_t}[1_K(W_t^2-f)].
$$

On the favorable half write $y=tu$ and $B=y^2/(4D)$. Then $W_t=(y-2B)/t$ and $f=2B/t^2$. On $y\le bL$ and $B\le(2b+2)L$, bound the two derivative factors by $O(L^2/t^2)$ and use that the total mass is bounded by one. For $y\le bL$ and larger $B$, the exponential $e^{y-B-L}$ absorbs every power of $B$. For $y\ge bL$, the upper bound on $D$ gives $B\ge y^2/(4CL)$, and for fixed $b$ large the factor $e^{y-B/2}$ decays exponentially in $y$. The remaining $e^{-B/2}$ absorbs powers of $B$. On the unfavorable side $W_t=(-y-2B)/t$ and $e^{-y-B}$ absorbs powers of $y$ and $B$, giving a contribution $O(t^{-3})$ including the Jacobian $1/t$. Hence $|m''(t)|=O(L^2/t^2)$.

If the favorable side satisfies $D_+(\ell)=\ell+\tfrac12\log\ell+O(1)$, the bound improves to $O(L/t^2)$. On $y\in[L/2,bL]$ use $z=(y-2L)/(2\sqrt L)$ and the identity

$$
e^{y-y^2/(4D)-L}=e^{D-L}e^{-(y-2D)^2/(4D)}.
$$

Here $e^{D-L}\lesssim L^{-1/2}$ and $D=L+O(\log L)$, giving a Gaussian majorant in $z$. On this region $|W_t|\lesssim(\sqrt L|z|+\log L)/t$ and $f\lesssim L/t^2$. For $y<L/2$, the factor $e^{y-L}$ makes the mass exponentially small in $L$ and the derivative factors are already bounded as above. The estimate for large $y$ is unchanged. Both curvature estimates require only the given bounds on $D_\pm$.

Choose $K_t$ to bound $\sup_{u\in[t/2,3t/2]}|p''(u)|$, with $K_t=CL^2/t^2$ or $CL/t^2$ as applicable. By $0<p(t)\le1$, the step $a=\sqrt{2p(t)/K_t}$ lies in the allowed Taylor range for large $t$. Taylor's theorem in the direction opposite the derivative and nonnegativity give

$$
|p'(t)|\le\sqrt{2K_t p(t)}.
$$

For each fixed $s$, a second Taylor estimate gives

$$
|p(t+s)-p(t)|\le |s|\sqrt{2K_t p(t)}+K_t s^2/2.
$$

Thus $p(t)\gg L^2/t^2$ or $p(t)\gg L/t^2$, respectively, suffices for both relative hypotheses. At smaller residuals the relative theorem still applies when its ratio and subexponential bound can be verified directly.
