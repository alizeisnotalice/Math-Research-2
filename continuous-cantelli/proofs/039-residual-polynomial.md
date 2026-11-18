# Relative translations of small origin residuals

Let $f$ be a bounded continuous exact graph. Write $M=\sup f$, $\gamma_t=N(t,1)$, $\nu_t=e^{-t^2(f-1)/2}\gamma_t$ and $W_t=X-tf(X)$. Fix a compact origin core $K=[-r,r]$ and write $m(t)=\nu_t(K)$ and $p(t)=1-m(t)$. These are bilateral core and residual masses. In particular $p(t)>0$, and for every real $s$,

$$
E_{\nu_t}[1_{K^c}e^{sW_t-s^2 f/2}]=p(t+s).
$$

The equality follows by expanding $(t+s)^2$. If the left denominator tends to infinity, the left core has mass $g(0)/t+o(1/t)$ when $t\to+\infty$. This follows from dominated convergence after substituting $y=t|x|$. A one-sided core thus gives a different residual at the scales studied below.

## All translate derivatives under a one-sided envelope

Write $L=\log t$ and $f(u)=u^2/[2D(\ell)]$ for $u>0$ small, with $\ell=\log(1/u)$ and $D>0$. Fix any positive integer $N$. Under $D(\ell)\le C\ell$, one has

$$
|m^{(N)}(t)|=O_N((L/t)^N).
$$

Under the narrower upper bound $D(\ell)\le\ell+\tfrac12\log\ell+C$, one has

$$
|m^{(N)}(t)|=O_N((\sqrt L/t)^N).
$$

Neither estimate requires derivatives or a matching lower bound on $D$. In fact, if $D(\ell+h)/D(\ell)\to1$ for every fixed $h$, then exactness gives $\limsup D(\ell)/\ell\le1$. Otherwise, choose $D_n=D(\ell_n)\ge(1+\eta)\ell_n$ and test $t_n=2D_n e^{\ell_n}$. At $x=v e^{-\ell_n}$ the exponent divided by $D_n$ tends pointwise to $2v-v^2$. On a sufficiently short fixed interval around $v=1$, bounded convergence of indicators yields a subset of fixed positive length where the ratio exceeds $1-\varepsilon$. If $(1-\varepsilon)(1+\eta)>1$, the resulting lower bound on total tilted mass diverges.

For the derivative estimates, parameter differentiation yields a sum of terms $W_t^{N-2j}f^j$ with coefficients $N!(-1/2)^j/[j!(N-2j)!]$. On the unfavorable half of $K$ set $y=t|x|$ and $B=t^2 f(x)/2$. Its exponential factor is $e^{-y-B}$ and $W_t=(-y-2B)/t$. Its contribution to each absolute derivative term is $O_N(t^{-N-1})$, without any extra hypothesis on the unfavorable germ.

On the favorable side $y=tx$, $B=y^2/(4D)$, $W_t=(y-2B)/t$. For the broad envelope split at $y=bL$ and $B=(2b+2)L$, with $b$ sufficiently large depending on $C$. On the region below both thresholds, all derivative factors are $O_N((L/t)^N)$ and the total tilted mass is at most 1. If $y\le bL$ but $B>(2b+2)L$, the exponent $y-B-L$ suppresses every polynomial factor in $B$. For $y>bL$, the upper envelope gives $B\ge y^2/(4CL)$, so $y-B/2\le-cy$ for sufficiently large $b$. The remaining exponential factor $e^{-B/2}$ controls the polynomial factors.

For the narrow upper bound, the same estimates handle $y>bL$, and for $y<L/2$ the factor $e^{y-L}$ gives exponential decay in $L$, with $e^{-B}$ controlling powers of $B$. On $y\in[L/2,bL]$ one has $D\le L-\tfrac12\log L+C'$. If $D<L/4$, then $B\le(2b+2)L$ permits the broad derivative bound and the factor $e^{D-L}\le e^{-3L/4}$ in the square completion
$$
e^{y-B-L}=e^{D-L}e^{-w^2},\qquad w=(y-2D)/(2\sqrt D).
$$
If $B>(2b+2)L$, use the original exponent as before. If $D\ge L/4$, set $q_*=\tfrac12\log L-C'$, $Q=L-D-q_*\ge0$ and $a=y-2L+2q_*$. Then $w=(a+2Q)/(2\sqrt D)$, $Q\le L$ for large $L$, and
$$
Q+w^2/2\ge c a^2/L.
$$
This follows from $a^2\le2(a+2Q)^2+8Q^2$, $D\le L$, and $Q^2\le LQ$. Also $|W_t|\le C\sqrt L\,|w|/t$ and $f\le CL/t^2$. The derivative integrand is therefore bounded by a constant times $(\sqrt L/t)^N L^{-1/2}e^{-c a^2/L}$, after absorbing powers of $w$ into $e^{-w^2/2}$. Integration in $y$ proves the narrower estimate.

For $N=2$ this verifies the previously claimed $O(L/t^2)$ curvature from a one-sided upper bound. The Glaeser-style inequality from nonnegative $p$ and a local bound $|p''|\le K_t$ yields $|p'(t)|\le \sqrt{2K_tp(t)}$ and relative freezing when $p(t)/K_t\to\infty$. At smaller scales the next argument replaces freezing by a classification of subsequential laws.

## Polynomial residual floors and positive interval mass

Suppose each tail has $\liminf f\ge1$, $f$ is positive away from the origin, and $\sup_{E_t}t(1-f)_+\to0$ for $E_t=\{|X-t|\le t/2\}$. Consider any sequence $t\to\infty$ along which $p(t)\ge c t^{-A}$ for fixed $c>0$ and $A<\infty$. It suffices here to have the broad derivative estimates above. For each fixed $s$, the tilted mass on $K^c\setminus E_t$ with parameter $t+s$ is $O_s(e^{-c' t^2})$. On compact annuli this follows from positivity of $f$. On tails, choose an eventual $f\ge1-1/16$ and use the Gaussian shift, since $|X-(t+s)|\ge t/3$ for $t$ large.

The positive residual admits a separate curvature identity:
$$
(\log p)''(t)=\operatorname{Var}_{\nu_t(\cdot\mid K^c)}W_t
-E_{\nu_t(\cdot\mid K^c)}f\ge-M.
$$
Consequently $q_t(s)=p(t+s)/p(t)$ obeys $q_t(s)^2\le e^{Ms^2}q_t(2s)$.

Choose an integer $N>A$. Taylor's theorem, the absolute derivative estimates, and $p^{(N)}=-m^{(N)}$ show that $q_t(s)$ is its degree-$N-1$ Taylor polynomial at 0 plus $o(1)$, locally uniformly in $s$. Its coefficients are bounded along the sequence. If their maximum $A_t$ diverged along a subsequence, the polynomials divided by $A_t$ would tend to a nonzero polynomial $Q$ with $Q(0)=0$ and $Q(s)\ge0$. Taking $s$ with $Q(s)>0$ makes $q_t(s)^2$ of order $A_t^2$, while $e^{Ms^2}q_t(2s)=O_s(A_t)$, a contradiction.

The bounded coefficients give a uniform upper bound on $q_t(s)$ for $s$ in any fixed compact interval. Log semiconvexity also gives $q_t(s)q_t(-s)\ge e^{-Ms^2}$, so the same compact interval has a uniform positive lower bound. In particular $q_t(s_t)\to1$ when $s_t\to0$. Every subsequence admits a further subsequence whose translate ratio tends locally uniformly to a polynomial $P$ with $P(0)=1$. The tilted transform on $E_t$ satisfies
$$
E_{\nu_t}[1_{E_t}e^{sW_t-s^2f/2}]=p(t+s)+O_s(e^{-c't^2}).
$$
In particular $\nu_t(E_t)\sim p(t)$. Under conditioning on $E_t$, $t(f-1)\to0$ in probability. The negative part is uniformly small by hypothesis. The positive event $t(f-1)>\delta$ has unconditioned probability at most $e^{-\delta t/2}$, negligible relative to the polynomial floor. Since $f\le M$, the shift formula at $2s$ also gives bounded conditional moments of $e^{2sW_t}$. Cauchy–Schwarz then gives convergence of the conditional moment generating function of $W_t$ to $e^{s^2/2}P(s)$. The conditional law of $X-t$ has the same weak limit.

Writing $P(s)=\sum_j b_js^j$ and using the probabilists' Hermite polynomials,
$$
e^{s^2/2}P(s)=\int_{\mathbb R} e^{sx}g(x)\sum_j b_j\operatorname{He}_j(x)\,dx.
$$
The limit is a probability law, so its density is nonnegative. The multiplier is a nonzero polynomial and thus its zeros are isolated. Every nonempty open interval has positive limiting mass. In particular, a fixed-height interval about $t$ contradicts the polynomial residual floor, since its unconditioned tilted mass is exponentially small.

Thus whenever $f\ge1+\delta$ on fixed-radius intervals about $t_n\to\infty$, it follows that $p(t_n)=o(t_n^{-A})$ for every fixed $A$. A uniform polynomial lower bound on $p(t)$ therefore gives tail return under uniform continuity. The same proof also works under the weaker, sequence-dependent condition $p(t)\gg (L/t)^N$ for some fixed $N$ in the broad case, or $p(t)\gg(\sqrt L/t)^N$ under the narrow upper envelope. The result does not assert translate ratios tend to 1 at these finer scales. For example, at a quadratic limiting polynomial $P(s)=1+a s+b s^2/2$, the weak limit can instead have density $g(x)[1+a x+b(x^2-1)/2]$, which is nonnegative precisely when $0\le b\le2$ and $a^2\le2b-b^2$.
