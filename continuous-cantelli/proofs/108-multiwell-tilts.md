# Spatial entropy, occupancy covariance and critical tilt shifts

Assume the exact bounded normalized graph
$$
P_t(dx)=e^{tx-t^2f(x)/2}\gamma(dx),\qquad P_t(\mathbb R)=1,\qquad 0\le f\le M,
$$
and set $W_t=X-tf(X)$ and $a(t)=E_t f$. Boundedness permits tilt-parameter differentiation without derivatives of $f$. The exact shifted likelihood and first two derivatives give
$$
\log\frac{dP_s}{dP_t}=(s-t)W_t-\tfrac12(s-t)^2 f,\quad
E_tW_t=0,\quad \operatorname{Var}_t(W_t)=a(t),\quad 0<a(t)\le1.
$$
The upper bound is $D(P_t\|N(t,1))=t^2(1-a(t))/2\ge0$, with $a(0)=1$. Direct covariance identities are
$$
a'(t)=\operatorname{Cov}_t(f,W_t),\qquad
2t\operatorname{Cov}_t(X,f)=\operatorname{Var}_t(X)+t^2\operatorname{Var}_t(f)-a(t).
$$

## Finite well tests

Take a finite measurable partition $B_i$ of the source into wells and their complements, and write $p_i=P_t(B_i)>0$. At the indicated tilt let
$$
x_i=E_t[X\mid B_i],\quad a_i=E_t[f\mid B_i],\quad
w_i=x_i-ta_i,\quad v_i=\operatorname{Var}_t(W_t\mid B_i).
$$
Then
$$
\sum_i p_i w_i=0,\qquad \sum_i p_i(w_i^2+v_i)=a(t),\qquad
\sum_i p_i e^{s w_i-s^2a_i/2}\le1\quad(s\in\mathbb R).
$$
The last inequality is conditional Jensen in the exact shifted likelihood. In particular, for $a_i>0$, $p_i\le e^{-w_i^2/(2a_i)}$. If $X\in[L_i,R_i]$ and $f\in[\ell_i,h_i]$ on a well, then at $t>0$ its score mean lies in $J_i=[L_i-th_i,R_i-t\ell_i]$. The covariance budget gives the test
$$
\sum_{i<j}p_i p_j\,\mathrm{dist}(J_i,J_j)^2\le a(t).
$$
Only pairs with useful score intervals need be retained. This distinguishes centers separated in space from centers whose different depths offset that separation.

For fixed cells the occupancy derivatives are
$$
(\log p_i)'=w_i,\qquad
(\log p_i)''=v_i-a_i,\qquad
\sum_i (p_i')^2/p_i\le a(t).
$$
The integrated form is a Hellinger speed bound
$$
\|\sqrt{p(t)}-\sqrt{p(s)}\|_2\le \frac12\left|\int_s^t\sqrt{a(u)}\,du\right|.
$$
Moving wells cannot be differentiated with these formulas. At any two tilts, the conditional likelihood identity gives
$$
D(P_t(\cdot\mid B_i)\|P_s(\cdot\mid B_i))
=(t-s)w_i+\tfrac12(t-s)^2 a_i+\log\frac{p_i(s)}{p_i(t)}\ge0.
$$
Thus each fixed well's mass at $s$ has a computable lower envelope from its score and depth at $t$. The categorical KL also satisfies $D(p(t)\|p(s))\le(t-s)^2 a(t)/2$ by data processing of the exact pairwise entropy.

Write $s_i^2=\operatorname{Var}_t(X\mid B_i)$, $\sigma_i^2=\operatorname{Var}_t(f\mid B_i)$, and $H(p)=-\sum p_i\log p_i$. Entropy decomposes over disjoint wells, so the exact spatial entropy and maximal entropy at a given support length or variance give
$$
\frac12\log(2\pi)+\frac12\operatorname{Var}_t(X)-\frac12t^2 a(t)(1-a(t))
\le H(p)+\sum_i p_i\min\{\log|B_i|,\tfrac12\log(2\pi e\,s_i^2)\}.
$$
An infinite-length cell uses only its variance bound. In applying this estimate, $s_i\le\sqrt{v_i}+|t|\sigma_i$, and $\sum p_i v_i=a-\sum p_iw_i^2$. Thus a separation of well score means reduces the variance available for the within-well spatial entropy. Directly on a single cell, its conditional density gives
$$
\log p_i\le \log |B_i|-\tfrac12\log(2\pi)
+t x_i-\tfrac12t^2a_i-\tfrac12 E_t[X^2\mid B_i]
$$
when its spatial length is finite.

## Critical spatial entropy budget

At a positive large tilt set $L=\log t$, $A=t^2 a(t)$, $Z=tX$, $D=tW_t$, and $U=(Z-2L)/\sqrt A$. Assume $A\le2L$. The entropy-floor estimates of round 85 give $A\sim2L$ and
$$
E_t|U+D/\sqrt A|=O(\log A/\sqrt A),\qquad
h(U)=\frac12\log(2\pi e)+\frac12(E_tX^2-\delta),\quad
\delta=A+\log A-2L+1.
$$
For clarity the first estimate uses $p_Z(z)=g(z/t)e^{(z-2L+D)/2}$. On $z-2L=-w<0$ split according as $D\ge w/2$ or not. The first contribution to $E_t(Z-2L)^-$ is at most $2\sqrt A$ and the other has an integrable exponential bound. Since $E_t(Z-2L)\le0$, a Laplace comparison centered at $2L$ with scale $1+E_t|Z-2L|$ has cross entropy $O(\log A)$. The elementary bound $\int_{p<q}p\log(q/p)\le1/e$ then gives $E_t(Z-2L+D)^-=O(\log A)$ and, by its nonpositive mean, the stated $L^1$ relation.

Let $Q_K$ have density proportional to $e^{-\psi_K/2}$, where $\psi_K(x)=x^2$ on $[-K,K]$ and equals $2K|x|-K^2$ outside. Take $K=\sqrt{\log A}$. Since $E_tD^2/A=1$, $\psi_K\le x^2$ and $\psi_K$ is $2K$-Lipschitz,
$$
D(\mathcal L_t(U)\|Q_K)
\le\frac12(\delta-E_tX^2)
+O\!\left(\frac{(\log L)^{3/2}}{\sqrt L}\right).
$$
Indeed the normalizing constant is $\sqrt{2\pi}+O(e^{-K^2/2}/K)$ and $E_t\psi_K(U)\le1+O(K\log A/\sqrt A)$. This gives a categorical KL upper bound for any partition on the $U$ scale, including moving finite well windows. In particular if a union of wells has law mass $p$ and comparison mass $q$, then
$$
q\ge \exp\!\left[-\frac{\delta/2+H_2(p)+o(1)}{p}\right].
$$
The reference density is at most $(2\pi)^{-1/2}$ and converges in total variation to the standard Gaussian. Bounded $\delta$ therefore forbids positive mass in wells of vanishing normalized total length. If $p\to1$ and the reference union masses tend to $q_*$, then $\liminf\delta\ge-2\log q_*$. Without bounded $\delta$, its upper bound $\log(2L)+1$ still gives a normalized width lower bound $L^{-1/2+o(1)}$ for unions whose mass tends to one.

## A critical neighboring-tilt escape test

Let $t_n\to\infty$ with $A_n\le2\log t_n$ and let $U_n$ converge in distribution to $\mu$ along a subsequence. Tightness follows from $E D^2/A=1$ and the preceding $L^1$ relation. For fixed real $\lambda$ let $u_n=t_n+\lambda t_n/\sqrt{A_n}$ and $L_n^*=dP_{u_n}/dP_{t_n}$. Under $P_{t_n}$,
$$
\log L_n^*=\lambda D/\sqrt{A_n}-\tfrac12\lambda^2t_n^2 f/A_n
\ \Longrightarrow\ -\lambda U-\lambda^2/2 .
$$
Here $t_n^2 f/A_n\to1$ in probability because $t_n^2 f=2\log t_n+\sqrt{A_n}U_n-D$. Fatou gives
$$
\Delta_\mu(\lambda)=1-e^{-\lambda^2/2}\int e^{-\lambda v}\mu(dv)\ge0
$$
for every real $\lambda$. If $A(u_n)/A_n=O(1)$ at this $\lambda$, reverse relative entropy is bounded:
$$
D(P_{u_n}\|P_{t_n})=
\tfrac12\lambda^2(t_n/u_n)^2 A(u_n)/A_n.
$$
It gives uniform integrability of the likelihood ratios. In particular $\Delta_\mu(\lambda)=0$. If this energy ratio is bounded for each $\lambda$ in any interval adjacent to zero, the Gaussian moment generating function is forced on an interval, hence $\mu=N(0,1)$. Subgaussianity ensures all exponential moments and permits analytic uniqueness even for a one-sided interval.

If a candidate well limit $\mu$ is nonnormal, some arbitrarily small $\lambda\ne0$ has $\Delta_\mu(\lambda)>0$. Then along the same subsequence the neighboring energy ratio $A(u_n)/A_n$ tends to infinity. There is also a finite test. Under $P_{u_n}$ the log-likelihood is $\lambda(t_n/u_n)D_{u_n}/\sqrt{A_n}+\tfrac12\lambda^2t_n^2 f/A_n$, so Chebyshev and Markov yield
$$
P_{u_n}(\log L_n^*>K)\le
\lambda^2(t_n/u_n)^2[A(u_n)/A_n](4/K^2+1/K).
$$
For each continuity threshold $K$, the limit inferior of the left side is at least $\Delta_\mu(\lambda)$ by convergence of the truncated likelihood expectations. This quantifies how much extra tilted energy a nonnormal multiwell shape necessitates at neighboring tilts. It leaves open whether such energy transitions can occur for an exact continuous graph.
