# Integer difference witnesses for negative time traces

Let $\mu$ be a finite signed measure on $[0,\infty)$ and $R(q)=\int e^{-qt}\mu(dt)$. Use the forward difference $\Delta R(j)=R(j+1)-R(j)$. Then
$$
(-1)^k\Delta^k R(j)
=\sum_{\ell=0}^k(-1)^\ell\binom{k}{\ell}R(j+\ell)
=\int e^{-jt}(1-e^{-t})^k\,\mu(dt).
$$
These identities use the actual time-Laplace transform. Multiplying it by a nonconstant factor before taking differences changes the test.

## Explicit orders from a negative interval

Suppose $\mu|_{[a,b]}\le0$, $0\le a<b<\infty$, with positive negative mass in the inner subinterval specified below. Set
$$
x_-=e^{-b},\quad x_+=e^{-a},\quad q=(x_-+x_+)/2,\quad
d=(x_+-x_-)/2,\quad A=q(1-q),\quad r=d/16.
$$
Let $m_0>0$ be any lower bound for the negative mass where $x=e^{-t}\in[q-r,q+r]$, and $M$ any upper bound for the total positive mass. If the time density is at most $-m<0$ on $[a,b]$, one may take
$$
m_0=(m)\log\frac{q+r}{q-r}\ge 2mr/q.
$$
Choose an integer $n$ satisfying
$$
n\ge \max\left\{\frac8d,\ \frac{5A}{d^2}\log\frac{2M}{m_0}\right\},
$$
with the second constraint omitted if $M=0$. Put $j=\operatorname{round}(nq)$, $k=n-j$, and $p=j/n$. Then $j,k\ge1$ and $|p-q|\le r$. The normalized kernel is
$$
\frac{x^j(1-x)^k}{p^j(1-p)^k}
=\exp[-nD_{\rm Ber}(p\|x)] .
$$

Indeed $D_{\rm Ber}(p\|x)\le(p-x)^2/[x(1-x)]$ by $\log u\le u-1$. On $|x-q|\le r$, both $x/q$ and $(1-x)/(1-q)$ are at least $15/16$, while $|x-p|\le d/8$. Thus the divergence is at most $B_{\rm in}=4d^2/(225 A)$. On the whole interval $[x_-,x_+]$, $x(1-x)\le2A$ since one factor is at most its value at $q$ and the other is at most twice. As $\partial_x D_{\rm Ber}(p\|x)=(x-p)/[x(1-x)]$, its values on either boundary are at least $225d^2/(1024 A)$. Divergence grows farther outward. The boundary bound exceeds $B_{\rm in}+d^2/(5A)$. Dropping other negative contributions gives
$$
(-1)^k\Delta^kR(j)
\le-\tfrac12 m_0\,p^j(1-p)^k e^{-n B_{\rm in}}<0 .
$$
This is a finite explicit witness. For approximate samples with a uniform absolute error $\epsilon$, the sum incurs error at most $2^k\epsilon$. Direct integration with the positive kernel avoids that alternating-sum cancellation.

For positive step $h$, replace $x$ by $e^{-ht}$ and sample $R(q_0+jh)$; absorb $e^{-q_0t}$ into $\mu$. A known density bound $-m$ on $[a,b]$ supplies $m_0\ge2mr\,e^{-q_0b}/(hq)$ for $q_0\ge0$.

For a heat remainder tested against a fixed nonnegative $\phi\in C_b^2$, a single negative value $Q_t(\phi)\le-\delta$ supplies a negative time segment. The semigroup relation and positivity of newly absorbed sources give
$$
Q_{t+s}(\phi)
\le Q_t(\phi)+\|Q_t\|_{\rm TV}\|P_s\phi-\phi\|_\infty
\le Q_t(\phi)+s\|\phi''\|_\infty ,
$$
since $\|Q_t\|_{\rm TV}\le2$. Thus the value remains at most $-\delta/2$ for $0\le s\le\delta/(2\|\phi''\|_\infty)$ within the horizon. One may take $M\le H\|\phi\|_\infty$. For a pointwise time density at fixed spatial $y$, the free kernel bounds the positive part, giving $M\le\sqrt{2H/\pi}$.

## Two-layer center calibration

For the R002 two-layer data $s=1/2$, source variance $v=1/8$ and weight $2/3$, the signed time density at spatial zero is
$$
Q_t(0)=p_t(0)-\tfrac23\,1_{t\ge1/2}\,p_{t-3/8}(0),\qquad 0<t<1,
$$
and is zero for $t\ge1$. Its sign is strictly negative for $1/2\le t<27/40$. On $[1/2,3/5]$ the source/free ratio is at least $\sqrt{32/27}>13/12$ and $p_t(0)>1/2$, so the density is below $-1/24$. The positive mass on the whole horizon is at most $\sqrt{2/\pi}<1$.

For an elementary large-order witness, take the smaller negative window $x\in[11/20,3/5]$. Then $q=23/40$ and $d=1/40$. Take $n=17600$, so $j=10120$, $k=7480$ and $p=q$ exactly. Because there is no rounding, the inner radius may be enlarged to $d/8$: its divergence bound is $d^2/(49A)$ and the exterior bound is $d^2/(4A)$, still leaving a gap greater than $d^2/(5A)=1/1955$. The inner negative mass is at least $1/2208$, and the exponent $17600/1955>9$ exceeds $\log(4416)$.

A smaller rigorous witness is $(j,k)=(80,70)$. The companion exact script bounds
$$
\frac{\sqrt{2\pi}}{p^j(1-p)^k}
\int e^{-jt}(1-e^{-t})^k Q_t(0)\,dt
< -0.0043567.
$$
It uses 512 equal time cells, rational alternating-series bounds for $e^{-t}$, integer-square-root bounds for the densities, and unimodality of the normalized kernel. On $[0,1/16]$ the density singularity is integrated separately. The script checks the upper sum in exact rational arithmetic. Separate two-method numerical integration gives about $-0.00556522$ for the same normalized expression. The smaller pair $(59,53)$ also gives a negative numerical value, while the certified statement uses $(80,70)$.
