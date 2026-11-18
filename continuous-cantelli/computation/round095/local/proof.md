# The reflected-wall semigroup and the early-time coefficient

Let $p_t$ be the centered normal density of variance $t$. Fix $\tau,h>0$, $b\ge0$, and $d>0$. The reflected-wall cap is
$$
q_b(x)=\mathbf1_{x<b}[p_\tau(x)-p_\tau(2b-x)].
$$
It is nonnegative for $x<b$. This is the cap for an old signed remainder supplied by the reflection comparison, not an assumed positive evolution of a general signed remainder.

Measure the exterior point at $y=b+d$. Put
$$
\sigma=\sqrt{\frac{\tau h}{\tau+h}},\quad
k=\frac{d\sigma}{h},\quad v=\frac{b\sigma}{\tau}.
$$
Substitute $x=b-\sigma w$ in the positive heat integral. The two reflected exponentials differ by $2\sinh(vw)$, so
$$
(P_hq_b)(b+d)=2v\sigma\,p_h(d)p_\tau(b)\,R(k,v),\qquad
R(k,v)=\int_0^\infty w e^{-kw-w^2/2}
                    \frac{\sinh(vw)}{vw}\,dw.
$$
Here $\sinh(z)/z=1$ at zero by continuity. When $b=0$, the heat term vanishes and $R(k,0)$ is the normalized small-wall limit.

## Exact derivative and two-sided bounds

Define
$$
M(z)=\int_0^\infty e^{-zw-w^2/2}\,dw,\qquad
L(z)=-M'(z)=\int_0^\infty we^{-zw-w^2/2}\,dw .
$$
The Gaussian identity is $M(z)=\sqrt{2\pi}e^{z^2/2}\bar\Phi(z)$. For $v>0$,
$$
R(k,v)=\frac{M(k-v)-M(k+v)}{2v}
      =\frac1{2v}\int_{k-v}^{k+v}L(z)\,dz.
$$
The function $L$ is strictly convex by twice differentiating the positive integral. It follows that
$$
L(k)\le R(k,v)\le\tfrac12[L(k-v)+L(k+v)].
$$
Integration by parts gives $L(k)=1-kM(k)$. This supplies the uniform small-wall derivative check without subtracting close Gaussian terms in Float32.

For $k>v$, put $a=k-v$ and $c=k+v$. Dropping the factor $e^{-w^2/2}$ gives an upper bound $1/(ac)$. The difference from it is bounded by replacing $1-e^{-w^2/2}$ with $w^2/2$. Evaluate the positive exponential integrals to obtain
$$
0\le 1-acR(k,v)\le \frac1{a^2}+\frac1{ac}+\frac1{c^2}.
$$
The bound extends by continuity to $v=0$. It proves the fixed-wall expansion: at fixed $b,\tau,d>0$, as $h\downarrow0$, $k\to\infty$, $v/k\to0$, and $k^2R\to1$. Consequently
$$
(P_hq_b)(b+d)\sim
\frac{2b\,p_\tau(b)}{\tau\sqrt{2\pi}\,d^2}
h^{3/2}\exp[-d^2/(2h)].
$$
The expansion also gives an explicit error criterion when parameters vary, provided $k-v\to\infty$ and $v/k\to0$. On the scan's far weak-wall box $k\ge16$, $v/k\le0.1$, the same rational estimates give the uniform finite bound $|k^2R-1|<0.013$.

For bounded $v$ as $k\to\infty$ one can expand the integrand at $w=0$. Uniformly for bounded $v$,
$$
e^{-w^2/2}\frac{\sinh(vw)}{vw}
=1+(v^2/6-1/2)w^2+(v^4/120-v^2/12+1/8)w^4+O(w^6)
$$
on $0\le w\le1$. The contribution of $w>1$ is exponentially small after using $\sinh(vw)/(vw)\le e^{vw}$; the same holds for the integrated polynomial tails. Integrating against $w e^{-kw}$ gives
$$
k^2R=1+\frac{v^2-3}{k^2}+\frac{v^4-10v^2+15}{k^4}+O(k^{-6}).
$$
At fixed $\tau,b,d$, substitution of the exact expressions for $k$ and $v$ yields
$$
k^2R
=1-\frac{3h}{d^2}
+h^2\left(\frac{b^2}{\tau^2d^2}-\frac3{\tau d^2}+\frac{15}{d^4}\right)
+O(h^3).
$$
These are coefficients of the reflected Gaussian cap. Applying the cap to an old signed remainder bounds its exterior heat contribution.

## Transition and a finite counterexample to a uniform fixed-wall factor

For fixed $a=k-v$, with $v\to\infty$,
$$
0<M(a)-2vR(k,v)=M(k+v)\le1/(k+v).
$$
In this transition, $k^2R\sim\tfrac12 k M(a)$, even for fixed positive $a$. Thus the fixed-wall limiting factor one does not give a uniform upper bound over all $k>v$. An explicit real-arithmetic example is $k=64,\ v=63$. On $w\in[0,1]$, $e^{-w-w^2/2}\ge e^{-3/2}>1/5$, where the last comparison follows from the Taylor polynomial of degree five. Since $M(127)\le1/127$, one has $64^2R(64,63)>6$. The rational arithmetic checks these constants; higher-precision integration records the stronger numerical size. This does not contradict the fixed-geometry expansion used for old-wall suppression.

If the exterior distance shrinks with $h$ so that $k\to k_0\in(0,\infty)$ and $v\to0$, then $k^2R\to k_0^2L(k_0)\in(0,1)$. This is a second way the same fixed-geometry equivalent loses uniformity.

## Cutoff and implementation bounds

The integrand defining $R$ is dominated by $w e^{-w^2/2-(k-v)w}$, because $\sinh(z)/z=(1/2)\int_{-1}^1e^{zu}\,du\le e^z$ for $z\ge0$. For $a=k-v\ge0$, the implementation truncates at $\min(14,48/a)$ when $a>0$, or at 14 when $a=0$. At 14 the omitted integral is at most $e^{-98}$. If $48/a<14$, it is at most $49e^{-48}/a^2\le49e^{-48}$.

If $a<0$, put $c=-a$ and normalize out $e^{c^2/2}$. On $[0,\infty)$, the normalized integrand is bounded by $w e^{-(w-c)^2/2}$. The implementation integrates $|w-c|\le14$. Since $c\le64$, the omitted normalized integral is at most $(1+2c/14)e^{-98}<11e^{-98}$, by a Gaussian tail estimate. The degree-eight truncation of $\sinh(z)/z$ used for $|z|\le1/4$ has real-arithmetic remainder bounded by its first neglected term times $1/(1-1/2496)$. The Fraction check verifies these small constants. They account for analytic truncation; the independent numerical audit addresses quadrature and binary32 rounding.
