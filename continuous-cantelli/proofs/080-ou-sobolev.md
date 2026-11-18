# Low Hermite spectrum, OU contraction, and fractional roughness

Let \(a=\varphi^2\) be a bounded nonnegative exact Cantelli variance, \(m=\mathbb E_\gamma a\), \(h=a-m\), and \(\epsilon=\|h\|_\infty>0\). The normalized transform identity is
\[
\mathbb E_\gamma[e^{zX-z^2/2}(e^{z^2h(X)/2}-1)]=0.
\]
It holds for real \(z\) by exactness and for complex \(z\) by analyticity. Local uniform domination for the Gaussian source and bounded \(h\) permits each finite Taylor computation below. Let \(e_n=He_n/\sqrt{n!}\) be the normalized probabilists' Hermites and write \(h_n=\langle h,e_n\rangle\).

## Exact coefficients and low-frequency bound

Divide the normalized transform by \(z^2/2\) and use the Hermite generating function. For each \(n\),
\[
h_n=-\sum_{k=1}^{\lfloor n/2\rfloor}
 \frac{\sqrt{n!/(n-2k)!}}{2^k(k+1)!}
 \langle h^{k+1},e_{n-2k}\rangle.\tag{1}
\]
The negative sign comes from moving all terms beyond the linear term in \(h\) to the other side. In particular \(h_0=h_1=0\) and \(h_2=-\sqrt 2\,\mathbb Eh^2/4\), agreeing with the fourth-moment identity.

Let \(P_{\le N}\) denote Hermite projection. For each \(k\) in (1), coefficient shift by \(2k\) with the displayed diagonal has operator norm at most \(\sqrt{N!/(N-2k)!}/[2^k(k+1)!]\) on the retained degrees. Apply the triangle inequality in coefficient \(\ell^2\), followed by \(\|h^{k+1}\|_2\le\epsilon^k\|h\|_2\) and the bound \(\sqrt{N!/(N-2k)!}\le N^k\). This gives
\[
\|P_{\le N}h\|_2\le G(N\epsilon/2)\|h\|_2,\qquad
G(u)=\sum_{k\ge1}\frac{u^k}{(k+1)!}=\frac{e^u-1-u}{u},\quad G(0)=0.\tag{2}
\]
The exact diagonal bound permits retaining just the norm of \(P_{\le N-2k}h^{k+1}\) in each corresponding term. If \(\eta=\|h^2\|_2/(\epsilon\|h\|_2)\), then replacing \(\|h^{k+1}\|_2\) by \(\epsilon^{k-1}\|h^2\|_2\) replaces \(G\) in (2) by \(\eta G\).

## Layer-cake Sobolev bounds

Define the spectral probability measure that gives mass \(|h_n|^2/\|h\|_2^2\) to \(\epsilon n\). By (2) its distribution function is bounded by \(\min(1,G(x/2)^2)\) for every \(x\ge0\), using \(N=\lfloor x/\epsilon\rfloor\). Let \(c_*=2u_*\), where the unique positive \(u_*\) satisfies \(e^{u_*}=1+2u_*\). For the nonnegative number operator \(\mathcal N e_n=n e_n\), the layer-cake formula gives, for every \(s>0\),
\[
\epsilon^s\|\mathcal N^{s/2}h\|_2^2 \ \ge\ C_s\|h\|_2^2,\qquad
C_s=s\int_0^{c_*}x^{s-1}[1-G(x/2)^2]dx.\tag{3}
\]
The statement allows infinite energy. A numerical value is unnecessary for an explicit constant: for \(0\le x\le2\), positivity of the series gives \(G(x/2)\le(e-2)x/2<3x/8\). Thus \(C_s\ge2^s[1-(9/16)s/(s+2)]\). In particular \(C_1\ge13/8\). If \(h\) has a weak Gaussian \(L^2\) derivative, Gaussian integration by parts and Parseval give \(\|h'\|_2^2=\sum_n n|h_n|^2\), since \(\langle h',e_{n-1}\rangle=\sqrt n h_n\). Conversely finite spectral energy gives the weak derivative by taking limits of the Hermite partial sums and their derivatives. Hence (3) includes \(\epsilon\|h'\|_2^2\ge(13/8)\|h\|_2^2\) without a pointwise differentiability assumption.

A check on the constants is \(C_2=2 C_1\). Differentiating \(uG=e^u-1-u\) gives \(uG'=u+(u-1)G\), and therefore \([u(1-G)^2]'=(1-2u)(1-G^2)\). Integrating over \([0,u_*]\) proves the relation.

## An independent direct OU contraction

Write \(T_t=e^{-t\mathcal N}\). Before taking \(N\to\infty\) in (1), apply \(T_t\). Each shifted diagonal is bounded by \(\sup_{n\ge0} e^{-tn}n^k\le(k/(et))^k\). Thus at \(t=\epsilon/2\),
\[
\|T_{\epsilon/2}h\|_2\le\|h\|_2\sum_{k\ge1}\frac{(k/e)^k}{(k+1)!}=(e-2)\|h\|_2.\tag{4}
\]
To evaluate the sum, the implicit series \(U(z)=ze^{U(z)}\) has coefficient \(k^{k-1}/k!\) by the elementary Lagrange inversion formula. For \(z<e^{-1}\), integrate \(uU'(u)\) from 0 to \(z\) and divide by \(z\): the resulting series is \(\sum_{k\ge1}k^kz^k/(k+1)!\), and the substitution \(u=Ue^{-U}\) gives \((e^U-1-U)/U\). As \(z\uparrow e^{-1}\), \(U(z)\uparrow1\), so monotone convergence proves the sum in (4). The finite projection estimate then permits passing \(N\to\infty\). Jensen applied to \(e^{-\epsilon n}\) in \(\|T_{\epsilon/2}h\|_2^2\) gives the independent direct derivative bound \(\epsilon\|h'\|_2^2\ge-2\log(e-2)\|h\|_2^2\). It has a weaker constant than (3). The stronger OU contraction from the distribution bound in (2) is
\[
\frac{\|T_{b\epsilon}h\|_2^2}{\|h\|_2^2}
\le 1-2b\int_0^{c_*}e^{-2bx}[1-G(x/2)^2]dx,\qquad b>0.
\]
Neither calculation assumes that an OU-smoothed graph remains an exact deterministic variance graph.

## Besov and physical increments

For a fixed \(0<c<c_*\), at every \(\epsilon\le t\le2\epsilon\) the high spectral mass \(\epsilon n>c\) and its lower bound give
\[
\|(I-T_t)^k h\|_2\ge(1-e^{-c})^k\sqrt{1-G(c/2)^2}\,\|h\|_2.\tag{5}
\]
For any \(s>0\), integer \(k>s/2\), and \(0<q\le\infty\), this yields an OU Besov increment lower bound
\[
\epsilon^{s/2}\left\|t^{-s/2}\|(I-T_t)^kh\|_2\right\|_{L^q((0,\infty),dt/t)}\ge c_{s,k,q}\|h\|_2.
\]
Integrate just over \([\epsilon,2\epsilon]\), or take a supremum there when \(q=\infty\). The same statement using the truncated time interval \((0,1)\) holds for \(\epsilon\le1/2\).

One may also obtain an ordinary centered difference seminorm with Gaussian midpoint weight. Couple \(X,Y\) with standard Gaussian margins and correlation \(e^{-t}\) at \(t=\epsilon\). Then \(M=(X+Y)/2\) and \(D=Y-X\) are independent, with variances \((1+e^{-t})/2\le1\) and \(2(1-e^{-t})\le2\epsilon\), respectively. The expectation \(E|h(Y)-h(X)|^2=2\langle(I-T_t)h,h\rangle\) is bounded below by a fixed positive constant times \(\|h\|_2^2\) from the same high spectral mass. The midpoint density is at most \(\sqrt2\phi(m)\). For \(0<s<1\), the difference density is at most \(C_s\epsilon^s |d|^{-1-2s}\), by maximizing the Gaussian times \(|d|^{1+2s}\). Consequently
\[
\epsilon^s\int\phi(m)dm\int\frac{|h(m+d/2)-h(m-d/2)|^2}{|d|^{1+2s}}dd
\ge c_s\|h\|_2^2.
\]
The positive comparison constants are universal for fixed orders. This conclusion is an increment bound and does not require pointwise differentiability.

## Executed verification

`hermite_check.py` verifies 49 raw transform coefficient identities through order 14 using SymPy, and independently calibrates the fractional constants and OU bounds with 50- and 75-digit mpmath. The analytic integral gives \(C_1\approx1.84746413904031077\), and at \(b=1/2\) the stronger squared OU contraction bound is approximately \(0.1854712448218525\). The derivative bound with constant \(13/8\) is exact by the analytic estimate above.
