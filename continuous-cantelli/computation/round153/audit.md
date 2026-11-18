# Quantitative negative heat test for a nowhere-dense source

Set \(\varepsilon=1/8\), \(s=7/8\) and \(c=15/16\). The round-128 interval is \(J=[2,2+\varepsilon^{3/2}/100]\). Take the closed subinterval
\[
 I=[a-r,a+r],\qquad a=2+1/5000,\quad r=3/20000 .
\]
Its right displacement from 2 is \(7/20000<\varepsilon^{3/2}/100\), as checked by squaring positive rationals.

Starting from \(I\), at stage \(n\ge1\) remove from every surviving interval its open central subinterval of relative length \(2^{-n}/256\). Let \(K\) be the intersection of the successive finite unions. It is compact. Every stage-\(n\) component has length at most \(2^{-n}|I|\), so \(K\) has empty interior and is nowhere dense. Because the total removed length is at most \(|I|\sum_{n\ge1}2^{-n}/256=|I|/256\), the positive Lebesgue measure satisfies \(|K|\ge255|I|/256\).

Use \(\eta=c\,1_K\gamma\) on layer \(s\), \(P_\varepsilon\eta\) on \(1+\varepsilon\), and \(\gamma-\eta-P_\varepsilon\eta\) on the middle layer at time 1. The layers sum to \(\gamma\). By completion of the square, for \(x\in I\) and all \(y\),
\[
 \frac{p_\varepsilon(y-x)}{g(y)}\le\varepsilon^{-1/2}e^{x^2/[2(1-\varepsilon)]}<3\cdot27.
\]
Here \(g\) is the standard Gaussian density, \(\varepsilon^{-1/2}<3\), \(x^2/[2(1-\varepsilon)]<3\) and \(e<3\). Since \(\eta(\mathbb R)<|I|/2\), the middle layer dominates \((1007/20000)\gamma\). For any real \(z\), with \(q=z^2/2\), the change in the exponential transform from the constant layer is the transform of \(\eta\) multiplied by
\[
 e^{-qs}+e^{-q(1+\varepsilon)}e^{q\varepsilon}-e^{-q}(1+e^{q\varepsilon})=0 .
\]
Thus this positive source preserves the spatial marginal and the complete transform identities. Positive mass splits among distinct layers and the source is non-graph.

## Smooth test and strict rational time

Use \(\psi=0\) outside \(I\); inside, with \(u=(x-a)/r\), set
\[
 \psi(x)=\exp\bigl(1-1/(1-u^2)\bigr).
\]
The smooth flat extension at the endpoints is standard from the exponential decay compared with powers of \((1-u^2)^{-1}\). The maximum is one, and \(\psi\ge e^{-1/3}\ge2/3\) on the middle half of \(I\).

For \(x\ge2\), use \(-\log(1-\varepsilon)\le \varepsilon/s\) to bound the free density ratio
\[
 p_s(x)/g(x)\le e^{-3\varepsilon/(2s)}\le (1+3\varepsilon/(2s))^{-1}=14/17<7/8.
\]
Write \(g_{\min}=g(a+r)\), \(g_{\max}=g(a-r)\). The rational Taylor checks give \(g_{\min}>1/30\) by \(g(0)>1/3\), \((a+r)^2/2<9/4\), and \(e^{9/4}<10\). This gives the explicit mass bounds \(\gamma(I\setminus K)<3/5120000\) and \(\gamma(K)>51/5120000\) from \(g_{\max}<1/2\) and \(g_{\min}>1/30\). Also \(g_{\max}/g_{\min}<5/4\), since the corresponding log ratio is \(\Delta=((a+r)^2-(a-r)^2)/2<1/100\) and \(e^\Delta\le(1-\Delta)^{-1}\). Hence
\[
 \langle Q_s,\psi\rangle
 \le -(1/16)(2/3)r g_{\min}+c\,|I\setminus K|g_{\max}
 \le -\left(1/24-75/8192\right)r g_{\min}
 <-r/960=-1/6400000.
\]

For the heat error put \(A=(1-u^2)^{-1}\). On the positive half the magnitude of the derivative with respect to \(u\) equals \(2u A^2e^{1-A}\), bounded by \(8/e<4\) from the maximum of \(A^2e^{1-A}\). This magnitude increases and then decreases because its logarithmic derivative times \(u\) equals \(-2A^2+6A-3\), which changes sign once for \(A\ge1\). Reflection and scaling consequently give \(\|\psi''\|_1<16/r\).

For \(0<h<\varepsilon\), \(Q_{s+h}=P_hQ_s\) as no new layer has occurred. The density of \(Q_s=p_s\,dx-\eta\) has \(\|Q_s\|_\infty<1/2\), since each of the two nonnegative densities is smaller than \(1/2\). Heat differentiation and \(L^1\) contraction give \(\|(P_h-I)\psi\|_1\le h\|\psi''\|_1/2\). Thus the pairing changes by less than \(4h/r\). In particular for the strict rational time
\[
 t=s+10^{-12}=875000000001/1000000000000
\]
we obtain
\[
 Q_t(\psi)< -1/6400000+1/37500000=-311/2400000000.
\]

## Reproduction

`certify.py` checks the constants with exact rational arithmetic. Its exponential enclosures use finite Taylor sums and geometric bounds on their remaining terms. It encloses \(\pi\) through alternating arctangent series and the Machin identity. The same script constructs and checks every rational surviving interval through a requested stage, independently checking the total removed length against both the product formula and the analytic length bound. The finite enumeration is a reproducibility check. The infinite nowhere-density and mass assertions follow from the proof above.

`bash run.sh` runs three bounded MacBook processes at stages 8, 12 and 15, with 32,768 surviving intervals checked at the final stage. Python 3.12.13 executed the calculation from 09:46:43.763619 to 09:46:44.224640 UTC. All three workers and the supervisor exited zero, and `final-ps.txt` confirms every owned PID was absent. `summary.json` contains the exact constants and timing; all 16 transferred file hashes were checked.
