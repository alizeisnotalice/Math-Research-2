# A two-sided tail test for fixed support

Let \(F_f(t)=\int e^{tx-t^2 f(x)/2}\gamma(dx)\). Suppose \(f(0)=3/4\), \(f(x)=1\) for \(|x|\ge R\), \(R\ge2\), and \(f\) has Lipschitz constant at most \(1.500001\). With \(r=1/13\), the fixed origin pin gives \(f(x)\le 7/8\) when \(|x|\le r\). Moreover
\[
\gamma([-r,r])\ge \frac{2}{13\sqrt{2\pi}}e^{-1/338}>\frac1{20}.
\]
The last estimate uses \(1/\sqrt{2\pi}>1/3\) and \(e^{-1/338}>39/40\).

Subtract the constant graph, then multiply the residual by \(e^{t^2/2}\). On the pinned neighborhood the residual integrand has the lower bound \(e^{-|t|r}(e^{t^2/16}-1)\). On the entire support, wherever the integrand is negative it has absolute value at most \(e^{|t|R}\). Consequently
\[
e^{t^2/2}(F_f(t)-1) >
\frac1{20} e^{-|t|r}(e^{t^2/16}-1)-e^{|t|R}.
\]
For \(|t|\ge16R+18\), use \(e^{t^2/16}-1\ge\tfrac12 e^{t^2/16}\) and \(t^2/16-|t|(R+r)>|t|\). The right side is at least
\[
e^{|t|R}\big(e^{|t|}/40-1\big)>0 .
\]
The same argument gives a lower bound \(|t|-\log 40\) for the log positive-to-negative residual ratio if a negative part exists. The formula applies on both tilt wings. At fixed \(R\), it provides a common held-out parameter with a uniform positive residual over the entire compact support class, independently of any optimization.

The parameter needed in this proof grows with \(R\). The lower bound on the unrescaled residual also deteriorates with \(R\). An exhaustion by increasing supports with a common local modulus remains open.

## Meaning of the heat tests

With \(\varepsilon>0\), the Gaussian-smoothed heat-residual density is
\[
(Q_t*p_\varepsilon)(y)=p_{t+\varepsilon}(y)-
\int_{f(x)\le t}p_{t-f(x)+\varepsilon}(y-x)\,\gamma(dx).
\]
The formula includes the \(f(x)=t\) layer. A strictly negative value means \(Q_t\) fails positivity, as certified whenever one has rigorous bounds on the integrals. For the approximate profiles here it tracks a departure from the constant graph. Exact normalization requires the full transform residual to vanish.

To retain the directed interval description, set \(a(x)=2-f(x)\) and \(b(x)=f(x)-1\). The independently tested potential is
\[
U(y)=-\int b(x)\int_0^1p_{1-zb(x)}(y-x)\,dz\,\gamma(dx).
\]
Every update changes one graph endpoint at each \(x\). The vertical source remains one directed interval at all stages.
