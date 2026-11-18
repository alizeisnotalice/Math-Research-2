# R117 small-set Hellinger and noncentral overlap

## Necessary center energy

For a bounded exact graph let \(P_t(dx)=e^{tx-t^2 f(x)/2}\gamma(dx)\), so \(P_t(\mathbb R)=1\) for every real \(t\). Fix \(t>1\), \(A=2\log t\), and a fixed union \(E\) of low-level windows. In particular \(E\) may be a subset of a small sublevel set. Write
\[
p=P_t(E),\quad r=1-p,\quad D=\frac{t}{\sqrt A}(X-tf(X)),\quad F=\frac{t^2}{A}f(X),\quad
\mathcal E=E_t F.
\]
For \(h=\lambda t/\sqrt A\), the two noncentral overlap masses are
\[
q_\pm=P_{t\pm h}(E)=p\,E_t[e^{\pm\lambda D-\lambda^2F/2}\mid E].
\]
The Hellinger affinity of the two restricted complements is \(E_t[1_{E^c}e^{-\lambda^2F/2}]\). Jensen from below and Cauchy–Schwarz from above give
\[
r e^{-\lambda^2 E_t[F\mid E^c]/2}
\le E_t[1_{E^c}e^{-\lambda^2F/2}]
\le \sqrt{(1-q_-)(1-q_+)}.
\]
As a result, \(0<p<1\) and \(\lambda\ne0\) require
\[
\mathcal E\ge
\mathcal H_\lambda :=p E_t[F\mid E]+\frac{r}{\lambda^2}
\log\frac{r^2}{(1-q_-)(1-q_+)}.
\]
The independent mass-energy bound is \(\mathcal E\ge p E_t[F\mid E]\). The infinitesimal limit of \(\mathcal H_\lambda\) is
\[
\mathcal V = p E_t[D^2\mid E]+\frac{p^2}{r}E_t[D\mid E]^2 .
\]
Indeed, expand the two exponentials to quadratic order and add the two logarithms. The score identities \(E_tD=0\), \(E_t D^2=\mathcal E\) also prove \(\mathcal E\ge\mathcal V\) directly. The mean-score constraint additionally yields each directional Jensen bound by conditioning the full shifted likelihood on \(E^c\). Their average is \(\mathcal H_\lambda\) and their maximum can strengthen it.

## Model and all-shift certificate

Use \(z=(t x-A)/\sqrt A\). The two windows in \(z\) have midpoint parameters \(m-d,m+d\) and radii \(w_L,w_R\), with \(w_L,w_R<d\). For the local coordinate \(y=z-(m\pm d)\), put \(\rho_j(y)=y^2/\sigma^2+\delta_j\), \(\delta_L=0\). Prescribe
\[
t^2 f(x)=A+2\sqrt A\,z+c+\rho_j(y),
\qquad x=(A+\sqrt A\,z)e^{-A/2}.
\]
The normalization
\[
c=\log A+2\log\left(\frac{\sum_j\int_{-w_j}^{w_j}e^{-(x^2+\rho_j(y))/2}\,dy}
{p\sqrt{2\pi}}\right)
\]
makes \(P_t(E)=p\) on this restricted profile. On the windows,
\[
D=-z-\frac{c+\rho_j}{\sqrt A},\qquad
F=1+\frac{2\sqrt A\,z+c+\rho_j}{A}.
\]
For each shift, the combined integrand is an exponential of a quadratic in \(y\). Its integral is a completed-square noncentral Gaussian probability. This also permits direct adaptive quadrature and validated composite midpoint integration.

A strict example uses \(A=16384\), \(m=0\), \(d=2.05\), \(w_L=w_R=.01\), \(\sigma=1\), \(\delta_R=0\), and \(p=.228\). The interval verifier obtains
\[
p E_t[F\mid E]<.228062,\qquad
\mathcal V<.958525,\qquad
\mathcal H_{1.25}>1.10171123522.
\]
For every real shift, \(P_{t+\lambda t/\sqrt A}(E)<.963650\). To see this last bound, the validated inputs give for each direction an active window with directional score at most \(a>0\), the other window's directional score at most \(-b<0\), depth lower bounds \(F_a,F_b>0\), and a bound \(w\) on each window's conditional weight. Put \(s=|\lambda|\). The overlap is at most
\[
pw[e^{sa-s^2F_a/2}+e^{-sb-s^2F_b/2}].
\]
On \(0\le s\le1\), drop the quadratic factors and use endpoint convexity to bound the bracket by \(\max(2,e^a+e^{-b})\). On \(s\ge1\), bound it by \(e^{a^2/(2F_a)}+e^{-b-F_b/2}\). All four directional/region bounds are recorded in the outward certificate. Thus at energy \(\mathcal E\le1\) the joint test excludes this local geometry even though its individual overlap caps at all real shifts, the moment lower bound and the mass-energy lower bound pass. The complement may contain additional low sources.

The limiting symmetric window model \(A\to\infty\) gives the simple condition \(p(E[\cosh(\lambda z)\mid E]-1)\le e^{\lambda^2/2}-1\) at energy at most one. For the same narrow wells, \(\mathcal V\le .9675408\) and the mass of every individual shift is less than .960274. At \(\lambda=1.25\), the left side exceeds the right by at least .05671033, from monotone bounds on \(\cosh\) over the wells. The finite certificate above retains this strict separation.

## Native run, verification, and scope

The M5 Pro production tests three distributions over \(A\in[16,16384]\) and well separation, asymmetry, widths, depth and shape. Each accepted profile evaluates eight shifts from .125 to 2.25 with 384 midpoint panels per window, using cancellation-safe residual and logarithm series. For summary ranking the shader retains shift pairs whose two masses are at most .98. It screened 191,889,408 inputs, accepted 183,919,396 with positive computed depth, and evaluated 1,471,355,168 candidate shift pairs. For 5,256,475 accepted inputs the numerical score-variance and mass-energy lower bounds were below one while at least one safe complement-Hellinger bound was above one.

The strongest ranked sample with both shifted masses below .9 has a completed-square replay of \(q_+=.60296463\), \(q_-=.89001707\), \(\mathcal V=.98230402\), and \(\mathcal H_{.85}=3.38083055\). It is certified outward with \(\mathcal H_{.85}>3.38077532\). A second targeted sample is certified at \(\mathcal H_{.85}>4.76209263\) with \(q_-<.978059\). The interval checks bound the composite midpoint error by \((2w)^3\sup |g''|/(24n^2)\), with \(g=e^{C+Ly-By^2/2}\) and \(n=1024\). Polynomially weighted central moments use the integration-by-parts recurrence \(B I_{j+1}=L I_j+j I_{j-1}-[y^jg(y)]_{-w}^w\). The CDF replay and double-precision adaptive quadrature agree by at most \(1.274\times10^{-11}\) for the saved output fields across eight selected parameter sets. Local and remote CPU JSON results match exactly. The strict conclusions follow from the outward certificates and the analytic inequalities.

The production interval was 06:28:30.626–06:30:01.377 UTC on 2026-11-18, 90.750558333 s elapsed and 68.219081209 native Metal command seconds. IOAccelerator samples were 161 nonzero readings from 89% to 100%, plus four zero readings. The main and smoke runs returned zero. The first diagnostic stopped at Metal compilation because `half` was a reserved type name and was corrected before the successful smoke. All three monitored PIDs were absent after completion.

### Reproduce

```sh
cache="$HOME/CodexRemoteMath/cache"
mkdir -p "$cache/swift-modules" "$cache/cc117-temp"
TMPDIR="$cache/cc117-temp" swiftc -O -module-cache-path "$cache/swift-modules" scan.swift -o overlap_scan
sh run.sh main 384 68 1048576 35575387
"$HOME/.venvs/math/bin/python" check_cpu.py main.json main.cpu.remote.json
```

The active command-time cutoff makes the final batch count depend on device throughput. The saved `main.json` reproduces the selected numerical checks. Model constraints are necessary for exact bounded graphs, and the certificate excludes the stated local shape under its center energy bound.
