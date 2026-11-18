# Positive sine sensitivity on a bounded future hull

Let \(I=[a,b]\), \(L=b-a\), \(u=(x-a)/L\), and \(F=\{f>s\}\cap I\). Under either hull hypothesis in R170, the future-transform identity and \(Q_s\le p_s\,dx\) imply that
\[
R_{\lambda,\vartheta}
=\frac{\int_F e^{\lambda^2 H(x)}\sin(\lambda u+\vartheta)\,g(x)\,dx}
{\int_I\sin(\lambda u+\vartheta)p_s(x)\,dx}\le1,
\qquad H(x)=\frac{f(x)-s}{2L^2},
\]
for \(0<\lambda=kL\le\pi\), \(0\le\vartheta\le\pi-\lambda\). Here \(g=p_1\), \(p_s\) has variance \(s\), and the physical height is \(h=f-s\). In physical frequency the exponential is \(\exp(k^2h/2)\).

For \(0<\lambda<\pi\) the two endpoint phases suffice, since
\[
\sin(\lambda u+\vartheta)
=\frac{\sin(\lambda+\vartheta)}{\sin\lambda}\sin(\lambda u)
+\frac{\sin\vartheta}{\sin\lambda}\sin(\lambda(1-u)).
\]
The coefficients are nonnegative, so the ratio lies between the two endpoint ratios, by weighting their positive denominators. For \(\lambda=0\) the constant test gives the mass budget. The normalized limits of the extreme phases give the linear weights \(u\) and \(1-u\). At \(\lambda=\pi\) only the first Dirichlet sine remains.

## Strict continuous example

Take \(s=1\) and \(I=[-1/2,1/2]\). Use two disjoint raised strips \(0.49<|x|<0.5\), with height \(h=3/2\) throughout \(0.4905\le|x|\le0.4995\), and linear ramps down to height zero at both endpoints of each strip. Off the strips, start at \(f=1\). It is possible to make \(f(0)=0\) and \(\int f\,g=1\) by subtracting \((1-|x|/\rho)_+\) with \(0<\rho<0.1\). Its weighted mass increases continuously from zero and at radius 0.1 exceeds the strips' total added Gaussian mass, by the interval certificate. This adjustment does not change the strict superlevel set at \(s=1\). The profile is nonnegative, continuous, bounded and even.

All admissible phases give the same ratio for this even profile, as the odd component integrates to zero. Thus the value can be tested by the centered weight \(\cos(kx)\). On \(I\), comparison of \(g\) with its center value loses a factor at most \(e^{1/8}\) in either direction. At the low-frequency endpoint the mass and both linear ratios are at most \(0.02e^{1/8}<0.022663\). At the first sine, use the raised strips as an upper envelope of the source:
\[
R_\pi\le e^{3\pi^2/4+1/8}(1-\cos(\pi/100))<0.916763 .
\]
At \(k=14/5\), keep only the plateau:
\[
R_{14/5}\ge
e^{3(14/5)^2/4-1/8}
\frac{\sin((14/5)\cdot0.4995)-\sin((14/5)\cdot0.4905)}{\sin(7/5)}
>1.483672 .
\]
These are outward interval bounds. Direct 65-digit integration of the full ramps gives the respective ratios \(0.0184339565465\), \(0.669751164198\) and \(1.558188909717\). The intermediate frequency supplies an exclusion beyond the mass and first-sine tests for two disconnected upper peaks.

The improvement can grow without bound if the allowed heights also grow. For two symmetric strips of width \(w\) with a fixed plateau fraction, on \(L=1\) set the physical height \(h=\alpha_{\rm phys}\log(1/w)\). The first sine ratio is \(O(w^{2-\alpha_{\rm phys}\pi^2/2})\) because the sine has a simple zero at the endpoints. At any fixed \(0<k<\pi\) the positive centered sine is bounded away from zero on the strips and the plateau gives a lower bound of order \(w^{1-\alpha_{\rm phys}k^2/2}\). Therefore for \(\pi/\sqrt2<k<\pi\) and \(2/k^2<\alpha_{\rm phys}<4/\pi^2\), the first and mass ratios tend to zero while the interior ratio diverges. Equivalently the normalized amplitude \(H\) uses \(\alpha_{\rm norm}=\alpha_{\rm phys}/2\) on \(L=1\). This comparison of necessary tests does not require an exact Gaussian-output graph.

## MacBook Metal scan

The numerical profiles use two strips attached to opposite endpoints, heights \(h_i=2L^2H_i\), and 5% linear ramps with a 90% plateau. The strata cover symmetric, asymmetric, varying \(s\), and a neighborhood of the fixed explicit witness. Each case uses 65 frequencies from 0 to \(\pi\) and both extreme phases. The \(\lambda=0\) ratios are the two linear limits, with the mass computed separately. Each numerator uses eight Gauss–Legendre nodes on each ramp or plateau, and the denominator uses eight on the hull.

Among 435,945,472 finite cases, 197,917,248 had all three \(\lambda=0\) and the \(\lambda=\pi\) ratios below .999. Of those, 78,777,350 exceeded 1.001 at an interior grid frequency. In 28,982,361 of that subset only one extreme phase exceeded the threshold and the other remained below .999 across the interior grid. These counts characterize the scanned strata and margins. The strongest asymmetric retained sample's rejecting phase has a locally refined CPU ratio of 29.46676 at \(\lambda\approx3.016717\), with first-sine ratio about .9878945. The ratios on all saved profiles were independently replayed by adaptive SciPy and 65-digit mpmath, agreeing within \(8.660\times10^{-14}\); the largest absolute Metal discrepancy over the replayed fields is \(3.484\times10^{-5}\). Local and remote numerical payloads match exactly.

The Apple M5 Pro main scan ran 11:03:44.534–11:04:21.544 UTC on 2026-11-18, using 33.119185040 native Metal seconds over 37.009681833 s elapsed. The device-wide utilization histogram is in `summary.json`. Both smoke and main returned zero and both PIDs were absent afterward. Thirteen downloaded files match their remote SHA-256 hashes.

### Reproduce

```sh
cache="$HOME/CodexRemoteMath/cache"
mkdir -p "$cache/swift-modules" "$cache/cc170-temp"
TMPDIR="$cache/cc170-temp" swiftc -O -module-cache-path "$cache/swift-modules" scan.swift -o sine_scan
sh run.sh main 37 131072 19093169 65
"$HOME/.venvs/math/bin/python" check_cpu.py main.json main.cpu.remote.json
```

The wall cutoff makes the number of batches system-throughput dependent. The saved `main.json` and `check_cpu.py` reproduce the reported CPU and interval checks.
