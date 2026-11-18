# Weighted Cantor resonance windows

Let \(X=2+\sum_{n\ge1}\epsilon_n3^{-n}\), where the signs are independent uniform, and let \(\mu\) be its law. Write \(\sigma=g\mu\) for the standard Gaussian density weight. Set \(d=\log2/\log3\), \(P(w)=\prod_{n\ge1}\cos(w/3^n)\), and \(\xi_m=2\pi3^m\). Fourier transforms here use the phase \(e^{i\xi x}\). Removing the unimodular center phase gives \(F(w)=e^{-2iw}\widehat\sigma(w)\).

For \(\beta=1,2,3\), evaluate the observable
\[
I_b=e^{-D_L^2/(2b)}M_b
=\frac1{2\pi}\int e^{-b(\xi-D_L/b)^2/2}|\widehat\sigma(\xi)|^2\,d\xi,
\]
where
\[
D_L^2=2c\beta A L^{1-\beta},\quad s_L=(\beta-1)A L^{-\beta},
\quad b=A L^{-\beta}+s_L=\beta A L^{-\beta}.
\]
Here \(D_L^2/(2b)=cL\). Choose \(L=(\xi_m\sqrt{\beta A/(2c)})^{2/(\beta+1)}\) for exact resonance. The rounded rows use the nearest positive \(j\log3\) and keep the resulting center offset \(D_L/b-\xi_m\) in the Gaussian. Write \(I_b^0\) for the same quadratic integral centered at zero and \(W=b^{-1/2}\) for its standard deviation in frequency.

## Energy and exponent

The centered window is
\[
I_b^0=\frac1{\sqrt{2\pi b}}\iint e^{-(x-y)^2/(2b)}d\sigma(x)d\sigma(y)
\asymp W^{1-d}.
\]
Indeed the weighted Cantor law remains \(d\)-Ahlfors regular. Use a ball of radius \(1/W\) for the lower bound and Gaussian dyadic shells for the upper bound.

For a shifted window, split off the first \(m\) signs and put
\(A_m(w)=E[g(2+X_m)e^{iw X_m}]\), with \(X_m=\sum_{n\le m}\epsilon_n3^{-n}\) centered at zero. The discarded tail has size at most \(3^{-m}/2\), so freezing \(g\) at the prefix has absolute error \(O(3^{-m})\), uniformly in the frequency. The prefix phase at \(\xi_m\) is exactly one. By independence,
\[
F(\xi_m+w)=P(2\pi+3^{-m}w)A_m(w)+O(3^{-m}),\qquad
F(w)=P(3^{-m}w)A_m(w)+O(3^{-m}).
\]
On \(|w|\le\epsilon3^m\) the second multiplier is bounded away from zero. Therefore the first expression equals \([P(2\pi)+O(3^{-m}|w|)]F(w)+O(3^{-m})\). Integrate the difference of squared moduli against the Gaussian window. The term with \(|w||F(w)|^2\) is \(O(W^{2-d})\) by a Gaussian with double variance, and the term with \(|F(w)|\) is \(O(W^{1-d/2})\) by Cauchy–Schwarz. The outside region is exponentially suppressed since \(W3^{-m}=O(L^{-1/2})\). It follows that \(I_b/I_b^0=P(2\pi)^2+O(L^{-1/2})\) for exact resonance.

Rounding changes the Gaussian center by \(\Delta=O(\xi_m/L)\), or \(|\Delta|/W=O(L^{-1/2})\). The mean value theorem bounds the difference of the two Gaussian windows by \(C|\Delta|W^{-1}\exp(-w^2/(4W^2))\). The shifted energy with doubled variance also has the same \(O(W^{1-d})\) upper bound by the real-space formula and \(|\cos|\le1\). This gives the same rate for the rounded windows.

In particular \(P(2\pi)\ne0\) and
\[
\log M_b=cL+\frac{\beta(1-d)}2\log L+O(1),\qquad
\frac1L\log(e^{-(1-d)L}M_b)\longrightarrow c-1+d.
\]
At \(c=1-d\) the compensated quantity grows polynomially as \(L\to\infty\). The finite numerical values of \(I_b\) at the sampled critical orders are still below one.

## Native window and independent check

The GPU evaluates the weighted characteristic function by the convolution
\[
F(w)=\frac1{\sqrt{2\pi}}E[e^{2iU}P(w+U)],\quad U\sim N(0,1),
\]
using 24 normalized Gauss–Hermite nodes. At resonance it uses \(P_m(w+U)P(2\pi+3^{-m}(w+U))\), so the large base phase is never added to a small offset. The unshifted product uses the same prefix. Prefix trig evaluations use split Float32 scaling and a rotation for each Hermite node. The Gaussian energies use a uniform mesh of spacing 0.5 with tails clipped at eight standard deviations, and 12 diagnostics use spacing 0.25. Since \(|\widehat\sigma|^2\) is the Fourier transform of a difference measure on \([-1,1]\), Poisson summation gives Gaussian suppression of the infinite-mesh aliases at distances at least \(2\pi/\Delta\xi-1\). The final metric is computed as \(\log M_b=cL+\log I_b\), without exponentiating \(cL\).

The CPU independently computes the spatial difference kernel
\[
\frac{W}{\sqrt{2\pi}}\iint e^{-W^2(x-y)^2/2}
\cos((D_L/b)(x-y))\,d\sigma(x)d\sigma(y).
\]
It builds prefix-weight autocorrelations by FFT on exact ternary addresses, and adds the tail differences on the ternary digits \(-2,0,2\) of probabilities \(1/4,1/2,1/4\). The first Gaussian weight correction uses \(E[U\mid U-V]=(U-V)/2\), \(E[V\mid U-V]=-(U-V)/2\) for the two independent symmetric tails, which gives the factor \(1-\delta R\Delta_{\rm prefix}/2\). Further ternary-depth and tail refinements give a direct convergence check. A separate 17-digit address sum checks the peak Fourier power.

The 18 saved CPU pair checks cover both exact and rounded windows at \(m=2,7,12\). The largest relative shifted-window difference versus the GPU is \(6.838\times10^{-7}\); the centered-window difference is \(2.763\times10^{-7}\). The largest depth-refinement change in the CPU shift result is \(3.809\times10^{-7}\), and the largest GPU half-step change is \(1.127\times10^{-7}\). Both hosts produce identical numeric results.

For \(c=1-d,A=0.5,m=12\) the exact-resonance rows are:

| \(\beta\) | \(I_b/I_b^0\) | \(L^{-1}\log(e^{-(1-d)L}M_b)\) |
| ---: | ---: | ---: |
| 1 | 0.13796588597 | \(-1.8710\times10^{-6}\) |
| 2 | 0.13797787921 | \(-1.7285\times10^{-4}\) |
| 3 | 0.13810883303 | \(-1.7520\times10^{-3}\) |

The persistent ratio reference is \(P(2\pi)^2=0.1379657099587948\ldots\). At \(m=12\) the calculated compensated exponents have the predicted signs for \(c=0.25\) and \(c=0.55\) across all three \(\beta\).

## Receipt and reproduction

The Apple M5 Pro main scan evaluated 72,640,380 frequency sites, each with both Fourier windows, across 444 commands. It covered \(m=1,\ldots,12\), \(c=0.25,1-d,0.55\), \(A=0.5,2\), and the three \(\beta\) values. The run was 11:37:40.634–11:37:41.327 UTC on 2026-11-18, with 0.390141000 s of native Metal commands and 0.693387292 s elapsed. The short run yielded two device-wide ioreg samples, 9% and 70%. Both smoke and main exited zero and both PIDs were absent afterward. All 14 recorded remote SHA-256 hashes match.

```sh
cache="$HOME/CodexRemoteMath/cache"
mkdir -p "$cache/swift-modules" "$cache/cc184-temp"
TMPDIR="$cache/cc184-temp" swiftc -O -module-cache-path "$cache/swift-modules" scan.swift -o cantor_window
sh run.sh main 12 .5
"$HOME/.venvs/math/bin/python" check_cpu.py main.json main.cpu.remote.json
```

The 24-node Hermite input in `hermite.json` uses SciPy's `roots_hermitenorm(24)`, with the node weights multiplied by \(e^{2iu}/(2\pi)\).
