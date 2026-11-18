# Gaussian time-Marchaud crossings

For \(0<\kappa<1\), the time kernel is
\[
F_\kappa(y)=\int_0^\infty [p_1(y)-p_{1+t}(y)]t^{-1-\kappa}\,dt ,
\quad p_v(y)=\frac{e^{-y^2/(2v)}}{\sqrt{2\pi v}} .
\]
At starting variance \(s>0\), scaling changes the kernel by \(s^{-1/2-\kappa}\) and the positive zero by \(\sqrt s\).

Fourier transformation gives \(2^{-\kappa}\Gamma(1-\kappa)|\xi|^{2\kappa}e^{-\xi^2/2}/\kappa\). Inversion and Kummer's transformation, with \(x=y^2/2\), yield
\[
F_\kappa(y)=
\frac{\Gamma(1-\kappa)\Gamma(1/2+\kappa)}{\kappa\sqrt2\,\pi}\,
e^{-x}{}_1F_1(-\kappa;1/2;x).
\]

## Uniqueness, monotonicity, and endpoint scales

The factor determining the sign for \(x\ge0\) is \(1-S_\kappa(x)\), where
\[
S_\kappa(x)=\sum_{n\ge1}
\frac{\kappa(1-\kappa)_{n-1}x^n}{(1/2)_n n!}.
\]
Every coefficient is positive and the linear term alone is unbounded, so there is exactly one positive simple zero \(y_\kappa\). At \(y=1\), \(p_{1+t}(1)<p_1(1)\) for \(t>0\), giving \(y_\kappa>1\).

For each fixed \(y>1\), \(p_1(y)-p_{1+t}(y)\) is negative for a first interval of positive \(t\) and positive after its single crossing \(t_*(y)\). This follows from the variance derivative of \(p_v(y)\), which changes sign from positive to negative at \(v=y^2\). At \(y=y_\kappa\), increasing \(\kappa\) and scaling the new integrand by \(t_*(y)^{\Delta\kappa}\) multiplies the former integrand by a factor greater than one in the negative interval and less than one in the positive interval. Thus the root strictly decreases with \(\kappa\).

Let \(e=1-\kappa\) and \(x_\kappa=y_\kappa^2/2\). Split
\[
S_\kappa(x)=2\kappa x+\kappa e\,U_\kappa(x),\qquad
U_\kappa(x)=\sum_{n\ge2}\frac{(2-\kappa)_{n-2}x^n}{(1/2)_n n!}.
\]
Continuity of the positive series at \((1,1/2)\), or the implicit function theorem, gives
\[
\lim_{\kappa\uparrow1}\frac{y_\kappa-1}{1-\kappa}
=\frac{1-U_1(1/2)}2
=0.41068613464244799769\ldots .
\]
The equality follows by writing \(x_\kappa=1/2+e\,r\) in the root equation, which reduces to \(2\kappa r+\kappa U_\kappa(1/2+er)=1\).

For \(\ell=\log(1/\kappa)\to\infty\), use
\(S_\kappa(x)=2\kappa\int_0^x{}_1F_1(1-\kappa;3/2;u)\,du\). Euler's beta integral near \(t=1\) gives, uniformly for small \(\kappa\),
\[
{}_1F_1(1-\kappa;3/2;u)
=\frac{\Gamma(3/2)}{\Gamma(1-\kappa)}e^u u^{-\kappa-1/2}
\left[1+\frac{\kappa(\kappa+1/2)}u+O(u^{-2})\right].
\]
Integrating up to \(x\) adds a relative term \((\kappa+1/2)/x\). The root equation becomes \(x=\ell+\tfrac12\log(x/\pi)-\tfrac1{2x}+O(x^{-2})\), since terms involving \(\kappa\) are exponentially small once \(x\sim\ell\). Solving by expansion gives
\[
y_\kappa^2
=2\ell+\log(\ell/\pi)
+\frac{\log(\ell/\pi)-2}{2\ell}
+O\!\left(\frac{(\log\ell)^2}{\ell^2}\right).
\]

## Native GPU and CPU check

The Metal kernel solves \(S_\kappa(x)=1\) by parallel bisection. In its logarithmic-small-\(\kappa\) stratum it retains \(\log\kappa\) separately, rescales the positive term sum by powers of two before overflow, and tests the sign of \(\log S_\kappa(x)\). In its near-one stratum it bisects \(r=(x-1/2)/(1-\kappa)\) using the cancellation-free equation above. The reported \((y-1)/(1-\kappa)=2r/(y+1)\) stays resolved even when the raw Float32 root rounds to one.

On both CPUs, the 108-digit root equation uses
\[
S_\kappa(x)=2\kappa x\,{}_2F_2(1,1-\kappa;2,3/2;x).
\]
Factoring out \(\kappa\) is important. The initial direct \({}_1F_1(-\kappa;1/2;x)\) CPU method returned 1 at a saved root near \(\log\kappa=-300.35\), effectively losing the small parameter. The revised method keeps direct Kummer and SciPy comparisons where their ranges are suitable. It also evaluates the independent positive Euler integral for \({}_1F_1(1-\kappa;3/2;u)\) and the direct time-kernel integral at selected roots after \(u=t/(1+t)\) and integration by parts:
\[
F_\kappa(y)/p_1(y)
=\frac1{2\kappa}\int_0^1 e^{xu}u^{-\kappa}(1-u)^{\kappa-1/2}
[1-2x(1-u)]\,du .
\]
These checks analytically remove the singular constant at \(u=0\). On the checked probes the normalized Euler residual is at most \(1.171\times10^{-97}\) and the normalized time residual at most \(3.453\times10^{-108}\).

| \(\kappa\) | First positive zero |
| ---: | ---: |
| 0.01 | 3.066350109536241 |
| 0.10 | 2.094081693507304 |
| 0.25 | 1.643800190452846 |
| 0.50 | 1.306929727719281 |
| 0.75 | 1.122285077063811 |
| 0.90 | 1.043840586888272 |
| 0.99 | 1.004132791121480 |

The main scan evaluated 786,444 finite probes in 49 native command buffers, with three dense strata and 12 anchors. It reached \(\log\kappa=-600\) and \(\log(1-\kappa)=-80\). No monotonicity reversal above the scan tolerance of \(2\times10^{-6}\) was observed. Across 63 saved CPU comparisons, the maximum absolute GPU error in \(y\) was \(2.388\times10^{-6}\). The maximum error in the near-one scaled gap was \(1.077\times10^{-7}\). Both host CPU payloads match exactly. At \(\ell=600\), the high-precision root is \(34.716781337583253\); the refined asymptotic reduces the error in \(y^2\) from about \(2.70667\times10^{-3}\) to \(3.49428\times10^{-6}\).

## Receipt and reproduction

The Apple M5 Pro production ran 12:06:42.368–12:06:42.441 UTC on 2026-11-18, using 0.042957417 s of Metal command time over 0.072726208 s elapsed. One device-wide utilization sample was 99%. Both smoke and main returned zero, both PIDs were absent at the post-run check, and 13 archived remote hashes match. The GPU launch initially used a remote directory that already held completed CPU work and overwrote that run's wrapper name. The GPU files were subsequently copied into the separate GPU archive, and the CPU worker retained its results and distinct wrapper copy.

```sh
cache="$HOME/CodexRemoteMath/cache"
mkdir -p "$cache/swift-modules" "$cache/cc192-temp"
TMPDIR="$cache/cc192-temp" swiftc -O -module-cache-path "$cache/swift-modules" scan.swift -o marchaud_scan
sh run.sh main 262144 16384
"$HOME/.venvs/math/bin/python" check_cpu.py main.json main.cpu.remote.json
```
