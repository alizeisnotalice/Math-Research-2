# R144 dyadic transfer scan for two oscillatory wells

Consider the frozen-mark transfer at \(d>0\),
\[
G_d(t)=t^2\int_M^\infty
p_{m/d}(t-m)m^{-2}[A+B\cos(2\sqrt m)]\,dm,
\qquad R=\sqrt t,\quad\phi=2R\bmod 2\pi,
\]
where \(p_v\) denotes the centered Gaussian density of variance \(v\). This is the local transfer in R144 after the realized mark \(x(m)-b\) is frozen at \(d=a-b\). The R144 bound \(x(m)-a=O(m^{-5/2}e^{-dm/2})\) makes this replacement exponentially accurate after normalization. On \(m\in[t/2,2t]\) this follows by differentiating the Gaussian in its mark, and on its complement the kernel tail is exponentially small. The numerical experiments set \(M=1\) and evaluate the transformed integral on a symmetric Gaussian window. At the sampled \(R\ge2\pi\), that window lies within the integration domain.

## Intermediate-phase and pair expansions

Set \(y=(m-t)/\sqrt m\), \(h=y/(2R)\), \(v=\sqrt{1+h^2}\), and \(u=v+h\). The transformed integrand relative to \(p_{1/d}(y)\) is
\[
\frac{A+B\cos(\phi+y+\eta)}{u^3v},\qquad
\eta=\frac{y^2}{2R(v+1)}.
\]
The paired weights at \(y,-y\) have even part \(1+4h^2\) and odd part \(-h(3+4h^2)/v\). Indeed \(u=e^{\operatorname{arsinh}h}\), so the triple-angle hyperbolic formulas give both identities. Consequently when \(M=0\) the nonoscillatory \(A\) part is exactly \(A(1+(1/d)/R^2)\).

With \(V=1/d\) and \(C=e^{-V/2}\), expand the weight \(1-\frac32y/R+y^2/R^2+O(R^{-3})\) and phase \(\eta=y^2/(4R)+O(R^{-3})\). Gaussian moments, from derivatives of \(\exp(-V\xi^2/2)\) at \(\xi=1\), yield
\[
G_d(t)=A+BC\cos\phi+\frac{BC}{R}\frac{5V+V^2}{4}\sin\phi
+\frac{AV+BC\,T(V)\cos\phi}{R^2}+O(R^{-3}),
\qquad T(V)=V+\frac{V^2}{32}-\frac{3V^3}{16}-\frac{V^4}{32}.
\]
Taylor expansion on a slowly growing symmetric window and the Gaussian tails give a uniform remainder in \(\phi\). On the two phase axes \(\sin\phi=0\), the remainder improves to \(O(R^{-4})\) since the paired cosine term is even in \(1/R\).

For \(d=1,3\), take \(A=1\) for each and \(B_1=1/2\), \(B_3=-e^{-1/3}/2\). The leading oscillation cancels and the remaining expansion is
\[
G_1(t)+G_3(t)=2+\frac{19e^{-1/2}}{36R}\sin\phi
+\frac{4/3+(313e^{-1/2}/1296)\cos\phi}{R^2}
+O(R^{-3}).
\]
More generally the leading cancellation cannot also cancel the \(R^{-1}\sin\phi\) coefficient for two distinct finite marks and nonzero amplitudes, since \(V\mapsto(5V+V^2)/4\) is strictly increasing. The conclusion concerns this fixed form of the center densities.

## Numerical stress and timing

The Metal kernel evaluates both \(d=1,3\) per profile. Each phase class uses the six dyadic meshes \(n=16,32,64,128,256,512\) on \([-12,12]\). It samples four fixed phases and random intermediate angles, with \(R=\pi k+\phi/2\) and \(k=2\) through 65536. The paired implementation integrates the correction to the leading Gaussian term before adding that term analytically. It exploits the exact even-weight increment \(4h^2\), odd weight, and \(\cos\eta-1=-2\sin^2(\eta/2)\). This retains residuals of order \(R^{-2}\) on the axes in Float32. The direct implementation subtracts 2 from the computed sum.

Among 60 sampled profiles replayed in 68-digit arithmetic, the largest pair errors were:

| Dyadic points | Stabilized pair error | Direct pair error |
| ---: | ---: | ---: |
| 16 | 0.0105875 | 0.174382 |
| 32 | \(1.2323\times10^{-5}\) | \(5.7936\times10^{-5}\) |
| 64 | \(7.8699\times10^{-9}\) | \(2.3838\times10^{-7}\) |
| 128 | \(1.5320\times10^{-8}\) | \(4.9249\times10^{-7}\) |
| 256 | \(7.0313\times10^{-9}\) | \(5.6251\times10^{-7}\) |
| 512 | \(2.1932\times10^{-8}\) | \(4.3075\times10^{-7}\) |

At \(\phi=0,k=65536\), the true pair residual is \(3.4909907000200\times10^{-11}\). The 128-point stabilized GPU value is \(3.490990\times10^{-11}\); direct Float32 gives \(-1.1920929\times10^{-7}\). For the same angle and \(k\), the 16-point stabilized value differs by \(2.854\times10^{-12}\), while the 64-point value differs by \(2.173\times10^{-18}\). This separates coarse-mesh loss from cancellation in direct Float32 accumulation. At the intermediate phase \(\phi=1.15,k=65536\), the pair residual is \(1.41919272098834\times10^{-6}\) and the error remaining after both analytic corrections is \(2.30729\times10^{-17}\). At \(\phi=0\) the corresponding remainder is \(1.09251\times10^{-22}\), consistent with the parity improvement.

The two CPU coordinate formulas agree within \(10^{-61}\) on a wider symmetric Gaussian window. The SciPy adaptive versus mpmath maximum absolute difference in the saved selected one-well results is \(2.518\times10^{-16}\). The local and remote numerical payloads match exactly. These are numerical checks of the frozen-mark transfer; the expansion and graph-comparison rate follow from the analytic calculation.

The M5 Pro production completed 10,538,188,800 paired profiles (21,076,377,600 individual well evaluations) in 10,050 command buffers. The interval was 09:12:16.937–09:13:02.946 UTC on 2026-11-18. Native Metal time was 31.632161587 s across 46.008758000 s wall. IOAccelerator utilization is sampled device-wide; its full histogram appears in `summary.json`. Both smoke and main returned zero, both PIDs were absent after collection, and the local copies match all 13 remote SHA-256 hashes.

### Reproduce

```sh
cache="$HOME/CodexRemoteMath/cache"
mkdir -p "$cache/swift-modules" "$cache/cc144-temp"
TMPDIR="$cache/cc144-temp" swiftc -O -module-cache-path "$cache/swift-modules" scan.swift -o wells_scan
sh run.sh main 46 1048576 72973211
"$HOME/.venvs/math/bin/python" check_cpu.py main.json main.cpu.remote.json
```

The wall cutoff makes the number of command batches depend on device throughput. The archived `main.json` reproduces the CPU audit.
