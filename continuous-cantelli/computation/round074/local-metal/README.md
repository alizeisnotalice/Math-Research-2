# One-step circular Gaussian modes with a Cantor minimum set

On the unit circle let $C_\lambda$ be the Cantor set of the two endpoint similarities with ratio $\lambda$, and take variance $v(x)=d\,\mathrm{dist}(x,C_\lambda)^p$. The Fourier coefficient for uniform input is
$$
F_n=\int_0^1 e^{2\pi i n x-2\pi^2 n^2v(x)}\,dx.
$$
The Gaussian convention is variance in the unit-circle coordinate. Reflection makes $F_n$ real.

The exact level compression and limit proofs are in proof.md. A level-$j$ gap has length $l_j=(1-2\lambda)\lambda^{j-1}$ and thermal scale $s_j=2\pi^2 n^2d(l_j/2)^p$. All gaps of the same level are integrated together through their factored phase product. The correction subtracts the zero-variance Fourier coefficient before summing. The coarse contribution sums levels with $n l_j\ge1$ or $s_j\ge1$; the fine contribution sums the rest.

## Parameter grid and Metal run

The deterministic grid sampled $\lambda\in[0.23,0.46]$, with a separate exact ternary branch $\lambda=1/3$, and powers $p\in\{0.5,0.75,1,1.25,1.5,2,2.5,3,4,6,8\}$. It varied $d$ logarithmically over $[2^{-10},32]$ and $n\le1024$. Frequency families included powers $3^k$, their neighbors $3^k\pm1$, and hashed $n$. The exact ternary branch uses modular phase reduction and an exact integer comparison for the gap-length threshold when it can be active.

Float32 Metal evaluates each remaining integral by 48-node Gauss–Legendre quadrature per panel. The degree-four near-endpoint coordinate smooths the half-power contacts. The panel splits constrain oscillatory integration. For larger thermal exponents the code truncates the positive heat factor and subtracts a sinc; for small exponents it directly evaluates the stable exponential difference. The cutoff stops either at the tail budget or at depth 96.

| Recorded main quantity | Value |
|---|---:|
| UTC interval on 2026-11-18 | 03:47:15.825–03:49:00.826 |
| Wall time | 105.001688 s |
| Summed Metal command intervals | 78.990880 s |
| Submission and wait intervals | 103.842881 s |
| Mode/profile combinations | 395,870,208 |
| Nonfinite outputs | 0 |
| During-run GPU utilization observations | 195 |
| Systemwide GPU utilization minimum, mean, maximum | 0%, 99.487%, 100% |

The hardware was Apple M6 running macOS 27.2 build 26B123 and Apple Swift 6.4.2. Main PID 13444 exited with code zero. The main, smoke PID 12770, and paired-precision PID 23009 were verified absent. Device-wide utilization includes other processes. The Metal timestamps identify the work of this process.

## Numerical audit and limitations

Independent CPU adaptive quadrature checked the total, coarse and fine values of all ten saved stratum representatives. The largest absolute discrepancy was $8.65\times10^{-7}$, in the coarse value for a high quadratic mode at a nonternary ratio. For the deepest 96-level representative, the total was $2.981\times10^{-7}$ on Metal and $3.818\times10^{-7}$ on CPU. Its relative difference is about 22%, despite the small absolute difference.

A matched 48/80/112-node Metal replay of all ten inputs changed the total and component values by at most $2.99\times10^{-8}$ between adjacent orders. The largest discrepancy to CPU remained $8.65\times10^{-7}$. The discrepancy is not removed by increasing the quadrature order in this binary32 implementation. Tiny coefficients require a higher-precision check.

Two representatives were independently integrated with 65-digit working precision. In the strongest coarse/fine cancellation stratum the contributions were approximately $-0.1123946460$ and $+0.0656200746$, giving $-0.0467745714$. The selected $n=27$ ternary resonance gave $+0.0821167870$. The high-precision calculations carry a separately recorded analytic series tail bound. Two explicit tests with 255 gaps each agreed with the compressed finite-level formula to $1.55\times10^{-13}$ absolute.

The exact real-arithmetic Fraction calculation bounds the 96-level geometric tail by $6.35\times10^{-13}$. Combining the geometric tail, small-level leading approximation, cutoff of the positive heat term, expm1 polynomial and sinc polynomial gives a pre-quadrature omitted bound below $4\times10^{-12}$. These are analytic error budgets for the real-arithmetic rules; they do not bound GPU rounding or quadrature.

## Resonant shape

For standard ternary ratio, $p<2$, and $n=3^k+r$ with fixed residue shift $r$, the ratio of $F_n$ to the total mass of the same heat weight tends to $\widehat\mu_C(1)\widehat\mu_C(r)$. At $p=0.75,d=0.3,n=19683$ the CPU ratio was $0.37143734777$, compared with the limit $0.37143735671$. At $p=1.5$ and the same terminal frequency it was $0.37000782940$. These are finite-frequency checks. The proof gives eventual positivity along the resonant and adjacent frequencies for the specified subquadratic ternary model. The proof also gives the circle nonuniformity conclusion for all ratios in the quadratic and superquadratic cases.

## Reproduce

Run on Apple Silicon with the data volume mounted. Timed replays can stop at a different prefix length.

```sh
./run.sh 0 36 smoke
./run.sh 105 48 main
PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python verify.py cantor-main.json cpu-main.json --precise
PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python explicit_audit.py
/Users/zzh/.venvs/main/bin/python certificate.py
./compare/run.sh
```
