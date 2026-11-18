# Periodic heat-kernel Fourier scan on Apple M6

The model pushes a uniform phase $\theta$ through Gaussian noise on the circle with phase-dependent variance
$$
v(\theta)=d\,|2\sin(g\theta/2)/g|^p[Q_0(\theta)+\eta B(\theta)].
$$
Here $Q_0$ has harmonics $g,2g,3g$, with both sine and cosine terms, and $B$ has harmonics 1, 2, and 3. Their coefficients are normalized so $Q_0\ge0.6$, $|B|\le1$, and $Q_0+\eta B\ge0.35$. The model's output coefficient is
$$
F_n=(2\pi)^{-1}\int e^{in\theta-n^2v(\theta)/2}\,d\theta.
$$
For an off-multiple $g\nmid n$, the unperturbed coefficient is exactly zero by symmetry. When $\eta>0$ the GPU stores the scaled leakage $F_n/\eta$ and its linear response at zero. Other stored amplitudes are raw Fourier coefficients. A positive common variance would add a known factor $e^{-n^2v_*/2}$.

## Profiles and numerical method

The primary run used deterministic counter hashes, $g\in\{1,2,3,4\}$, $n\in\{1,\ldots,96\}$, and valley exponents $p\in\{0.5,0.75,1,1.5,1.75,2,2.25,3,4,6,8\}$. It varied $d$ logarithmically over $[0.01,16]$. One third of the counter families used $\eta=0$, and the others used logarithmic $\eta\in[2^{-18},1/4]$.

The Float32 kernel integrates each valley cell in local coordinates. Panels resolve the valley scale $(2/(dn^2))^{1/p}$ and bound the oscillatory width in phase. The first panel uses a fourth-power substitution at the valley. The kernel used 72 Gauss–Legendre nodes per panel and an analytic tail cutoff. For scaled leakage it evaluates the heat-weight difference with a stable small-shift expm1 polynomial. The symmetry-zero cases return zero directly. The main count includes these structurally decided profiles.

## Actual computation and verification

Apple M6, macOS 27.2 build 26B123, Apple Swift 6.4.2. Main PID 79547 ran from 03:21:12.754 to 03:22:57.763 UTC on 2026-11-18 and exited with code zero. The 105.011284 s run processed 523,632,640 mode/profile combinations. Summed native Metal command-buffer intervals were 99.996543 s, and command submissions and waits totaled 103.322529 s. No values were nonfinite. There were 192 during-run systemwide GPU utilization readings, from 99% to 100%, averaging 99.990%. Before and after readings were 99% and 100%. These are device-wide readings. The project-specific times come from the Metal timestamps. Both main PID 79547 and smoke PID 77555 were verified absent.

| Stratum | Combinations |
|---|---:|
| All | 523,632,640 |
| Isolated subquadratic, finite-to-limit comparison | 25,459,426 |
| Isolated quadratic, finite-to-limit comparison | 5,090,622 |
| Isolated superquadratic | 45,252,720 |
| Symmetric, forbidden mode | 83,632,828 |
| Broken symmetry, forbidden mode | 167,280,145 |
| Weak breaking, linear response | 84,010,263 |
| Multiple valleys, allowed mode | 141,811,505 |
| Isolated sharp, high mode | 12,270,692 |
| Isolated smooth, high mode | 12,271,954 |

The strata overlap. Their saved representatives maximize the stratum's score. For the isolated subquadratic and quadratic asymptotic strata, the score is the relative difference from the appropriate leading expression. For weak symmetry breaking it is $|F_n/\eta-L_n|$. In the other strata it is the magnitude of the stored amplitude.

Independent SciPy integration of the ten saved representatives agreed with Metal by at most $5.06\times10^{-8}$ in absolute complex amplitude. For linear response the largest discrepancy was $6.06\times10^{-9}$. Independent 65-digit integration additionally checked the subquadratic comparison and weak-breaking representatives. At the weak-breaking representative, $\eta\approx9.98\times10^{-4}$ and the GPU observed $|F_n/\eta-L_n|\approx2.264\times10^{-5}$.

On five fixed profiles, CPU controls extended to $n=385$, preserving residue modulo $g$. At $n=385$, an isolated $p=0.75$ profile's rescaled error from the predicted limit was about $2.64\times10^{-8}$. An isolated quadratic profile's error was about $2.70\times10^{-5}$. The rescaled $p=4$ coefficient had magnitude about $4.41\times10^{-8}$. The largest observed relative deviation for the subquadratic main stratum was about 0.695 at $p=1.75,n=24$. These finite-mode differences are compatible with the stated limits.

## Analytic conclusions and bounds

The local limits, exact forbidden-mode cancellation, and incompatibility of any studied profile with a uniform output are proved in proof.md. Below quadratic order the scaling $n^{2/p}F_n$ has an explicit valley-weighted limit on each residue class. At quadratic order the weights contain the Fourier factor of a Gaussian. Above quadratic order the scaled coefficient tends to zero by Riemann–Lebesgue and $L^1$ convergence. The latter assertion supplies no rate. The exclusion refers to this circle model.

The Fraction checks establish real-arithmetic omission bounds $10^{-16}$ for raw tails and $10^{-14}$ for scaled leakage and its derivative. They also bound the small-shift normalized polynomial error by $1.4\times10^{-10}$ and the linear-response remainder by $9\eta/4$. Independent quadrature treats numerical rounding separately.

## Reproduction

Run from this directory on Apple Silicon with the data volume mounted. Timed replays can stop at different batch counts.

```sh
./run.sh 0 48 smoke
./run.sh 105 72 main
PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python verify.py heat-main.json cpu-main.json --precise
/Users/zzh/.venvs/main/bin/python certificate.py
```
