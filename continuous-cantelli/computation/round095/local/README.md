# Early reflected Stefan-wall stress test

This task tests the positive reflected Gaussian cap from the old-wall comparison. It is independent of the MacBook scan of a moving Laplace saddle and far-tail coefficients.

For fixed old variance $\tau$, wall $b$, new heat duration $h$, and exterior distance $d$, set $\sigma^2=\tau h/(\tau+h)$, $k=d\sigma/h$, and $v=b\sigma/\tau$. The exterior heat of the reflected cap equals $2v\sigma p_h(d)p_\tau(b)R(k,v)$, where
$$
R(k,v)=\int_0^\infty we^{-kw-w^2/2}\frac{\sinh(vw)}{vw}\,dw.
$$
At $v=0$ this is the continuous normalized derivative limit. The underlying unnormalized cap vanishes there.

## Local Metal run

The deterministic screen covered $2^{-12}\le k\le256$, $0\le v\le64$, with distinct slices near $v=0$, near $k=v$, in the interior-saddle regime $k<v$, and in the far fixed-wall regime $k\ge16,\ v/k\le0.1$. The 64-node-per-panel Metal integrator uses the positive normalized integrand. It completes the square in the interior regime and accumulates the result in logarithms.

Apple M6, macOS 27.2 build 26B123, Apple Swift 6.4.2. Main PID 27432 ran from 05:25:01.313 to 05:26:41.315 UTC on 2026-11-18 and exited with code zero. Wall time 100.002626 s, summed native Metal command intervals 80.854477 s, submission and waits 93.311899 s. There were 4,028,366,848 parameter combinations with zero nonfinite outputs. During 185 systemwide GPU samples the range was 0%–94% and mean 86.778%. Before and after readings were 7% and 93%. The systemwide figure includes other work. Main PID 27432 and smoke PID 27196 were verified absent.

## Accuracy and wall geometry

Nine stratum representatives were checked independently using adaptive double quadrature and the exact scaled Gaussian-tail formula. Their maximum GPU discrepancy in $\log R$ was $8.12\times10^{-5}$, at an extreme $\log R$ near 2044. The two CPU methods agreed within $1.47\times10^{-10}$. Three representatives were also checked at 70-digit working precision. Direct quadrature of the unnormalized physical heat kernel at three separate input geometries agreed with the normalized formula to $8.11\times10^{-15}$ relative.

- For $k>v$, the analytic bound $R\le(k^2-v^2)^{-1}$ is exact, and the companion lower estimate in proof.md quantifies the finite error. The GPU reported no overshoots.
- On $k\ge16,\ v/k\le0.1$, the exact estimate gives $|k^2R-1|<0.013$. The maximal sampled absolute logarithmic deviation was about 0.01156. Fixed physical geometry gives $k^2R\to1$; a sequence with $d=\sqrt h,\ b=\tau=1$ instead gives $k^2R\to0.34432046$.
- A finite exact certificate gives $k^2R>6$ at $k=64,v=63$, while 70-digit computation gives approximately 21.0588368. This demonstrates the geometry restriction on extending the fixed-wall equivalent to a uniform finite upper bound.
- For fixed positive $\tau,b,d$, the expansion of the cap has relative factor $1-3h/d^2+h^2[b^2/(\tau^2d^2)-3/(\tau d^2)+15/d^4]+O(h^3)$ after removing its leading $h^{3/2}e^{-d^2/(2h)}$ and constant. The proof supplies the leading constant and the analytic remainder argument. The coefficients apply to this reflected cap.

The exact cutoff and small-argument truncation checks are in certificate.py. They bound nominal real-arithmetic tails below $10^{-18}$ in the endpoint regime and below $10^{-40}$ after normalization in the interior regime. Independently measured discrepancies include the effects of binary32 and quadrature.

## Reproduce

```sh
./run.sh 0 48 smoke
./run.sh 100 64 main
PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python verify.py wall-main.json cpu-main.json
/Users/zzh/.venvs/main/bin/python certificate.py
```
