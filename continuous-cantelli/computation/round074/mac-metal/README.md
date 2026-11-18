# MacBook low-frequency Haar tests for digit Cantor profiles

Write (x=\theta/(2\pi)) and use the unit-circle distance.
For base (B), a digit set (D\subset\{0,\dots,B-1\})
containing (0,B-1), and depth (m), let (K_m) be the union
of closed retained digit cylinders. The scan uses
[
a(2\pi x)=h\operatorname{dist}(x,K_m)^\alpha,qquad
C_n=\int_0^1 e^{2\pi inx-n^2a(2\pi x)/2}\,dx
]
for (1\le n\le6). The objective is
(max_{n\le3}) or (max_{n\le6}) of
(|C_n|/\int_0^1[1-e^{-n^2a/2}]dx).
It measures the relative low-frequency imbalance of the lost
weight, and all raw values are retained for the selected profiles.
The GPU grid is an empirical ranking of sampled candidates.

## Factorization and infinite-level error

Because the uniform-input nonzero Fourier modes vanish, each
(C_n) is the negative sum of the Fourier contributions of
(1-e^{-n^2a/2}) on the complementary gaps.
At level (k) a gap inside a cylinder is a scaled copy of an
initial gap, and the prefix phase factor divided by the number
of cylinders is
[
 \prod_{j=1}^{k}\frac1{|D|}\sum_{d\in D}e^{2\pi i n d/B^j}.
]
The combined cylinder mass factor is ((|D|/B)^k).
Each original gap is split at its center, and the symmetric
weight requires only one half-gap integral. The Metal kernel
uses a 48-panel midpoint rule on that half.

The endpoints of all removed gaps remain in (K_\infty), so
the distances to (K_m) and (K_\infty) coincide outside
(K_m). The latter distance on (K_m) is at most (B^{-m}/2).
For each (n), the exact finite-to-infinite coefficient error obeys
[
 |C_n(K_m)-C_n(K_\infty)|\le
 (|D|/B)^m\frac{n^2 h}{2}(B^{-m}/2)^\alpha.
]
This controls truncation separately from numerical quadrature.

## Numerical scan and checks

Seven bases (3,4,5,7,9,13,19), five exponents (1/4,1/2,1,2,3),
and three interior-digit retention probabilities (1/4,1/2,3/4)
give 105 strata. The two extreme digits always appear. A full
draw has its middle digit removed. The pseudo-random height is
log-uniform on ([0.03,30]), and the depth is uniform on
(1,\dots,10). Some strata have only a few distinct masks,
particularly base three, while the height and depth vary.

The run completed **98,703,360** finite profiles, covering
six harmonics per profile. It took 459 full strata cycles and
940,032 profiles per stratum. The two minima retained per
stratum, for the first three and first six modes, were checked
on both Macs. The transformed 128- and 256-node Gauss–Legendre
rules agreed to at most (6.12\times10^{-14}) in the checked
integrals. Explicit enumeration of 24 finite-level gap lists
agreed to at most (3.55\times10^{-14}). Metal versus Float64
using its same midpoint mesh differed by at most
(5.81\times10^{-7}) in the reported scores. The largest
midpoint-versus-refined score difference was
(3.348\times10^{-4}), so quadrature resolution dominates those
cases. The two CPU calculations gave identical numerical
verification outputs.

The sampled minimum for both objectives was the base-19
mask (D=\{0,8,10,18\}), (m=10), (\alpha=1/4),
and (h=29.682938). Its Metal loss-normalized score was
0.00095333514. A separate 58-digit CPU calculation integrated
the positive Gaussian weight directly over the gaps and retained
cylinders. It found (|C_2|\approx0.00097823335987108)
and (|C_2|/Z_2\approx0.8455290564), where
(Z_2=\int_0^1 e^{-2a}dx). The small raw harmonic reflects damping while the
positive-weight normalized mode remains large.
The finite-to-infinite bound at that frequency is
(5.426\times10^{-9}). The high-precision check is in
`positive.local.json`.

## Exact controls

**Middle thirds distance.** Let (B=3,D=\{0,2\}), and
(0<h\le1,\alpha\ge3/4). On the primary gap
([1/3,2/3]), (cos(2\pi x)\le-1/2), and
[
1-e^{-h d^\alpha/2}\ge(7/8)(h/2)d^\alpha,
\qquad d=\operatorname{dist}(x,K).
]
Indeed (d^\alpha\le\sqrt{1/6}<1/2).
If (I) is the integral of (d^\alpha) over the primary gap,
the sum of these moments over the other gaps is at most
(I\rho/(1-\rho)), where (\rho=2/3^{\alpha+1}).
The rational inequality (3^7>(46/7)^4) gives
(\rho/(1-\rho)<7/16). The maximum absolute contribution of
all other gaps is therefore strictly smaller than the magnitude of
the primary gap's negative first lost-weight harmonic. Thus
(mathrm{Re},C_1>0) for this family, at every finite
depth and at the infinite Cantor limit.

**Two-level digital masks.** Give a finite cylinder union
variance zero and its complement any fixed positive variance.
For (n\ne0), its (C_n) equals a nonzero variance contrast
times
[
 \frac{e^{2\pi i n/B^m}-1}{2\pi i n}
 \prod_{j=1}^{m}\sum_{d\in D}e^{2\pi i n d/B^j}.
]
For any finite number (N) of low modes, choose distinct
primes (p,q>N), (B=2pq), and
(D=\{2qj:0\le j<p\}\cup\{1+2pk:0\le k<q\}).
The two disjoint digit orbits make the first digit factor
zero for every (n=1,\dots,N). At (n=p,q) exactly one
orbit survives and the later digit sums are confined to open
semicircles. The integral factor is also nonzero at these
two frequencies. Since (gcd(p,q)=1), the step profile has
fundamental period one even though it passes the selected low
frequency tests. The exact guard and 85-digit example checks
use (p=7,q=11,B=154). This measurable step control is a
separate model from the continuous distance profile.

## Run and reproducer

The MacBook Apple M5 Pro production PID 9160 started at
03:46:53 UTC on 2026-11-18 and exited zero at 03:49:37.
Native command buffers totaled **100.1559445 seconds**.
Across 273 device-wide sampled values, GPU use averaged
99.634%; the initial device reading was 100%.
Both pilot and production PIDs were absent after completion.

From the authenticated MacBook directory
`~/CodexRemoteMath/cc68-metal-nov18`:

```sh
export TMPDIR="$HOME/CodexRemoteMath/cache/cc68-temp"
xcrun swiftc -O -module-cache-path "$HOME/CodexRemoteMath/cache/swift" \
  -Xcc -fmodules-cache-path="$HOME/CodexRemoteMath/cache/clang" scan.swift -o cantor_scan
sh run.sh smoke 24 0 32 613
sh run.sh sync 48 100 2048 3751
"$HOME/.venvs/math/bin/python" check_cpu.py sync.json sync.cpu.json
"$HOME/.venvs/math/bin/python" exact_controls.py exact.remote.json
```

The seed fixes the candidate stream and the runtime stops at
the first full cycle reaching the GPU target. The MacBook
environment was macOS 27.3 build 26C5035r, Swift 6.4.2,
Python 3.12.13, SciPy 1.18.1 and mpmath 1.3.0.
The direct positive reference was run locally with
`positive_reference.py sync.json positive.local.json`.
