# Round 75 MacBook Cantor-distance scan and exact low-mode cancellation

This run reuses the validated native Cantor gap kernel in an
independent directory. It continues the finite-union model
(a(2\pi x)=h,\mathrm{dist}(x,K_m)^\alpha), with retained digit
set (D) in base (B) containing the two extreme digits. For
(n=1,\dots,6) the kernel estimates
[
 C_n(K_m)=\int_0^1 e^{2\pi inx-n^2a(2\pi x)/2}\,dx.
]
It ranks the first three and first six modes by the maximum
magnitude divided by each mode's lost-weight mass
(int(1-e^{-n^2a/2})). Complex raw values at the CPU-replayed
minima are recorded separately. The seed starts at 700001, beyond
the earlier seed range, and the midpoint rule uses 144 panels
per half-gap. The profiles are drawn from seven bases, five
exponents, and three interior-digit probabilities, giving 105 strata.

## Exact phase-sensitive tail formula

Let (
ho=|D|/B), let (d(y)=\mathrm{dist}(y,K_\infty)) on
the unit circle, and define
[
\sigma_m(n)=\prod_{j=1}^m\frac1{|D|}\sum_{d\in D}
             e^{2\pi i n d/B^j},\qquad
g_m=\frac{n^2h}{2}B^{-m\alpha},\qquad
\omega_m=2\pi nB^{-m}.
]
All removed-gap endpoints belong to the infinite set. Within
a retained level-(m) cylinder the distance is the scaled copy
(B^{-m}d(y)), because its own endpoints already belong to the
set. Hence
[
 C_n(K_\infty)-C_n(K_m)=
 -\rho^m\sigma_m(n)\int_0^1 e^{i\omega_m y}
   (1-e^{-g_m d(y)^\alpha})\,dy .
]
If (g_j) are the lengths of the first-level gaps, the exact
radial moment is
[
 M_\beta=\int_0^1 d(y)^\beta\,dy
 =\frac{\sum_j g_j^{\beta+1}}
 {2^\beta(\beta+1)(1-|D|/B^{\beta+1})}.
]
Therefore the tail is bounded by
(\rho^m|\sigma_m|g_m M_\alpha), and the error when replacing it
by (-\rho^m\sigma_m g_m M_\alpha) is at most
[
 \rho^m|\sigma_m|g_m
 (|\omega_m|M_\alpha+g_m M_{2\alpha}/2).
]
These bounds follow from the linear exponential loss bound and
the quadratic error of its linearization. If (1\le n<B^m/2)
the integral has strictly positive projection along the center
of its phase arc. In that range the exact omitted coefficient
vanishes precisely when the prefix factor vanishes.

## A continuous full-period Cantor profile that passes prescribed low modes

Given an integer (N\ge1), take distinct primes (p,q>N),
set (B=5pq), and omit the digits
[
 E=\{1+5qj:0\le j<p\}\cup
     \{3+5pk:0\le k<q\}.
]
All other digits are retained. The two extreme digits are retained
and no two omitted digits are adjacent, since their residues
modulo five are 1 or 3. Thus the largest gaps are exactly the
omitted first-level cells, all of width (1/B).

For (1\le n\le N), the omitted-digit phase sum (S_E(n))
vanishes as the sum of two complete prime orbits. All primary-gap
contributions then cancel because the equal gap shapes differ
only by phase. The full base-digit sum at these frequencies is
zero, so the retained digit factor also vanishes. The exact tail
identity therefore gives (C_n(K_\infty)=0) at all the prescribed
low frequencies, for every (h,\alpha>0).

The indicator of the union of largest gaps has nonzero Fourier
coefficients at (n=p,q), where one prime orbit survives at
each frequency. Their coprimality forces every translation which
preserves the largest-gap set to be an integer. A translation
period of the continuous distance profile preserves its zero
set and these largest gaps. Hence the distance profile has
fundamental period one. This is a finite-test construction.

For the explicit example (p=7,q=11,B=385), one can also prove
that its actual Gaussian output modes 7 and 11 are nonzero when
(alpha\ge1) and (0<h\le1). Denote by (I_\alpha) the
moment of the distance power on one primary gap and divide
the primary and tail magnitudes by
(|S_E(n)|(n^2h/2)I_\alpha). The primary is at least
[
 (1-8n^2/B^2)(1-n^2/(8B)),
]
using the cosine lower bound and the linear loss lower bound at
distances at most (1/(2B)). The exact radial tail bound is at
most (|E|/(B^{\alpha+1}-|D|)\le3/24643).
The rational margins checked by `continuous_control.py` are
strictly positive at both frequencies. That script's 65-digit
primary-gap checks and rational outputs match across both Macs.

## Scan and execution receipts

MacBook Metal PID 22042 ran from 04:01:18 to 04:03:55 UTC on
2026-11-18, exited zero, and was absent from the post-run process
query. Its native command buffers accumulated **120.476727
seconds**. The run completed **47,308,800** finite profiles
across 220 full cycles. Each of 252 sampled GPU-device values
was 100%, and the prescheduling device reading was 100%.

Both Macs reproduced the 210 retained candidates, including
24 direct enumerations of the finite gap lists. The maximum Metal
versus Float64 discrepancy at the same midpoint mesh was
(5.10\times10^{-7}). The maximum midpoint-versus-refined
quadrature discrepancy was (1.26731\times10^{-4}).
The direct positive 58-digit check of the grid's minimum uses
base 19, (D=\{0,8,10,18\}), depth 10, (\alpha=1/4) and
(h=29.439157). Its second raw mode has magnitude 0.00099539546
and its positive-weight normalization is 0.84552082. The refined
phase-sensitive truncation bound there is (2.42744\times10^{-9}),
with a leading-tail approximation error of (2.30452\times10^{-11}).

There is a conservative four-workload shared window from
04:02:41.431400 to 04:03:01.002 UTC, or **19.570600 seconds**,
using the observed post-run MacBook clock offset range
([-0.036053,+0.219373]) s if stable across the run.
The local GPU ended before the longer MacBook jobs.
Twenty-nine clock-contained MacBook monitor snapshots captured
the GPU and relevant CPU PIDs. The samples included all 15
CPU workers and all 29 corresponding GPU-device readings were 100%.

MacBook cycle wall times are monotonic relative to a start
captured after initialization. Its elapsed-zero local-clock time
lies between 04:01:17.779627 and 04:01:19.371535. The lower bound
uses the lifecycle start marker and clock margin. The upper bound
uses the result file's exact modification time after assembly of
the recorded elapsed duration, and the clock margin.
This places the completed-command checkpoint pair at cycles
120 and 144 inside the common window for every admissible anchor.
The command difference is 13.1192165 seconds. Enclosing local
GPU checkpoints bound local unscheduled time by 3.020636 seconds,
including rounding. Thus the two GPUs simultaneously executed
native commands for at least **10.09 seconds** while both CPU
worker windows were active.

## Reproduce

In the saved MacBook directory
`~/CodexRemoteMath/cc75-metal-nov18`:

```sh
export TMPDIR="$HOME/CodexRemoteMath/cache/cc75-temp"
xcrun swiftc -O -module-cache-path "$HOME/CodexRemoteMath/cache/swift" \
  -Xcc -fmodules-cache-path="$HOME/CodexRemoteMath/cache/clang" scan.swift -o cantor_scan
sh run.sh sync 144 120 2048 700001
"$HOME/.venvs/math/bin/python" check_cpu.py sync.json sync.cpu.json
"$HOME/.venvs/math/bin/python" continuous_control.py continuous.remote.json
```

The production run reused the already compiled Apple M5 Pro
binary. The host environment was macOS 27.3 build 26C5035r,
Swift 6.4.2, Python 3.12.13, SciPy 1.18.1, and mpmath 1.3.0.
The fixed seed reproduces the stream; the time-bounded cycle
count can vary. Selected sources, statistics, and receipt hashes
are recorded separately.
