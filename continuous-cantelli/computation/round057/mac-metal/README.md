# MacBook extended signed-window scans

The valley and outer compensation band are centered at (a).
Their half-supports are (0le yle w/a) and (rle yle r+L),
with (w/a<r). Let (P_V(s)) and (P_B(s)) denote the
Gaussian probabilities of their full symmetric supports at mean
(a+s). Choose (0<b<1), (A=bP_B(0)/P_V(0)), and let
(q=(1+s/a)^2). The signed normalized transform difference is
[
Delta(s)=ig[(1+A)^q-1ig]P_V(s)
  -ig[1-(1-b)^qig]P_B(s).
]
This is a specific local piecewise-constant variance perturbation
with the rest of the graph kept at one. It makes (Delta(0)=0).

The first sweep extends the number of profile parameters, offset
points, and midpoint panels from the earlier calculation. The next
two sweeps use exact finite-(a) derivative normalization
[
D(a)=(1+1/A)log(1+A)+(1-b)log(1-b)/b>0,qquad
s=
rac{
u D(a)}{a m_2},quad
m_2=mathbb E[Y^2\mid rle |Y|le r+L].
]
The positivity follows because the valley amplitude has normalized
derivative greater than one at (q=1), and the band amplitude
has normalized derivative smaller than one.

## Analytic check of the derivative scale

Fix (w,r,L,b) in their nondegenerate ranges and let (a) tend
to infinity. Since (Aasymp a), (D(a)simlog a). On a
bounded set of (
u), expanding the two amplitude ratios about
(q=1) gives a difference
(D(a)(q-1)+O((log a)^2(q-1)^2)).
The calibrated probability ratios follow from
(e^{-s^2/2}mathbb E[cosh(sY)\mid V]) and its band counterpart.
Their difference is (-m_2s^2/2+O(a^{-2}s^2+s^4)).
As (q-1=2s/a+s^2/a^2), the normalization by
(a^2m_2/D(a)^2) gives
[

rac{a^2m_2}{D(a)^2}
rac{Delta(s)}{bP_B(0)}
   =2
u-
rac{
u^2}{2}
      +O!left(
rac{(log a)^2}{a^2}
ight).
]
The error bound is uniform on bounded (
u)-intervals for the fixed
profile parameters. In particular there is a positive-to-negative
zero at (
u=4+O((log a)^2/a^2)). This follows by comparing the
signs at (4) plus and minus a sufficiently large multiple of the
error estimate.

The profile also produces **eventual positive** signed differences
at fixed (a) as (t=a+s) tends to positive infinity. The band loss
amplitude is at most one. The valley gain contributes
(exp(t^2log(1+A)/a^2)) up to a constant factor. The ratio of
the Gaussian interval probabilities has logarithm bounded below
by (-O(t)) using a fixed interior subinterval of the valley and
an upper bound on the band. The quadratic positive gain therefore
dominates. Thus for large (a) at a fixed profile shape, positivity
immediately right of (a), a negative intermediate value from the
scaled limit, and eventual positivity force at least two right
zeros. This argument describes the synthetic model.

## Grids and numerical recovery

All integrations use 2560 midpoint panels per positive half-support.
The nominal full-grid work is listed separately from targeted repairs.

| Run | Centers (a) | Shift scale and span | Evaluated or targeted finite values | GPU command seconds |
| --- | --- | --- | ---: | ---: |
| extension | 32 to 1024 | (log a,
uin[-2,8]) | 352,112,313 | 30.9250 |
| derivative | 64 to 65536 | (D(a),
uin[-2,8]) | 355,335,309 | 27.8205 |
| finite | 8 to 64 | (D(a),
uin[-2,16]) | 355,113,377 of 355,335,309 | 27.7949 |
| far repair, finite | 8 to 64 | same mesh, (|s|>4) | 36,544,180 | 12.1087 |
| far repair, large | 64 to 65536 | same mesh, (|s|>4) | 2,047 | 0.6247 |

The original finite sweep produced 221,932 nonfinite values at large
physical shifts. Separate Gaussian attenuation, cosh growth, and
amplitude growth can underflow or overflow before their products
are formed. The repair combines the amplitude and the Gaussian
action before exponentiating and uses a common exponent shift.
For physical results below Float32 range it stores a finite
rescaled sign, so those repaired magnitudes should be read solely
as signs. Both repair runs reported zero unexpected nonfinite
values. They skip the already computed small-shift part of the
corresponding mesh. Original raw outputs are included to show the
issue. Their high-shift crossing counts are superseded by the
targeted repair checks and are kept as raw execution data.

The MacBook and local CPU use independent stable Float64
quadrature or log-CDF differences, plus 80-digit mpmath roots.
Twelve recorded brackets resolved to an actual crossing in the
expanded CPU bracket. Four original large-(|s|) crossing brackets
were rejected by the CPU reference. All four occurred in the
targeted repair range. The repaired finite grid yielded four
selected CPU-confirmed crossings, with derivative-scaled values
between 5.232019948 and 15.992001952. The repaired high-(a)
far-shift subset had no crossing. On the original log-scaled
extension mesh, two selected roots at (a=1024) were
0.650509204 and 5.168919369, changing to 4.000005669 and
4.000516477 under the exact derivative normalization.
Two selected checks at (a=65536) gave 4.000001861 and
4.000001311 on the derivative scale, with 80-digit paths
confirming the values. At one grid edge the original Float32
literal bracket ended at 4 and CPU places the root just beyond
that edge. Such selected values do not locate global mesh extrema.

## Execution and reproduction

The Apple M5 Pro ran five bounded jobs, each exiting zero.
Aggregate native GPU command time was **99.274 seconds**.
The jobs occupied 114 seconds of total wall time with gaps for
new-code staging and diagnosis. This aggregate represents four
phases of computation and one further targeted recovery across
the saved run timeline. A post-run query for all five Metal PIDs
found none still running. Utilization is a device-wide metric
and the initial reading was 99%. Per-run metrics and timestamp
bounds are in `summary.json`.

Use the corresponding sources in the root and subdirectories.
Each directory is independent so its Swift program loads the
correct same-directory Metal kernel. For each directory compile on
the MacBook with the same cache paths:

```sh
export TMPDIR="$HOME/CodexRemoteMath/cache/cc52-temp"
xcrun swiftc -O -module-cache-path "$HOME/CodexRemoteMath/cache/swift" \
  -Xcc -fmodules-cache-path="$HOME/CodexRemoteMath/cache/clang" scan.swift -o window_scan
```

Run the following respective commands in the root, `derivative`,
`finite`, `far_repair`, and `large_far_repair` directories:

```sh
sh run.sh extension 37 23 21 19 17 61 2560 8192
sh run.sh derivative 37 21 21 21 17 61 2560 8192
sh run.sh finite 37 21 21 21 17 61 2560 8192
sh run.sh far_repair 37 21 21 21 17 61 2560 8192
sh run.sh large_far_repair 37 21 21 21 17 61 2560 8192
```

From the root, execute
`"$HOME/.venvs/math/bin/python" check_cpu.py . cpu.remote.json`.
The MacBook environment was macOS 27.3 build 26C5035r, Swift 6.4.2,
Python 3.12.13, SciPy 1.18.1, and mpmath 1.3.0. The search uses
no random seed. Selected sources and outputs match the hashes
in `remote.sha256`.
