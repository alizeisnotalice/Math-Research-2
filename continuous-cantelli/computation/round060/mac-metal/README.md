# MacBook far-right return window

The input has a symmetric valley of halfwidth (v=w/a) at (a)
and a compensating band (rle |x-a|le r+L). Write (R=r+L).
With Gaussian mean (a+s), let (P_V(s)) and (P_B(s)) denote
the probabilities of the valley and band. Calibrate
(A=bP_B(0)/P_V(0)), (ell=log(1+A)>0), (0<b<1).
For (q=(1+s/a)^2), the signed excess is
[
Delta(s)=(e^{qell}-1)P_V(s)-[1-(1-b)^q]P_B(s).
]
This model has a centered zero and, for sufficiently large (a),
a near-right negative excursion. The present calculation resolves
the far-right return on the scale (s=sigma R a^2/ell).

## Endpoint expansion and return location

For an interval ([u-d,u]), after extracting the common
(e^{-s^2/2}/sqrt{2pi}) factor, its Gaussian probability is
[
 K(s;u,d)=e^{su-u^2/2}int_0^d e^{-(s-u)h-h^2/2},dh.
]
For (s>2|u|), the latter integral is
[

rac{1+O(s^{-2}+e^{-(s-u)d})}{s-u}.
]
This follows by bounding (|e^{-h^2/2}-1|le h^2/2) and completing
the exponential integral to infinity. For the far-right window
near (sigma=1), the valley has (u=v, d=2v), and (sv) tends
to infinity. The positive side of the outer band has (u=R,d=L);
the negative side is exponentially smaller. Both amplitude
factors have exponentially small relative errors here.

Let (H(s)) be the log of the positive-to-negative mass ratio in
(Delta(s)) and set
[
 k=R-v-2ell/a,qquad C=ell+(R^2-v^2)/2,qquad
 Q(s)=ell s^2/a^2-ks+C.
]
Then in a fixed sufficiently narrow band around (sigma=1),
[
 H(s)=Q(s)+log
rac{s-R}{s-v}
       +O(s^{-2}+e^{-c s v}).
]
For fixed (w,r,L,b) and large (a), the larger root of (Q) is
[
 s_Q=
rac{a^2}{2ell}left(k+sqrt{k^2-4ell C/a^2}
ight).
]
It lies in the stated band. The derivative of each log probability,
after removing the common quadratic factor, is the corresponding
Gaussian tilted conditional mean of (x-a). Endpoint concentration
gives the means (v+O(1/s)) and (R+O(1/s)). Hence
(H'(s)) is bounded below by a positive constant in the chosen
band and tends to (R) at (s_Q+o(1)). Sign comparison then gives
a unique zero (s_*) in that band and
[
 (s_*-s_Q)s_Qlongrightarrow1.
]
This identifies a specific far return root in the given model.

## Numerical method

The grid uses 47 logarithmic (ain[16,65536]), 29 logarithmic
(win[1/32,1.5]), 27 logarithmic (rin[1/4,2.5]), 23 logarithmic
(Lin[1/8,3]), and 19 uniform (bin[0.02,0.98]).
For each core the kernel samples 73 factors (sigmain[0.5,1.5]).
It computes the signed log odds divided by (sR), removing the
common large Gaussian exponent before evaluation. Each positive
half-support is integrated in the endpoint coordinate
(z=s(u-y)), with 3072 midpoint panels and an upper cap of 48.
The negative half-support is included by its exact multiplying
factor (1+exp(-2sy)). On this grid (s>2R), so the omitted
positive endpoint integral tail is controlled by the
(exp(-z/2)) envelope beyond the cap.

The independent Float64 check computes the exact interval
probabilities from scaled complementary error functions.
The 95-digit check uses separately rescaled Gaussian tails and
solves the crossing without Float32. For the fixed profile
(w=1/2,r=1,L=1,b=1/2), the computed
((s_*-s_Q)s_Q) values at (a=32,128,1024,65536,10^8)
are approximately 1.093144, 1.031503, 1.005780, 1.000153,
and 1.000000173. These are numerical checks of the analytic
limit above.

## Main outcome and checked runtime

The GPU evaluated **1,173,988,701** finite signed-log values over
16,082,037 profile cores. There were 15,973,670 sampled
negative-to-positive crossing brackets inside (sigmain[0.5,1.5]).
All final offset scores were positive. The other 108,367 profiles
started with a nonnegative score; these lack a crossing bracket on
the stated mesh. Selected brackets were replayed on both machines
with Float64 and 95-digit integration results agreeing. The
earliest selected factor was about 0.507075964. At (a=65536)
the selected factor was 0.9993056214. The large-center result
identifies a sampled profile and does not establish a global
parameter extremum.

The Apple M5 Pro process PID 71769 started at **03:04:09 UTC**
and exited zero at the 03:08:36 UTC lifecycle receipt on
2026-11-18. The post-run process query verified it was absent.
It accumulated **257.278997 seconds** of native GPU command time.
The initial device-wide GPU reading was 100%. The 442 production
samples had mean 99.321%, minimum 0%, and maximum 100%.
The GPU computation completed its bounded full grid.

The conservative four-device interval from the first shared run
is **03:04:22.329480–03:05:04.326 UTC**, totaling 41.996520
seconds. The MacBook CPU timer interval is first translated using
the post-run remote-minus-local offset bounds
([-0.024974,+0.116470]) s, if stable during the run. The local
CPU full-worker activity covers that interval. The local Metal
cap uses its last completed checkpoint before the first wrapper
stopped. In the interval, 67 MacBook monitor snapshots captured
both the Metal PID and relevant remote CPU workers, with a peak
count of all 15 workers and aggregate sampled CPU use of 826.6%.
Its 67 corresponding device GPU readings were 100%.

Both GPU queues submitted sequential commands. The integer-second
MacBook GPU lifecycle markers bound its whole wall time above by
268 seconds. Subtracting the measured command total bounds
potential unscheduled MacBook time above by 10.721004 seconds.
The local Metal checkpoint pair at 03:04:22.219 and
03:05:04.326 encloses the conservative common interval. Its
cumulative command values, 16.564656 and 57.061458 seconds,
bound local unscheduled GPU time in that enclosure above by
1.612198 seconds including timestamp granularity. Subtracting
these two upper bounds from the common interval gives at least
**29.66 seconds** of simultaneous native GPU command activity
while both CPU group windows were active. The corrected local
GPU follow-up was a separate run.

## Reproduce

From the authenticated MacBook directory
`~/CodexRemoteMath/cc58-metal-nov18` with the saved sources:

```sh
export TMPDIR="$HOME/CodexRemoteMath/cache/cc58-temp"
xcrun swiftc -O -module-cache-path "$HOME/CodexRemoteMath/cache/swift" \
  -Xcc -fmodules-cache-path="$HOME/CodexRemoteMath/cache/clang" scan.swift -o window_scan
sh run.sh sync 47 29 27 23 19 73 3072 8192
"$HOME/.venvs/math/bin/python" check_cpu.py sync.json sync.cpu.json
```

The MacBook environment was macOS 27.3 build 26C5035r,
Swift 6.4.2, Python 3.12.13, SciPy 1.18.1 and mpmath 1.3.0.
No random seed is needed. Selected sources, outputs, and recording
files match their remote SHA-256 hashes.
