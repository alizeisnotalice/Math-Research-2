# Independent Metal segment and matched quadrature-order checks

This run continues the repaired paired-quartic model. On the two sides, compare $cu^4(1+b(u/r)^2\pm h u/r)$ with the symmetric baseline $cu^4(1+b(u/r)^2)$, at $y=a=kh/r$. The Gaussian-normalized core density difference divided by $h^2$ is the full scalar. The leading scalar is its second-order approximation evaluated with the same source normalization. The GPU counts signs as separated only if both scalar magnitudes are at least $2\times10^{-6}$. This is a floating-point reporting threshold.

## Actual execution

The stated 02:58:00 UTC gate had already passed when the task began locally. The first actual Metal attempt started at 03:04:04.975, PID 45448. Its launcher erroneously inserted the index argument into the sampling condition and then exited, while the kernel process continued. After attaching a separate sampler, the last completed checkpoint at 03:05:04.326 recorded 983,040,000 profiles and 57.061458 Metal command seconds. The process disappeared before writing a final JSON. The partial logs are retained as quartic-partial.log, ioreg-partial.txt, and monitor-partial.log. They do not supply a completed 80-second run.

After correcting the sampling condition, the replay ran **03:06:30.407–03:08:00.410 UTC** on 2026-11-18. It is a separate completed run with actual synchronization evidence bounded by these timestamps. It does not establish an 80-second overlap at the past gate. The replay retained the same new start counter as the partial attempt and used 128 nodes per panel.

| Quantity | Completed run |
|---|---:|
| Wall time | 90.002955 s |
| Summed native Metal command intervals | 86.528803 s |
| Command submission and waits | 87.256048 s |
| Parameter profiles | 1,487,405,056 |
| Counter interval | [2,147,483,648, 3,634,888,704) |
| Nonfinite profiles | 0 |
| IOKit during-run GPU samples | 168 |
| Systemwide GPU utilization minimum, average, maximum | 0%, 99.274%, 100% |

The counter interval ends before $2^{32}$; it is disjoint from the archived order-96 run's $[0,1{,}950{,}875{,}648)$. The launcher recorded PID 48168 exiting with code zero after all 2,837 batches. The before and after device samples were 100% and 0%, respectively. Systemwide utilization includes other work. Project attribution comes from command-buffer timestamps.

Independent CPU quadrature checked the eight maxima of the full-minus-leading gap, one for each reported stratum. The full scalar agreed to $5.51\times10^{-7}$ maximum absolute error. The gap maximum discrepancies were at most $8.69\times10^{-7}$. The raw GPU log contains 12,934,969 separated sign differences overall and 430 in the lowest-$h$ tuned stratum. These are numerical counts.

## Matched inputs and sensitivity

The paired check fixes the actual binary32 profile, including $h$ and $k$. It compares quadrature orders without refitting $k$ between them. Its 73 distinct saved inputs combine current and previous stratum maxima with the previous sparse/rare audit and the formerly unstable profile. Every paired input was independently integrated with SciPy double quadrature.

| Comparison | Maximum absolute difference over both scalars |
|---|---:|
| Order 96 versus CPU | $4.06\times10^{-7}$ |
| Order 128 versus CPU | $6.19\times10^{-7}$ |
| Order 160 versus CPU | $4.75\times10^{-7}$ |
| Order 96 versus 128 | $5.28\times10^{-7}$ |
| Order 128 versus 160 | $6.48\times10^{-7}$ |

All three orders agreed with the CPU on both signs for every one of the 55 paired inputs whose CPU magnitudes both met the separation threshold. The paired diagnostic used three Metal command intervals with combined duration 0.779 ms. Main PID 48168, the partial PID 45448, and paired PID 51509 were all verified absent. Native completion metadata are in the respective logs.

## Real-arithmetic bounds for the repair

Let $t=u/r$, $V=\varepsilon+cu^4(1+bt^2)$, $Z=u^2/V$, and $\ell=cu^4t/V$. The relative variance change on each side is $\pm x$, where $x=h\ell$. In the prescribed parameter domain $|x|<3/10$, $V^{-1/2}\le2^{27}$, $h\ge2^{-14}$, and $|a u|\le3$. The common source factor $C_h=2\sqrt{(M+\varepsilon)/(2\pi)}\exp(-a^2(1-1/(M+\varepsilon))/2)$ is below 3.

On the region $Z>192$, the corrected and baseline Gaussian terms are bounded by $2V^{-1/2}e^{-Z/(2(1.3))}$ and $V^{-1/2}e^{-Z/(2(1.3))}$. After both sides are combined and divided by $h^2$, the omitted full scalar is bounded in absolute value by
$$
3\cdot3\cdot27\cdot2^{27}\cdot2^{28}e^{-960/13}
< 10^{-13}.
$$
Here the positive-$u$ integration interval has length at most one. For the leading scalar the integrand is a linear combination of $1,Z,Z^2$ times $e^{-Z/2}$, with coefficients bounded using $\ell\le8/5$ and $|k|\le49/4$. Each $Z^j e^{-Z/2}$ decreases for $Z\ge192$ and $j\le2$. Substitution of the endpoint and a rational exponential upper bound give an omitted leading scalar below $10^{-27}$.

For $|z|\le1/4$, the degree-nine odd polynomial used in place of sinh has full series tail at most $(1/4)^{11}/[11!(1-1/2496)]<10^{-14}$. The Fraction script verifies the constants in all three bounds. These estimates cover analytic truncation. The paired GPU and CPU comparisons report numerical error at their listed inputs.

## Reproduction

Run in this directory with the data volume mounted. The original M6 used macOS 27.2 build 26B123 and Apple Swift 6.4.2. A timed replay may run a different number of batches. A future scheduled launch can use QUARTIC_START_EPOCH; QUARTIC_START_INDEX optionally selects an aligned counter offset.

```sh
QUARTIC_START_INDEX=2147483648 ./run.sh 90 128 main
PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python verify.py quartic-main.json cpu-main.json
./paired/run.sh
/Users/zzh/.venvs/main/bin/python cutoff_certificate.py
```
