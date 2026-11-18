# Repaired Metal scan of paired quartic near-cancellation

The model compares the quartic core with variances $cu^4(1+b(u/r)^2\pm h u/r)$ against its symmetric baseline $cu^4(1+b(u/r)^2)$. Both observations are evaluated at $y=a=kh/r$. Write $p_s$ for the centered Gaussian density of variance $s$. Denote the two core densities divided by $p_{M+\varepsilon}(a)$ by $R,R_0$. The full scalar is $(R-R_0)/h^2$.

Let $t=u/r$, $V=\varepsilon+cu^4(1+bt^2)$, $Z=u^2/V$, $\ell=cu^4t/V$, $x=h\ell$, and $W=e^{-(u^2+Z)/2}/\sqrt V$. Pairing the two halves gives
$$
R-R_0=C_h\int_0^r W\left[e^{L_e}\cosh(kht-L_o)-\cosh(kht)\right]du
$$
with $C_h=2\sqrt{(M+\varepsilon)/(2\pi)}e^{-a^2(1-1/(M+\varepsilon))/2}$ and
$$
L_e=-\tfrac14\log(1-x^2)-\frac{Zx^2}{2(1-x^2)},\qquad
L_o=\tfrac12\left[\frac{Zx}{1-x^2}-\operatorname{atanh}x\right].
$$
The leading scalar is $C_h(A-kB)$, where $A=\int_0^r W(Z^2-6Z+3)\ell^2/8\,du$ and $B=\int_0^r Wt(Z-1)\ell/2\,du$. The tuned family samples offsets around $k=A/B$.

The binary32 formula factors the leading cancellation before dividing by $h^2$. The repair uses the explicit odd Taylor polynomial for sinh at $|x|\le1/4$. This preserves small arguments in the sinh-over-$h$ factors. The unmodified study and its initial 32-lane diagnostic remain in the prior cc43 project directory. The new main scan reused the deterministic generator and extended its prefix with the corrected kernel. It used 96 Gauss–Legendre nodes per panel and a bound of 192 on $Z$ for negligible kernel evaluations. CPU quadrature has independent panels and uses a higher $Z$ limit.

## Actual main scan

The scan sampled logarithmic $r\in[2^{-6},1]$, $M-1\in[2^{-10},8]$, $\varepsilon\in[10^{-16},10^{-2}]$, and $h\in[2^{-14},2^{-2.5}]$, before applying the cap $|a|\le3$. It used $b\in[-3/8,3/4]$ and a logarithmic quartic edge scale normalized by $r^2$, with a variance cap. The tuned cases required $|B|>10^{-4}$ and $|A/B|<12$, and displaced $k$ by a signed log-sampled offset between $2^{-12}$ and $1/4$. Untuned cases sample $k$ uniformly in $[-12,12]$.

Apple M6, macOS 27.2 (26B123), Apple Swift 6.4.2. Main PID 11253 ran from 02:43:35.314 to 02:45:05.329 UTC on 2026-11-18 and exited with code zero. The 90.014699 s run evaluated 1,950,875,648 profiles. Metal command-buffer intervals totaled 85.369066 s, with command submission and waits totaling 86.384774 s. The 168 during-run systemwide GPU utilization samples had minimum 93%, average 99.827%, and maximum 100%, with 99% before and after. Utilization includes other work. There were zero nonfinite values.

The ledger's **separated signs** use the GPU rule that both absolute scalars are at least $2\times10^{-6}$. This is a numerical reporting threshold. Sampled CPU agreement establishes checks at the verified points. The strata overlap.

| Stratum | Profiles | Separated sign differences | Either scalar below threshold |
|---|---:|---:|---:|
| All | 1,950,875,648 | 16,964,119 | 46,277,627 |
| Tuned | 925,911,786 | 16,567,925 | 32,197,518 |
| Tuned, offset below 0.01 | 495,949,413 | 14,639,771 | 31,099,151 |
| Tuned, $h<0.01$ | 592,777,658 | 722,375 | 17,995,774 |
| Tuned, $h<0.001$, offset below 0.01 | 173,995,822 | 520 | 9,394,235 |
| Untuned | 1,024,963,862 | 396,194 | 14,080,109 |
| Tuned, negative $b$ | 298,950,606 | 7,825,086 | 9,980,826 |
| Tuned, nonnegative $b$ | 626,961,180 | 8,742,839 | 22,216,692 |

The inherited fields ratio and baseline denote the full-minus-leading absolute gap and the full scalar, respectively. The inherited above field counts all strict-sign comparisons, including values below the new reporting threshold. The robust and nearZero fields use the new threshold.

## Independent verifications

- SciPy double quadrature verified the representative with largest absolute gap in each of the eight strata. The maximum absolute discrepancy in the full scalar was $3.50\times10^{-7}$. For the lowest-$h$ stratum's representative it was $2.73\times10^{-7}$. The 65-digit verifier also retains the two exactly tuned CPU calibration profiles.
- A separate dispatch replayed 8,192 fixed counters and the previously unstable index. It selected 46 unique profiles with 90 stratum memberships, covering every stratum. The maximum absolute discrepancy across full and leading scalars was $1.84\times10^{-7}$, with no separated sign mismatches. At the previously unstable index, corrected GPU and CPU scalars were $-9.14985\times10^{-6}$ and $-9.10617\times10^{-6}$. The fitted $k$ uses order-96 quadrature, so it can differ from the previous order-68 diagnostic.
- The rare low-$h$ separated-flip class required a targeted early-filtered replay of 67,108,864 counters from the main prefix. It retained 17 cases in 2.185329 s of Metal commands. Independent double quadrature reproduced both signs in all 17. Fifteen retained the $2\times10^{-6}$ threshold in double precision and two were slightly below it, with no sign changes. Maximum absolute discrepancy was $2.49\times10^{-7}$.

Audit and rare-replay PIDs 14621 and 15975 exited successfully and were verified absent along with the main process. The replays are checks on the main prefix.

## Reproduce

Run from this directory on Apple Silicon. A repeat of the time-bounded main run can evaluate a different prefix length. The sparse audit's full input is generated on replay.

```sh
./run.sh 90 96 main
/Users/zzh/.venvs/main/bin/python verify.py quartic-main.json cpu-main.json
./audit/run.sh
./rare/run.sh
```
