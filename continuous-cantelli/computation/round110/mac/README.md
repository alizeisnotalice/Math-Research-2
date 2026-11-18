# Gaussian real transform for four-well profiles

Write \(F_f(\lambda)=\mathbb E\exp(\lambda X-\lambda^2 f(X)/2)\) for \(X\sim N(0,1)\). A normalized exact graph has \(F_f=1\). The family uses \(u>0\), \(s=e^{-u x^2/2}\), and
\[
P(s)=s\prod_{j=1}^4(s-r_j),\qquad f(x)=1+\epsilon P(e^{-u x^2/2}).
\]
Four distinct roots in \((0,1)\), \(\epsilon>0\), and \(f>0\) give two negative lobes on each half-line and \(f(0)>1\). Between successive roots, \(P'/P=1/s+\sum_j(s-r_j)^{-1}\) decreases from \(+\infty\) to \(-\infty\). It has exactly one zero in each such interval, and each negative lobe has a nondegenerate minimum.

## Exact moments and finite Jensen test

Let \(g_j=(1+ju)^{-1/2}\). For three input roots \(a,b,c\), set \(A=a+b+c\), \(B=ab+ac+bc\), \(C=abc\). Solving \(\mathbb E P(s)=0\) gives
\[
r_4=\frac{g_5-Ag_4+Bg_3-Cg_2}{g_4-Ag_3+Bg_2-Cg_1}.
\]
With \(P(s)=\sum_j p_js^j\), write \(J=\sum_j p_j(-ju)g_j^3\) and \(V=\sum_{i,j}p_ip_jg_{i+j}\). The fourth coefficient of \(F_f-1\) is \(-\epsilon J/4+\epsilon^2V/8\). Therefore \(\epsilon=2J/V\) matches it exactly whenever \(\epsilon>0\). The odd coefficients vanish by symmetry. SymPy also verifies the sixth coefficient
\[
-\epsilon\mathbb E[P(s)H_4(X)]/48
+\epsilon^2\mathbb E[P(s)^2H_2(X)]/16
-\epsilon^3\mathbb E[P(s)^3]/48.
\]
Here \(H_j\) are probabilists' Hermite polynomials.

Set \(t=\lambda^2/2\), \(\kappa=u/(1+u)\), and
\[
J_\lambda=\sum_{j=1}^5p_jg_j\,e^{-t(ju/(1+ju)-\kappa)}.
\]
This is \(e^{t\kappa}\mathbb E_{N(\lambda,1)}P(s)\). Jensen gives
\(\log F_f(\lambda)\ge-t\epsilon e^{-t\kappa}J_\lambda\). The exact rational-decimal input
\((a,b,c,u)=(.020179268,.52300256,.97971076,.14436331)\)
has two positive radial minima, in \((.79997,.80000)\) and \((.94263,.94265)\). At \(\lambda=6\), a 55-decimal-digit outward interval encloses \(J_\lambda\) near \(-.01904465000028334\) and proves \(\log F_f(6)>1.05511970393650\). The verifier encloses the unique minimum in each negative lobe with interval tests on the monotone logarithmic derivative. The exact Gaussian-moment formulas impose normalization and fourth-moment matching.

## Whole-family tail obstruction and late zero

With \(\gamma_\lambda=N(\lambda,1)\), split \(F_f-1=S_\lambda-D_\lambda\), where \(S\) integrates \(e^{t(1-f)}-1\) on \(f<1\) and \(D\) integrates \(1-e^{-t(f-1)}\) on \(f>1\), both against \(\gamma_\lambda\). Let \(\alpha=\min f\), and let \(x_*>0\) be the rightmost global minimum. Then
\[
S_\lambda\sim C_S\lambda^{-1} e^{-\alpha\lambda^2/2+x_*\lambda},
\quad C_S=p_1(x_*)\sqrt{4\pi/f''(x_*)}\,e^{1/f''(x_*)},
\]
\[
D_\lambda\sim C_D\lambda^2 e^{-\kappa\lambda^2/2},
\quad C_D=\frac{\epsilon\prod_j r_j}{2\sqrt{1+u}}.
\]
For \(S_\lambda\), set \(x=x_*+y/\lambda\) in a small fixed neighborhood of each minimum. Quadratic localization controls the rescaled integral by \(C e^{|y|-c y^2}\). All other minima are exponentially smaller in \(\lambda^2\) if their heights differ, or exponentially smaller in \(\lambda\) if their heights agree and they lie farther left. Away from the finitely many minima there is a positive height gap on the bounded support of \(f<1\). The subtracted constant on \(f<1\) is exponentially smaller since \(\alpha<1\).

For \(D_\lambda\), \(f-1\sim \epsilon(\prod_jr_j)e^{-u x^2/2}\) as \(x\to+\infty\). Replace \(1-e^{-v}\) by \(v\) using \(0\le v-(1-e^{-v})\le v^2/2\), and complete the square in the shifted Gaussian. Terms with \(s^2\), including the error from linearization, have strictly larger Gaussian decay exponent. The bounded region omitted when replacing the positive part by the tail gives \(O(\lambda^k e^{-\lambda^2/2+O(\lambda)})\). These estimates are uniform on a compact parameter neighborhood with separated roots, positive curvature and a positive gap to the second radial minimum.

It follows that
\[
\log(S_\lambda/D_\lambda)=-\tfrac12(\alpha-\kappa)\lambda^2+x_*\lambda
-3\log\lambda+\log(C_S/C_D)+o(1).
\]
For \(\alpha>\kappa\), the transform gap is eventually negative. For \(\alpha\le\kappa\) it is eventually positive, including equality by \(x_*>0\). This excludes each profile in the stated family. The argument extends the single-minimum computation of R007 to these four-well profiles, including tied heights.

In a compact separated branch, let \(q=\alpha-\kappa\downarrow0^+\), \(M=2x_*/q\). The ratios at \(\lambda=M/2\) and \(2M\) are respectively above and below one for small \(q\), by the uniform expansion. The transform is real analytic and eventually strictly below one, so there is a largest positive zero of \(F_f-1\). That zero diverges and obeys \(\lambda_{\rm last}/M\to1\) by the same expansion. Substitution into the zero equation first gives \(\lambda_{\rm last}-M=O(\log M)\), then
\[
\lambda_{\rm last}=M-\frac3{x_*}\log M+\frac1{x_*}\log(C_S/C_D)+o(1).
\]
The interval certificate in the CPU output verifies such a branch with fixed \((a,b,c)=(.19346052,.6506234,.71849674)\). At \(u=1.8654002890\) it encloses \(q\) near \(+7.07532\times10^{-11}\); at \(u=1.8654002892\) it encloses \(q\) near \(-8.39616\times10^{-11}\). Both radial minima remain positive and separated throughout that interval. Continuity supplies \(q\to0^+\) on this family. This is a delay in the far sign crossover of these profiles; moderate parameters can still detect their failure.

## Metal scan

The scan samples deterministic input triples by three distributions and stratifies \(u\) over \([1/8,16]\) or \([1/4,8]\). The GPU screens the rounded roots, normalization, amplitude and apparent minimum. It evaluates the full shifted-Gaussian kernel over \([-14,14]\) about \(\lambda\), with 224 midpoint panels and 12 values of \(\lambda\) from 0.25 to 12. The numerical diagnostic is the scaled residual \((F_f-1)e^{t\kappa}/(t\epsilon)\), splitting out the exact-form linear term and integrating the nonnegative exponential remainder.

The run screened 707,788,800 proposals and evaluated 109,519,357 accepted profiles, totaling 1,314,232,284 transform evaluations. It recorded 1,772,233 negative shifted-Jensen diagnostics. Metal command buffers took 73.048479416 s across 88.732079167 s wall. IOAccelerator system-wide utilization samples were 155 readings at 100% and seven at 99%. The executable ran on Apple M5 Pro from 05:58:44.330 to 06:00:13.062 UTC on 2026-11-18 and exited zero. The smoke test exited zero. Both PIDs were absent in the post-run check.

The native screening uses Float32 and its moment sums are cancellation-sensitive. In particular the two smallest reported \(|q|\) replay as approximately \(3.1840\times10^{-4}\) and \(6.9972\times10^{-5}\), both positive. The best reported normalized grid score \(1.3219\times10^{-5}\) replays as \(3.7555\times10^{-4}\). The reported six grid sign changes for a selected profile reduce to two in the CPU replay. Six selected profiles were recalibrated from the saved input triple and \(u\) in 65-digit arithmetic, followed by independent adaptive integration on the wider shifted window \([-18,18]\). SciPy and the selected 65-digit quadratures differed by at most \(3.854\times10^{-13}\). The MacBook Python 3.12.13 and local Python 3.13.18 produced identical replay outputs. Grid counts describe numerical screening. The explicit interval witnesses and the analytic tail argument carry the strict conclusions.

### Reproduce

On the authorized MacBook, in the source directory:

```sh
cache="$HOME/CodexRemoteMath/cache"
mkdir -p "$cache/swift-modules" "$cache/cc101-temp"
TMPDIR="$cache/cc101-temp" swiftc -O -module-cache-path "$cache/swift-modules" scan.swift -o multiwell_scan
sh run.sh main 224 73 524288 33718319
"$HOME/.venvs/math/bin/python" check_cpu.py main.json main.cpu.remote.json
```

The command-time cutoff makes batch count depend on device throughput. Check `main.stdout.log` for exact batch and sample counts. Use the saved `main.json` to replay the listed numerical hits with `check_cpu.py`.
