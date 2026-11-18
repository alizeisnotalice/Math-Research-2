# Gaussian transfer at the critical square-root oscillation scale

The model is the density \(m^{-2}[1+(1/2)\cos(2\sqrt m)]\) on \(m\ge1\), and \(K(t,m)\) is the Gaussian density at \(t\) with mean \(m\) and variance \(m/d\), where \(d=2\). The computed quantity is \(I(t)=t^2\int_1^\infty K(t,m)m^{-2}[1+(1/2)\cos(2\sqrt m)]dm\).

Let \(R=\sqrt t\), \(y=(m-t)/\sqrt m\), \(h=y/(2R)\), \(v=\sqrt{1+h^2}\), and \(u=v+h=\sqrt m/R\). The map sends \([1,\infty)\) to \([1-t,\infty)\) and has Jacobian \(dm/dy=m/(Rv)\). The exact integral becomes
\[
 I(t)=\sqrt{d/(2\pi)}\int_{1-t}^\infty
 \frac{e^{-dy^2/2}}{u^3v}\left[1+\tfrac12\cos\left(2R+y+\frac{y^2}{2R(v+1)}\right)\right]dy .
\]
The phase is evaluated as \(\pm\cos[y+y^2/(2R(v+1))]\) at \(2R=2\pi k\) or \((2k+1)\pi\). For \(R\ge1\), the amplitude \((u^3v)^{-1}\) is at most one if \(y\ge0\), and at most \((1+|y|)^3\) if \(y<0\). Thus dominated convergence against the Gaussian gives \(I(t)\to1\pm\tfrac12e^{-1/(2d)}\), with limits 1.3894003915357024 and 0.6105996084642976.

| \(k\) | Phase \(2R\) | Measured \(I(t)\) | Measured minus limit |
|---:|---|---:|---:|
| 10 | \(2\pi k\) | 1.390097362068790 | \(6.96970533\times10^{-4}\) |
| 10 | \((2k+1)\pi\) | 0.610886450960795 | \(2.86842496\times10^{-4}\) |
| 50 | \(2\pi k\) | 1.389428269303558 | \(2.78777679\times10^{-5}\) |
| 50 | \((2k+1)\pi\) | 0.610612009901869 | \(1.24014376\times10^{-5}\) |
| 100 | \(2\pi k\) | 1.389407360969431 | \(6.96943373\times10^{-6}\) |
| 100 | \((2k+1)\pi\) | 0.610602739757744 | \(3.13129345\times10^{-6}\) |

The two quadratures use the exact \(y\) transform and the independent shift \(s=\sqrt m-R\). In the latter, \(y=s(2R+s)/(R+s)\), and the amplitude becomes \(2(R/(R+s))^4\). Tanh–sinh in \(y\) and Gauss–Legendre in \(s\) were compared at 90 decimal digits. For cutoffs \(|y|\le a\), the omitted amount is bounded by \((27/4)(a^2+1)e^{-a^2}\) when \(a\ge1\). This follows by summing the bounds on both wings and using \((1+y)^3+1\le9y^3\) for \(y\ge1\). The cutoff parameters are \(a=20\) and \(a=18\), with tail bounds below \(5.19\times10^{-171}\) and \(4.27\times10^{-138}\), respectively. Agreement between quadrature values is a numerical cross-check.

At these phases, the term of order \(1/R\) integrates to zero by oddness. With \(V=1/d=1/2\), the coefficient of \(1/t\) is
\[
 V\pm\tfrac12e^{-V/2}\left[(V-V^2)+\tfrac38(3V^2-V^3)-\tfrac1{32}(3V^2-6V^3+V^4)\right].
\]
The bracket is \(247/512\) at \(d=2\). It follows by expanding \((u^3v)^{-1}=1-3y/(2R)+y^2/R^2+O(R^{-3})\) and the phase by \(y+y^2/(4R)+O(R^{-3})\), then using the cosine and sine derivatives of the Gaussian Fourier transform. The program records \(t[I(t)-\text{limit}]\) and the result after subtracting the correction term. This is an additional numerical comparison.

Use `python transfer.py 10 1 sample.json` to reproduce one parameter, or `bash run.sh` to run all six simultaneously on the MacBook's existing Python environment. Independent local calculations used Python 3.13.18 and mpmath 1.3.0. The files under `local/` record all six executed parameters.

The MacBook run completed with Python 3.12.13 and mpmath 1.3.0, with endpoints stamped 09:05:25–09:05:26 UTC. The six-worker timing intersection was 09:05:26.014788–09:05:26.361999 UTC, and their cumulative CPU time was 2.107607 seconds. Every supervisor and worker exit code was zero, and a separate PID query found none of them. All 22 transferred hashes were checked. The MacBook results agreed with the independent local output to over 60 reported digits. Both quadratures agreed at the working precision on all six inputs. `summary.json` contains their full decimal output and the associated process receipts.
