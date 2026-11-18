# Certified fractional-kernel crossings

For \(0<\kappa<1\), define
\[
 G(z)=\int_0^\infty t^{-\kappa}e^{zt-t^2/2}\,dt,
 \qquad D(z)=zG(z)-G'(z).
\]
Integration by parts in \(\partial_t(t^{1-\kappa}e^{zt-t^2/2})\) gives \(G''=zG'+(1-\kappa)G\), hence \(D'=\kappa G>0\) on the real line. Also \(D(0)=-\mu_1<0\) with \(\mu_n=\int_0^\infty t^{n-\kappa}e^{-t^2/2}dt\). For \(z\ge0\), \(D'\ge\kappa\mu_0\), so the crossing is unique and strictly positive.

| \(\kappa\) | Certified crossing, displayed rounded |
|---:|---:|
| \(1/4\) | 1.364492202705106986072317622127523 |
| \(1/2\) | 0.764950867389759566027391419412589 |
| \(3/4\) | 0.341463306056180180661343392823726 |
| \(99/100\) | 0.012571817875339041413802094038992 |

Each certified open interval has width \(3\times10^{-90}\). See `summary.json` for the exact decimal endpoints and outward sign enclosures.

## Rational certificate

The moment recurrence is \(\mu_{n+2}=(n+1-\kappa)\mu_n\). Thus
\[
 D(z)/\mu_0=-r+\kappa\sum_{n\ge0}(\mu_n/\mu_0) z^{n+1}/(n+1)!,
 \qquad r=\mu_1/\mu_0 .
\]
For positive \(z\), the series terms are nonnegative. Each parity chain has ratio \((n+1-\kappa)z^2/[(n+2)(n+3)]\) from index \(n\) to \(n+2\). After retaining indices \(0,\ldots,207\), the later ratios are bounded by \(z^2/211\). The first omitted term of each chain divided by \(1-z^2/211\) is therefore a bound for that chain's remaining sum. The odd chain is evaluated with its base factor \(r\) separated. The largest displayed pre-multiplier tail bound was less than \(1.565\times10^{-172}\).

For \(a_0=(1-\kappa)/2\), \(a_1=(2-\kappa)/2\), the identity \(r=\sqrt{2}\,\Gamma(a_1)/\Gamma(a_0)\) follows by changing variable to \(t^2/2\). For each shape \(a\) between zero and one put \(M=512\). Integrating the odd and even Taylor sums for \(e^{-u}\) on \([0,M]\) encloses the part of \(\Gamma(a)/M^a\) on that interval by the explicit sums \(\sum_{n=0}^{N}(-M)^n/[n!(a+n)]\), with \(N=2047,2048\). The remaining positive tail is bounded by \(e^{-M}/M<2^{-M}/M\). Since \(a_1-a_0=1/2\) and \(\sqrt{2M}=32\), the ratio of the resulting finite rational enclosures multiplied by 32 encloses \(r\) without any transcendental power evaluation. The widest ratio interval was less than \(2.082\times10^{-155}\).

These bounds are propagated with exact Fraction arithmetic to both endpoints of each root interval. The derivative sign converts the certified negative lower-endpoint sign and positive upper-endpoint sign into the asserted crossing enclosure. An approximate gamma value and numerical root merely select test endpoints. An independent comparison integrates \((z-t)t^{-\kappa}e^{zt-t^2/2}\) using \(y=t^{1-\kappa}\) near the singularity and evaluates it at the approximate crossing. For \(t\ge40\), \(z\le2\) and the tested \(\kappa\), the absolute integrand is at most \(2t e^{-0.45t^2}\). Its contribution after normalization by \(\mu_0\) is less than \((40/9)e^{-720}\). At 122 digits the direct numerical comparison differed from the moment evaluation by at most \(3.16\times10^{-123}\) in normalized \(D\).

## MacBook run

`bash run.sh` runs the four `certify.py` jobs using Python 3.12.13 and mpmath 1.3.0. The corrected production run had all-worker intersection 11:27:42.374991–11:27:44.114741 UTC and used 8.433488 cumulative worker CPU seconds. All final processes exited zero, a separate process query found no owned PIDs, and 18 transferred hashes matched. In the initial check, two independent gamma comparisons attempted to resolve a \(10^{-155}\)-wide enclosure at 122-digit working precision. This comparison was rerun at 182 digits. The exact rational certificate is independent of those comparisons. `run-certified/`, `summary.json`, `remote-sha256.txt`, and `final-ps.txt` hold the reproducibility receipts.
