# A joint time–space infeasibility certificate on a fixed-amplitude compact class

Consider continuous \(f:\mathbb R\to[0,H]\), \(H=5/2\), with the common modulus \(h\mapsto h\) on \([-3/4,3/4]\). Fix the rational points \(x_1=-1/2\), \(x_2=1/2\) and impose \(f(x_1)\in[3/4,13/16]\), \(f(x_2)-f(x_1)\ge1\). To obtain a specific compact class in the topology of local uniform convergence, also require the common modulus \((n+2)h\) on \([-n,n]\) for every positive integer \(n\). Arzelà–Ascoli and closedness give compactness. Only the modulus on \([-3/4,3/4]\) is used by the certificate. In particular, the calculation accounts for all values in \([0,H]\) beyond that interval and does not prescribe either tail.

The normalization constraint \(\mathbb Ef=1\) is attainable in this compact class. The tent \(f_{\rm tent}(x)=(7/4-|x-1/2|)_+\) and ramp \(f_{\rm ramp}(x)=3/4+\min(1,(x+1/2)_+)\) both attain the stated anchor values and meet the moduli. The tent has a certified mean \(<0.92\), while the ramp has exact mean \(5/4\) by reflection. An appropriate convex combination therefore has mean one.

## Tests and strict conclusion

For \(s>0\) use \(k_v(u)=v^{-1/2}e^{-u^2/(2v)}\); this is \(\sqrt{2\pi}\) times the Gaussian density of variance \(v\). If
\[
 F(t)=\mathbb E_\gamma e^{tX-t^2f(X)/2}=1\quad\hbox{for all real }t,
\]
the entire extension gives \(F(i\xi)=1\), and differentiation at zero yields \(\mathbb Ef=1\). The characteristic function of \(X+\sqrt{H+s-f(X)}G\), with independent standard Gaussian \(G\), is \(e^{-(H+s)\xi^2/2}F(i\xi)\). Consequently
\(\mathbb E_\gamma k_{H+s-f(X)}(y-X)=k_{H+s}(y)\). This establishes all the finite constraints below without assuming positivity of any intermediate signed heat remainder.

Use a moment test \(f(x)\), four real-time tests \(e^{t x-t^2f(x)/2}\) at \(t=-2,-1,1,2\), and six smoothed space-time tests \(k_{H+s-f(x)}(y-x)\). Their integer weights, divided by \(10^6\), are given in order.

| Test | Parameters | Weight numerator |
|---|---:|---:|
| Mean | | 1000000 |
| Real | \(-2\) | 34939 |
| Real | \(-1\) | 361498 |
| Real | \(1\) | \(-15824\) |
| Real | \(2\) | 128578 |
| Heat | \((1/4,-3/2)\) | \(-588210\) |
| Heat | \((1/4,-1/2)\) | \(-1000000\) |
| Heat | \((1/4,1/2)\) | 1000000 |
| Heat | \((1/4,3/2)\) | \(-964943\) |
| Heat | \((3/4,-1)\) | \(-182022\) |
| Heat | \((3/4,1)\) | 1000000 |

Write \(q(x,f)\) for their signed weighted sum and \(b\) for the weighted sum of exact targets. The refined rational certificate uses only \(f(1/2)\ge7/4\), hence applies to the original fixed-amplitude class and to this larger one-anchor class. It proves
\[
 \mathbb E_\gamma q(x,f)-b>0.01289603,\qquad
 \max_j|\mathbb E_\gamma q_j-b_j|>1/500 .
\]
The second claim divides the strict lower bound by the sum of absolute rational weights. The strict tail error in this test is less than \(1.35\times10^{-12}\).

## Directed computation

On each spatial cell \([x_0,x_1]\) inside the Lipschitz window, the right-anchor requirement yields the lower bound \(f\ge\max(0,7/4-\max(|x_0-1/2|,|x_1-1/2|))\). Elsewhere use \(0\le f\le H\). The first certificate uses both anchor envelopes. The stronger one uses this single lower envelope. Round the lower and upper bounds outward on a dyadic \(f\)-grid. For each covered product box, lower bound \(e^{-x^2/2}q\) by signed endpoint intervals. Positive coefficients use termwise lower bounds and negative coefficients use upper bounds. For the time terms the exponent uses
\[
 -x^2/2+tx-t^2 f/2=t^2/2-(x-t)^2/2-t^2 f/2.
\]
For heat terms, bound the positive prefactor \(1/\sqrt v\) and the negative exponent with their variance endpoints. The source factor \(e^{-x^2/2}\) is also enclosed. Take the minimum over all \(f\)-boxes for each spatial cell and sum with cell lengths. This yields a lower bound even when \(f\) varies arbitrarily within a cell.

Every numerical operation in this certification is directed integer or rational arithmetic. For the exponential, monotonicity rounds the nonnegative argument to endpoints on the grid \(k/4096\). The dyadic interval for \(\exp(-1/4096)\) comes from adjacent seventh and eighth alternating Taylor sums. Outward-rounded multiplication gives the table for all subsequent arguments. For arguments beyond 40, take zero below and the upper enclosure at 40 above. Reciprocal exponentials and inverse square roots use directed division and integer square roots. The Gaussian normalizer is enclosed by the alternating arctangent series in Machin's identity.

To bound the omitted region \(|x|>8\), each nonnegative term with positive weight may be discarded. A negative time weight pays at most its absolute weight times \(2e^{t^2/2-(8-|t|)^2/2}\) by completing the square and Chernoff. A negative heat weight pays at most its absolute weight times \(2s^{-1/2}e^{-32}\). The moment weight here is positive. The sum of those constants accounts for all possible tails. Enclose the target \(b\) from above and subtract it and the tail bound from the integrated lower bound.

The 12-worker discretized LP used to find the weights is exploratory. On the selected case, the sampled relaxation gave zero positive margin to each test group separately and a positive joint margin. That observation makes no claim about separate feasibility of actual continuous graphs. The proof uses only the subsequent rational checks. Both real-time and smoothed space-time terms occur in the verified certificate.

## Reproduction

Run **python search.py --duration 33** for the exploratory LP, **python certificate.py** for the two-anchor verification and **python certificate.py --nx 512 --nf 256 --right-anchor-only --output-prefix certificate-right** for the wider class. Then run **python audit.py** to recheck the exact rational residual arithmetic, the normalized tent-mean comparison and the process exits. The first-run telemetry is also available under the current filename **certificate-telemetry.jsonl**. Both certificates concern the given modulus and rational anchor regime. They leave the general continuous Cantelli question open.
