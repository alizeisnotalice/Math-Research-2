# Strict positivity at the lowered rational right anchor

The original eleven coefficients, variance range \(0\le f\le5/2\), and \(1\)-Lipschitz requirement on \([-3/4,3/4]\) are fixed. Let \(q(x,u)\) be the original signed combination of one mean, four time tests and six heat-kernel tests, and let \(b_*\) be the exact target for that combination. As before, for a rational anchor \(A\) at \(1/2\), bound the relaxed margin
\[
 M(A)=\int \inf_{\ell_A(x)\le u\le 5/2} q(x,u)\,\gamma(dx)-b_*,
 \qquad
 \ell_A(x)=
 \begin{cases}
 \max(0,A-|x-1/2|),& |x|\le3/4,\\
 0,& |x|>3/4.
 \end{cases}
\]
A strictly positive margin rules out exact graphs within this anchored class and implies a strictly positive uniform gap among the original tests. Gaussian tail bounds allow arbitrary admissible values on both tails.

The first targeted refinement uses spatial and height denominators 4096 and the original exponential argument denominator 4096. Its enclosure of \(M(171/100)\) crosses zero. The successful second calculation doubles the spatial and height denominators to 8192 and increases the exponential denominator to 32768. Directed rational results, with endpoints here displayed outward, are

| Right anchor | Directed enclosure of \(M(A)\) |
|---:|---:|
| \(3417/2000=1.7085\) | \([0.0004295266507,\ 0.0014985875900]\) |
| \(1709/1000=1.709\) | \([0.0007049127540,\ 0.0017743483990]\) |
| \(171/100=1.71\) | \([0.0012561769191,\ 0.0023259157728]\) |

The signed tail allowances in the fine run are below \(1.349\times10^{-12}\) for the lower bound and \(3.684\times10^{-8}\) for the upper bound. If the left anchor remains at least \(3/4\), the new threshold \(3417/2000\) corresponds to an amplitude requirement of \(1917/2000\). The threshold is the lowest verified by this fine sweep, while smaller anchors remain unsettled by the directed calculation. The floating-point check uses a dense height grid, local ternary refinement, and separate spatial integration across the modulus boundaries. It is a cross-check of the sign and remains separate from the integer certificate.

## Box pruning

On a fixed spatial cell, round the weakest admissible lower envelope down to the \(u\)-grid for the lower bound. Round the strongest admissible envelope up to construct a feasible constant height for the upper bound. The incumbent upper bound is an interval upper evaluation at such a feasible height over the whole spatial cell. A height box is discarded only if its directed lower bound is at least the current incumbent. Since subsequent incumbents are nonincreasing, every discarded box remains bounded below by the final incumbent. Taking the smaller of the final incumbent and all undiscarded terminal lower bounds therefore gives a valid lower bound. For the upper bound the final incumbent still comes from a feasible height. All integrand bounds use the original rational exponent and square-root machinery. Only the grid resolution and enumeration strategy changed.

## Runs

Rerun the initial calculation with **python refine.py** and the finer calculation with **python refine.py --nx 8192 --nf 8192 --exp-grid-log2 15 --prefix fine 3417/2000 1709/1000 171/100**. Run **python calibrate.py** for the independent floating calculation and **python audit.py** to reproduce exact arithmetic from the worker numerators and check process exits. The 12-worker fine interval was 07:23:11.193526–07:25:15.872163 UTC. The initial 12-worker interval was 07:20:40.336060–07:21:30.959892 UTC. Both started after the reported synchronized GPU runs ended.
