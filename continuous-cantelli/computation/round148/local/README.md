# Local Metal two-well phase test

The two positive flanks from project round 144 have gap parameters d=1,2, baseline weights A=1/2 each, and cosine weights B1=1/4, B2=-exp(-1/4)/4. The scan samples t=π²(n+u)² for integer n from 16 to below 2³⁰ and exact quarter turns or pseudorandom dyadic u. Each GPU evaluation does both rescaled Gaussian quadratures. It retains the largest raw sine error, raw scaled-residual error and centered second-order error in each exponent and phase stratum.

The normalized two-flank sum has expansion

    Z(t) = 1 + (13 exp(-1/2)/64) sin(2 sqrt(t))/sqrt(t)
             + [3/4 + (169 exp(-1/2)/2048) cos(2 sqrt(t))]/t + O(t^(-3/2)).

At the exact zero-sine phases the second coefficient is positive on both sequences. The asymptotic values are about 0.800050626 and 0.699949374. The opposite quarter turns give scaled limits ±0.123201540. The derivation in expansion.md also controls the exponentially small change from the constant-gap transfer to the actual round 144 flank. This calculation concerns the prescribed flank outputs.

The Metal stable path takes the phase turn separately and integrates the dephased second remainder in a form that avoids subtracting order-one quadrature sums. The CPU comparator uses 65-digit mpmath working precision, independently integrates the uncompensated oscillating kernel on [-25,25] and uses the exact untruncated constant-gap formula for the mass. At sampled n>=16 the omitted parts of the constant-gap oscillating integral have exponentially small tails. Across 60 saved production extrema the absolute discrepancy of the scaled second remainder is at most 3.695e-7.

To isolate quadrature, the verifier computes the imperfect leading cosine cancellation using high-precision nodes and weights before float32 arithmetic. It is about 6.998e-16 at order 32 and -8.900e-30 at orders 64 and 72 on the shader's interval [-11,11]. The GPU evaluation of this cancellation at order 64 has a residual up to 1.490e-8. As t increases, direct subtraction can magnify even a small constant floor. The centered remainder removes that floor analytically.

The supplemental phase run records the float32 angle passed to the Metal sine. Using the high-precision sine of the recorded angle, input-angle representation accounts for differences as large as 1.98208 from the intended sine among the saved samples. Evaluating Metal sine and high-precision sine on the same recorded angle also differs by up to 0.998926 there. The raw Metal sine returns zero in the saved largest-band cardinal cases, including intended sine +1 and -1. The actual angle is available in the supplemental file; the main run audits its reconstruction from epsilon.

## Measured execution

The Apple M6 production ran 2026-11-18 09:11:52.684–09:12:47.684 UTC. It completed 2,506,555,392 two-flank evaluations in 55.000384 s, with 49.306723 s of attributable Metal command time. The 105 systemwide IOAccelerator samples had median 96% and maximum 97%. PID 90777 exited with status zero and was absent.

The 64-node supplemental run started at 09:14:34.099 UTC and ended at 09:14:48.099. It completed 733,347,840 evaluations in 14.000114 s, with 12.226666 s Metal command time. The device utilization median was 95%, maximum 97%, over 28 systemwide samples. PID 92307 exited with status zero and was absent.

## Reproduction

Run from this directory on Apple Silicon with the command line developer tools and mpmath. The launcher respects the external data volume and uses its shared Swift module cache.

    ./run.sh 0 32 smoke
    ./run.sh 55 72 main
    ./phase/run.sh 14 64 phase
    TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python verify.py wells-main.json audit-main.json ioreg-main.txt
    TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python verify.py phase/wells-phase.json phase/audit-phase.json phase/ioreg-phase.txt

The runtime budgets can change the number of evaluated rows and retained extrema in a rerun.
