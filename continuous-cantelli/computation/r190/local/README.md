# Local Metal alignment check for complex Cantor heat averaging

The experiment uses the middle-thirds Cantor probability on [2,3]. It compares the positive Gaussian Fourier energy with the same Gaussian-weighted spatial kernel averaged over genuine middle gaps of depth-j cylinders, and also with the midpoint substitute. Each comparison holds the heat variance b and complex shift D fixed. All outputs are normalized by g(2)^2 p_b(0) exp(D²/(2b)). The proof in coupling.md gives explicit first-order and symmetric second-order spatial error bounds.

The main GPU run sweeps 240 cases at resonances D/b=2π3^m, for m from 0 to 9, gap depths j through 12, six values of sqrt(b)3^m, and Cantor depth 6 through 14. It integrates each actual gap with 16 nodes. Resonant phases use integer triadic differences modulo the exact period. This preserves phase at fine depths even when physical coordinates have few remaining Float32 digits.

In the sampled sweep, the real gap value is negative in 46 cases, of which 35 have j=m and 11 have j=m−1. Midpoint and true-gap real signs differ in 25 cases. At m=j=7 and sqrt(b)3^m=3.5, the normalized energy is +0.00074423737, the gap real part is -0.00165034284 and the midpoint real part is -0.00199753849. The maximum relative true-gap errors at j−m=1,2,3,4 are 0.42513, 0.03449, 0.003726 and 0.0004082. The conservative Lipschitz screen succeeds in six cases and the second-order screen in 57, after subtracting a bound on the truncation error from the measured energy. The maximum observed error as a proportion of either conservative bound is below 4.31%.

A 10-case companion fixes the scalar variance for each case at the midpoint radial age, b=(β−1)L^-β+(L+log 6)^-β, with L=j log 3 and β=2,3,4. It keeps feasible exact resonances with 1−log 2/log 3<b(D/b)^2/(2L)<1. This tests constant-variance sections at the theorem's scales. All 10 real gap values are positive and pass the symmetric second-order screen. The five finest of those evaluated cases also pass the Lipschitz screen.

Independent SciPy and NumPy direct spatial integrals check selected resonances with 24-node gap quadrature and one extra Cantor depth. The sampled maximum absolute Metal discrepancy is 1.461e-8 for the main run and 2.154e-7 for the companion. At the moderate resonances, independently integrating the nonnegative Fourier spectrum and taking the direct double spatial sum agree within 4.17e-17. The finest differences in the companion are close to its Float32 accumulation error. The separate interval certificate directly bounds the infinite-Cantor coarse-gap reversal and the fine-gap positive lower bound.

## Execution

Apple M6 main PID 38453 evaluated 894,336 outer spatial rows at 2026-11-18 11:35:33.821–11:35:33.841 UTC. Native GPU time was 0.017861 s within 0.020258 s elapsed. The theory-shaped companion PID 40097 evaluated 8,192 rows at 11:40:30.346–11:40:30.350 UTC, with 0.002879 s native time. Both production processes and the separate smoke process exited zero and are absent.

## Reproduction

Run on Apple Silicon with the existing Swift tools, NumPy, SciPy and mpmath. The launch script builds with the shared Swift cache on the external data volume.

    ./run.sh smoke
    ./run.sh main
    ./run.sh theory
    TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python audit.py report-main.json audit-main.json
    TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python audit.py report-theory.json audit-theory.json
    TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python certify.py
