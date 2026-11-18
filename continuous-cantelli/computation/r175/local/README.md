# R170 disconnected-peak budget calibration

The Metal scan tests centered or displaced pairs of disconnected trapezoidal strict superlevels. It varies interval width and center, selected level, the relative widths and separating gap, asymmetry, and plateau heights and shoulder slopes. Six strata cover broad parameters, a cap at graph height 2, symmetric peaks, a neighborhood of a rational design, narrow tall peaks, and a dyadic grid. For each candidate the program integrates each affine shoulder and flat plateau separately. It searches for strict positivity of the mass, determinant of the degree-two moment matrix, and concave quadratic slack, together with a negative first Dirichlet sine slack at the same level.

At the rational certificate, s=5/4 and I=[-1/2,1/2]. The left peak occupies normalized [0,3/32] with excess height 3/4 and plateau occupying fractions 1/8 through 7/8 of its support. The right peak occupies [27/32,1] with excess height 1/2 and plateau fractions 1/8 through 3/4. The interval certificate gives sine slack below -0.07384 while all nonnegative quadratic budgets pass strictly. It also certifies a scaled family at each 0<L<=1/4, with the same fractions, physical heights multiplied by L², and normalized sine slack below -0.09594. A continuous nonnegative mean-one completion exists with sup f<=2. See proof.md for the structural comparison and precise signs.

An adaptive CPU integration independently confirms the direction of every sign for all six saved Metal extrema. The largest absolute GPU differences are 6.53e-8 in mass, 2.52e-8 in the concave quadratic slack, 1.52e-8 in the determinant and 7.99e-6 in the sine slack. The large Float32 scan provides candidate and error information. The strict certified family uses interval Riemann boxes with rational endpoints, rather than Float32 test counts.

The Apple M6 production ran from 2026-11-18 11:02:09.687 to 11:02:34.688 UTC and executed 4,252,893,184 profile evaluations at order 16. Metal command time was 20.208375 s out of 25.000744 s elapsed. Systemwide IOAccelerator samples during the run had median 96% and maximum 98% in 47 samples. PID 26560 exited with status zero and was absent in the recorded exit check. The order-8 smoke run also exited normally.

## Reproduction

Run on Apple Silicon with the existing Swift tools, NumPy, SciPy and mpmath. The build uses the shared Swift module cache on the mounted data volume.

    ./run.sh 0 8 smoke
    ./run.sh 25 16 main
    TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python audit.py peaks-main.json audit-main.json ioreg-main.txt
    TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python certificate.py

The duration budget can change the grid count and saved extrema in a rerun. The interval signs apply to the stated exact rational data and the normalized width family.
