# Section 8 local Metal experiment

The profiles are piecewise linear, range in [1/2, 3/2], have f(0)=3/4 and f(1)=1, and equal 1 outside a bounded interval. They represent single directed vertical intervals for the heat potential at every stage. Four runs use fixed Lipschitz caps 1, 1, 3/4, and 3/2 and knot spacings 1/4, 1/8, 1/4, and 1/4. The modulus is uniform within each run. The radius grows through 2, 3, 4, 6, 8, and 10. Each stage refits all previous real constraints along with newly added positive and negative parameters. The fitter imposes a stagewise sup-norm update budget and records the measured changes of every previously fitted batch. The candidates are feasible numerical iterates subject to reported floating arithmetic tolerances. Stages that reach the iteration limit supply feasible residuals only.

On the last 3/2-Lipschitz candidate, 56 fitted points have maximum absolute residual 7.591e-7 in the fit quadrature. Held-out real tilts at 7.680591 and -7.680598 have residuals +3.970e-4 and +3.736e-4 in independent double integration. Separate interval integration proves residuals below -0.01 at both t=9.361495 and t=-9.380263 for the shader input rounded to float32. The two pinned values are binary-exact. A common analytic estimate proves every profile in the fixed-support class has a strictly positive residual on both far wings for |t| >= 16R+18. Its uniform threshold grows with R.

Pure imaginary and mixed parameters were all outside the fit set. The complex response reported in the validation JSON is the transform residual multiplied by exp(-Im(z)^2/2), which preserves the zero criterion. The double audit also evaluates the directed heat potential U at y=0, ±1/2, ±1, ±2, ±3. For the last 3/2-Lipschitz candidate the largest absolute U over those held-out arguments is 1.285e-5. Its smoothed Q at t=1.0941689014434814, y=-0.4735383987426758, and epsilon=1/32 is negative. An interval Riemann enclosure gives [-0.086715,-0.073380] with outward decimal rounding. A negative Q is a nontriviality diagnostic irrespective of whether the transform vanishes. Here the independent residual tests disqualify the fitted candidate as an exact solution. Existence for sequences with increasing support and a common local modulus remains open.

The Apple M6 production process ran from 2026-11-18 06:52:07.602 to 06:52:54.605 UTC. It made 4,622,303,232 evaluations of 24 fitted graphs over nine modes. This is a computation count. The command-buffer timestamps total 44.3403991 s in a 47.0028399 s run. The system IOAccelerator utilization median is 97%, maximum 99%, over 89 in-process samples. The device counter includes other activity. Process 56542 exited with status zero and was absent during the recorded check.

The CPU audit checks 216 saved extrema. Maximum discrepancies are 3.20e-8 for held-out real residuals, 1.015e-5 for scaled pure imaginary residuals, 2.08e-7 for scaled mixed residuals, and 2.87e-7 for smoothed Q. For the most extreme low-mass tails, eight-node quadrature errors reach 0.183 in logarithmic mass. The tested positive-to-negative margins are orders of magnitude larger and the class's far-tail sign follows independently from the analytic bound.

## Reproduction

The archived binary profiles are the actual Metal inputs. They can be regenerated from the fit script. Run the commands below on an Apple Silicon Mac with the existing Python, NumPy, SciPy and mpmath environment. The launch script uses the shared Swift module cache on /Volumes/data.

    TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python fit.py 6
    ./run.sh 0 8 smoke
    ./run.sh 47 64 main
    TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python audit.py scan-main.json audit-main.json ioreg-main.txt
    TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python certify_q.py

The runtime budget can change the pseudorandom validation extrema in a rerun. The enclosures and reported extrema apply to the saved production inputs and results.
