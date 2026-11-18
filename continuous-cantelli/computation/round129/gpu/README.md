# Section 8 noncompact-wing Metal calibration

The four saved terminal cores of the previous Section 8 correction experiment are kept on [-10,10]. On each wing, at radius r>=10, the offset from the constant coefficient is a (10/r)^p (1-10/r). It joins continuously at the endpoints. The amplitudes have either sign, with magnitudes from 0.025 to 0.22, and powers range from 0.55 to 3.5. The two original anchors remain f(0)=3/4 and f(1)=1. Each added wing has Lipschitz constant at most 0.176. The earlier real-parameter fitting was done on the cores before the wings were added.

Metal scans nine groups per core. It evaluates positive and negative real holdouts, pure imaginary and mixed complex holdouts, and large real tilts in subcritical, critical and supercritical powers. For a complex argument z, the saved residual of the transform is multiplied by exp(-Im(z)^2/2). Large real tests use log F after completing the Gaussian square. The production shader truncates centered favorable-wing quadrature to ±20 standard deviations, retains the core for |t|<64 and omits the opposing wing in large real tests. The independent CPU audit integrates both complete wings and the core.

Among 36 retained extrema the maximum absolute CPU-GPU discrepancy was 1.682e-8 for real holdouts, 8.813e-6 for scaled imaginary holdouts, 4.878e-7 for scaled mixed holdouts, and 4.797e-5 for the large-tilt logarithm. All 20 retained large-tilt signs agree with the favorable-wing asymptotic prediction. At t=7.978767, the real residual is 0.006345785, of which 0.006286868 comes from the added wings. At t=-7.979012 in a different retained graph, the corresponding values are 0.006366383 and 0.006307705. Near the largest pure imaginary frequency the wings have much smaller influence because the Gaussian is centered at zero.

The asymptotic theorem in tail-limit.md proves a nonzero favorable-wing coefficient excludes exact normalization in this family. For power p<2 the log of the transform grows in magnitude at scale |t|^(2-p). At p=2 the transform tends to exp(-a·10²/2), and at p>2 the transform residual decays with scale |t|^(2-p) and a nonzero leading coefficient. The proof covers any locally positive compact core.

Apple M6 production started at 07:18:29.473 UTC on 2026-11-18 and ended at 07:19:59.475. The planned 07:12:00 UTC gate was already past when the wrapper's prelaunch sample was recorded at 07:18:28. The experiment yielded 5,687,746,560 evaluations, with 67.791114 s of Metal command time over 90.002083 s elapsed. IOAccelerator readings had median 87% and maximum 97% over 166 in-process samples. This device counter is systemwide. PID 64141 exited with status zero and was absent in the recorded exit check. The peer MacBook GPU lifecycle enveloped the whole local run. The subsequent local CPU worker launch was at 07:20:40.276 UTC, after this run ended.

## Reproduction

Run on Apple Silicon using the existing Python, NumPy and SciPy installation. The shared compiler module cache is on /Volumes/data. The saved profiles.bin already contains the actual input.

    /Users/zzh/.venvs/main/bin/python make_profiles.py
    ./run.sh 0 8 smoke
    ./run.sh 90 160 main
    TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD" PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python verify.py scan-main.json audit-main.json ioreg-main.txt

The launch may be gated by setting SCAN_START_EPOCH on the production invocation. An already elapsed epoch permits immediate launch. The time budget determines how many samples are evaluated and can change the extrema in a rerun.
