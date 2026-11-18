# R101 local Metal calibration of Gaussian windows

The shader integrates 12 continuous trapezoidal dips in a constant background against the full Gaussian window and the unit window. Parameters vary across geometric spacings, exponential width powers (1, 3/2, 2), depth, and window offsets. Log-domain accumulation retains widths too small to represent as ordinary float32 values. Each potentially contributing dip is split into a plateau and two sloping sides. The sides are integrated on a scale set by their exponential decay. The dispatch grid uses deterministic hashed parameters and reports evaluations rather than a count of distinct profiles.

Production ran on Apple M6 under macOS 27.2 (26B123), from 2026-11-18 06:01:48.046 to 06:03:38.047 UTC. It covered 38,967 batches of 131,072 evaluations each, for 5,107,482,624 evaluations. Metal command timestamps summed to 94.9994997 s of the 110.001684 s production interval. The IOAccelerator device utilization counter median was 92% (206 in-process samples, range 0–95%). That counter is systemwide. The command timestamps belong to this executable. PID 38472 exited with status 0 and was absent in the subsequent process check.

The independent CPU audit covers all nine saved representatives. The largest absolute difference in a logged integral divided by lambda squared was 1.621e-8. The strongest positive-full, negative-unit representative gave CPU logarithms 131257400.8259 and -131291753.7365. The closest-to-zero log reported by the GPU was 0. The CPU recomputation gave -3.92165e-7 for that sample. The scan counts budget and helps locate robust cases. A rounded near-zero value does not certify a sign. Unit-window boundaries can also move relative to extremely thin dips if float32 center rounding is comparable to the distance to a window edge.

The separate fixed-profile proof gives full rates 9/400 and -9/128 along two infinite sequences, with unit-window rate -1/8 along the first. At lambda=320, the full integral has log greater than 2295 while the unit integral is below 1. At lambda=512, the full integral is less than 6 exp(-15000). The construction proves these statements for the stated continuous profile. Its realization by a Cantelli graph remains open. The fixed-profile controls in the JSON are from the CPU calculation.

## Reproduction

Run from this folder on an Apple Silicon Mac with the command-line developer tools. The launch script uses the shared Swift module cache on /Volumes/data. The CPU audit uses the existing Python environment with SciPy.

    ./run.sh 0 24 smoke
    ./run.sh 110 40 main
    PYTHONDONTWRITEBYTECODE=1 /Users/zzh/.venvs/main/bin/python verify.py window-main.json audit-main.json --ioreg ioreg-main.txt
    /Users/zzh/.venvs/main/bin/python certificate.py

The production duration is a stopping budget. Timestamp, evaluation count, and extreme samples can vary when the command is rerun. Hashes document the actual saved files.
