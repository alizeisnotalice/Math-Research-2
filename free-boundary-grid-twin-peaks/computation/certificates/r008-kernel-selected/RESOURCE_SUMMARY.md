# Resource and source record

All reported mathematical signs are exact-integer signs, not floating-point
tests.  Decimal values are explanatory only.

| run | status | observed resources |
|---|---|---|
| n=7, 30 sequential D4 orbit representatives | all exit 0; all satisfy the auxiliary candidate | 249.42 s wall, 433.07 s user, maximum RSS 3,052,814,336 B; driver time log was summarized into n7-critical-d4/compact-manifest.json |
| n=8 extended-L original certificate | exit 0; strict positive gap | 1.77 s wall, 2.93 s user, maximum RSS 374,603,776 B, peak footprint 374,194,824 B |
| n=8 frozen-source rerun | byte-identical JSON to the original certificate | 1.76 s wall, 2.97 s user, maximum RSS 374,554,624 B, peak footprint 374,145,672 B |

The source files are frozen in this package.  The rerun was compiled with
c++ -O3 -std=c++20 -pthread; current local tool inspection reported Apple
clang 21.0.0 and Python 3.15.0.  The rerun executable had SHA-256
259ee8be23912a25678486db8064ba82817e2ab742057fdb47c3cdf34fe60b38;
executables are intentionally excluded from the package.  The rerun JSON
SHA-256 was
6d082139b0f1f778d125e9ee40ce096a8ad09a920bc9531321f52fa462d5cfec,
matching the saved result.

The original n=8 source/output relationship is supported both by its saved
run record and by this new frozen-source rerun.  The package does not claim
that a binary hash alone proves a source correspondence.
