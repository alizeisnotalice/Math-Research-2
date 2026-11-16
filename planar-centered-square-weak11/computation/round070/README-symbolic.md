# R70 root symbolic cross-check

Wolfram 15.0.0 on MacOSX-ARM64 actually checked 14 exact scalar statements for normalized L=1: the rationalized derivative, square lower-bound identity, both global-minimum branches, inverse strict inequalities, the radical-to-polynomial receiver comparison, and the clamped mass-domain endpoint. Every result is `true`. These checks support, but do not replace, the separate geometric source-union proof. They certify no area estimate or general C2 value.

The original process PID 45246 exited 0; its absence was checked after wait. Kernel time was 0.065888 seconds and whole-process time 1.0232760829967447 seconds. No MATLAB was launched for this root computation. The source and actual output are byte-bound by the manifest.

Reproduce by copying `computation/round070/r070-captured-mass-symbolic.wl` into a fresh directory under `/Volumes/data/CodexTemp/T2`, then run `/usr/local/bin/wolframscript -file /absolute/fresh/path/r070-captured-mass-symbolic.wl`. The script writes its JSON beside itself and exits; inspect all 14 booleans rather than only the printed PASS. No floating values enter the checks.
