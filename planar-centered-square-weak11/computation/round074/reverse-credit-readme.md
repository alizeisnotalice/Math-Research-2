# R74 reverse residual credit: exploratory probe

This tests a separate candidate statement: for a finite positive atomic measure ν supported on Q=[−1/2,1/2]² and M_square ν(0)≤1, is ν(Q)≤4 |H(ν)∩Q|? This round proves no such theorem.

The fixed batch has 200 cases at each source count 2,3,4,6,8,12. Positions lie on the 1/20 grid; mode, integer weight and seed rules are in the final script. The common mass scaling is exact Fraction arithmetic over every complete tied distance prefix at the origin, so the origin is outside the strict halo. All closed bounding-box closures are enumerated. Strict eligibility is tested rationally, but receiver rectangle endpoints and clipped union area use floating point. No ratio exceeded 4+10^(−6); the displayed numerical maximum 4.000000000000001 is attained by opposite corner atoms. Negative search is not evidence of a general upper theorem.

Python 3.13.18 actually ran the 1200 cases with three worker processes in about 0.45 s. The original run receipt records exit 0, worker join and PID absence. No MATLAB was launched. The 1200-row raw output is reproducible and retained in the local project archive; the public summary retains top 20, per-count maxima, source and raw-output hashes. The old stdout key `number_above_4_plus_1e6` is a label typo: the source expression and the public summary both use 4+1e−6. The executed source is retained byte-for-byte.

From the topic directory, reproduce with an output path on the data volume:

```sh
/Users/zzh/.venvs/main/bin/python3 -B computation/round074/probe_reverse_credit.py --geometry computation/round072/actual_halo_geometry.py --output /Volumes/data/CodexTemp/T2/r074-reverse-credit-replay.json --workers 3 --per-count 200
```

The receipt and result contain original execution paths and PIDs as provenance, not portable path assumptions. The raw timing/PID fields will differ on replay. This helper does not test the full four-corner deletion theorem and has no effect on the general C₂ interval.
