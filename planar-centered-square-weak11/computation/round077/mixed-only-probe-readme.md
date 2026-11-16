# R77 mixed-only finite screen

This single exploratory batch tested the candidate that the union of the three receiver labels containing the old atom and at least one of two outside atoms, after subtracting the complete outside halo, has southwest-orthant area at most half the old mass. The old singleton is not added as a fourth full label. Intersections of the mixed receivers with the old singleton are still included; this is NOT a subtraction of the old halo.

There are1536 fixed seeds starting770000. The old mass is1, outside masses are the dyadic choices specified in the final source, and source positions are rational multiples of1/20. Input mass and strict subset eligibility are checked with Fraction, while square roots and all area arithmetic are floating. The search covers neither all positions nor all masses. All three mixed labels and all three residual labels are retained. It imposes no complete-grid capacity or effective-source-window restriction.

No value exceeded0.50001. The largest floating result is about0.48965808745, at seed770712. This negative batch is not a proof. Subsequent analytic work on a nearby configuration supplies a counterexample in a separate component; it does not retrospectively turn this batch into a positive search result. The record preserves exactly the top20 and counts that this run emitted. The entire batch can be regenerated from its fixed source and seeds.

Actual Python3.13.18 execution used parent51264 and workers51266,51267,51268. The pool joined, exit0 and all four PID absences were verified. Internal time about0.135s, wrapper time about0.197s. No MATLAB was launched. A first wrapper command had a Python syntax error before any Popen could execute; the probe itself was launched exactly once after that wrapper was corrected. That command-assembly error is not counted as a failed mathematical run.

From the topic directory, replay with:

```sh
python3 -B computation/round077/probe_mixed_only.py --geometry computation/round072/actual_halo_geometry.py --output /Volumes/data/CodexTemp/T2/r077-mixed-only-screen.json --count 1536
```

The receipt includes original execution paths and PIDs as provenance. This component proves no mathematical inequality and changes no general C2 bound.
