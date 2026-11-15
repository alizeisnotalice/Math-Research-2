# R69 exact fixed heavy-core payment checks

From the topic directory with an existing Python/SymPy environment
(actual:Python3.13.18,SymPy1.14.0):

```sh
python3 -B computation/round069/check_heavy_core_payment.py \
  --input computation/inputs/r069-heavy-core-fixed-geometry.json \
  --out /Volumes/data/CodexTemp/T2/r069/heavy-core-replay.json
```

The three inputs are already fixed R68/R69 fixtures. The new container
preserves their exact mathematical data and records original file hashes;
it does not claim their different JSON serializations are byte-identical.
No old workspace is consulted by the replay script.

The checks are exact Fraction domain/branch/near-source predicates and
rational rectangle Gamma, with independent predicate-cell integration.
Every high label is tested against the closed protected K via squared
rational comparisons. Gamma always subtracts C intersect K, not all K.
The theorem is proved separately; the code does not measure full halo,W
or Wprime or invoke any LP/MATLAB/random search. The generic touching-
receiver union in the theorem is not the experimental internal halo I.
