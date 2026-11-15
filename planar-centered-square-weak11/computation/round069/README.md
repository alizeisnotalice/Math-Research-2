# R69: fixed central/peripheral state batch

This producer evaluates exactly six global maximum-block states from five pinned R63 inputs and the final R68 four-corner six-source and tree-barrier nine-source controls. It recomputes the complete closed-bbox catalog, all global maximum ties, eligible components, full/internal/residual halos, and every eligible label meeting the selected S. The six-source input is the final four-corner version.

For L=m(S), K=[xmax(S)−sqrt(3L)/2,xmin(S)+sqrt(3L)/2] times the corresponding y interval. The reported sets are Hfull, R=H(A\S), I=H(S), full-bbox erosion state W, tight-residual-bbox erosion state Wprime, Delta=Hfull\R, P=Hfull\(I union R), I intersect R, W intersect R, Wprime intersect R, and (Wprime intersect R)\I. Every quantity has a global, central (intersection K), and peripheral (outside K) interval. W/Wprime use all eligible intersecting labels, retain union overlap, and receive no additional complete-residual-halo discount in their state cost.

All source coordinates and masses remain rational in their original physical normalization. Multiplication of coordinates by c, the lcm of coordinate denominators, and masses by c² permits use of the pinned dyadic root/rectangle helpers. Integer square-root brackets at 64 bits beyond that grid enclose each radical; rational area endpoints are divided by c² 2^128. The denominators 5 and 3 in the controls are represented exactly. K uses a separate bracket for sqrt(3L)/2. Lower central areas clip inner sets to inner K; upper central areas clip outer sets to outer K. Lower peripheral areas remove outer K, and upper peripheral areas remove inner K. Full rectangle unions, rather than sums of label areas, determine all areas.

From the repository root, with the archived files at their listed destinations:

```sh
python3 -B planar-centered-square-weak11/computation/round069/central_peripheral_state.py \
  --archive-root . \
  --spec planar-centered-square-weak11/computation/round069/fixed-eight-job-spec.json \
  --job s00-c00-n16-S48623 \
  --out /Volumes/data/CodexTemp/T2/r069-one-replay.json
```

The whole fixed batch can be reproduced with:

```sh
python3 -B planar-centered-square-weak11/computation/round069/run_fixed_batch.py \
  --archive-root . \
  --spec planar-centered-square-weak11/computation/round069/fixed-eight-job-spec.json \
  --producer planar-centered-square-weak11/computation/round069/central_peripheral_state.py \
  --out-dir /Volumes/data/CodexTemp/T2/r069-eight-replay --workers 8
```

The runner requires mounted /Volumes/data, records actual Python/version/command, monitors each owned math child (768 MiB sampled RSS threshold and 120-second wall threshold), waits for exit, and verifies its PID is absent. The saved original execution used Python 3.13.18 at /Users/zzh/.venvs/main/bin/python3; receipt commands preserve their original absolute provenance. No MATLAB, LP, random sample, or additional family size was run.

To inspect the archived results without another geometry run:

```sh
python3 -B planar-centered-square-weak11/computation/round069/summarize_fixed_batch.py \
  --archive-root . \
  --spec planar-centered-square-weak11/computation/round069/fixed-eight-job-spec.json \
  --results planar-centered-square-weak11/computation/results/round069-state \
  --out /Volumes/data/CodexTemp/T2/r069-ledger-replay.json
```

`build_fixed_inputs.py --source-archive-root . --out-root /Volumes/data/CodexTemp/T2/r069-input-rebuild` regenerates the eight minimal inputs and spec from the pinned historical sources. Their provenance is not a new historical run. The runtime helpers are the already archived R64 `mass_polytope_variation.py` (4ef61b...) and R66 `erosion_ring_diagnostic.py` (4dc4fb...), with full hashes in the spec.

`eight-state-ledger.json` binds the saved output/receipt/input chain and records exact endpoints. Decimal tables in the computation scope are displays of these rational enclosures. All eight selected net mixed-cost intervals lie below L/5; this is finite evidence only. R68 gives general residual exclusion and central exactness; coincident numerical intervals alone do not certify equality of sets. Five selected states (from four distinct inputs) have a strictly positive central pure-mixed area, so central exactness does not identify Hfull with the internal halo I.
