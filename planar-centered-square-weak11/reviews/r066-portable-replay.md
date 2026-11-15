# R66 portable archive replay supplement

2026-11-15. This adds a portable execution entry to the frozen independent R66 audit. It does not modify or replace the original seven-file audit manifest, either independent geometry program, any producer source or any proof pin.

From a clean archive checkout, run with Python3.10 or later, using only its standard library:

```sh
python3 -B planar-centered-square-weak11/computation/review_r066_portable_replay.py \
  --output /Volumes/data/CodexTemp/T2/r066-portable-new/result.json \
  --workers 12
```

`--output` is mandatory, names a new `.json` file below `/Volumes/data/CodexTemp`, and must not already exist. A sibling `.receipt.json` is written. The data volume must be mounted. `--archive-root` optionally names either the repository root or its `planar-centered-square-weak11` directory; the default is the repository root inferred from this entry's location. Workers may be set from1 to12. Run without Python `-O` so the frozen core's exact assertions remain active.

The entry embeds and checks SHA256 pins for32 archive-relative dependencies: the frozen independent R66 grid core, its frozen independent R64 geometry helper, and15 budget/erosion JSON pairs. These pairs cover R64 final, all12 R63 strata, the five-source thin control and the separate three-source fat control. The locations are the actual compute integration manifest destinations: `computation/results/round066/theta` and `computation/results/round066-erosion`. After checking hashes, producer-source metadata and pair linkage, the entry calls only the frozen grid core's `init_worker` and `evaluate_case` interface. It does not invoke producer geometry, historical source scripts, or old `.local` input paths. The frozen modules contain historical path constants, but these are not used by the imported core functions in this entry. There is no fallback to another workspace.

The actual replay on2026-11-15 read the root integration workspace's staged archive, not the original production locations. Python3.13.18, twelve workers,4.79056seconds, exit0. It replayed75 input records/202 proper maximum blocks,1106 individual high labels, and4,296,100 compressed grid cells. Every exact saved budget/halo/ring interval and classification matched. The original R63 subset remains199 budget passes,181 stronger-condition passes and18 strict failures. The new unique-fat-maximum three-source control has B_lower>3.5 and Δ_upper<1. The original201-block result remains unchanged; this run combines that coverage with the extra fat control.

The entry was executed from the review clone with `--archive-root` pointing at the staged clean integration clone; its result records the exact entry hash, all relative dependency hashes and the invocation. The owned entry PID97193 was checked absent after exit0. The archived `round066-portable-replay` result, runtime receipt and process receipt preserve this actual execution; no historical receipt has been relabeled as a new run. Output hashes vary when a replay is rerun because runtime metadata changes.

This is an exact finite archive replay, not an additional mathematical theorem. The fixed θ=3/4 filled-container envelope selection is refuted by the saved three-source example. The original true-Δ selection problem and general C2 bound remain open.
