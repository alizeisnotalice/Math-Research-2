# R79 one-batch production ledger: pending independent audit

One authorized batch ran on the existing Python 3.13.18, macOS arm64, using
the pinned source and parameters. Parent PID 15433 exited 0; helper wall
time was 3.847 seconds and script time 3.719 seconds. The configured worker
limit was eight. Five worker PIDs actually executed tile tasks:
15435,15436,15437,15438,15439. All six recorded PIDs were subsequently
observed absent. The executor context had joined the pool before writing
the summary. Idle pool-worker PIDs were not recorded, so this ledger does
not claim an observed list of eight worker PIDs.

No second quadtree batch, parameter search or automatic extra refinement
was launched. The separate 0.307-second postprocess only reads production
certificates and frozen actual-area outputs, verifies their hashes and
observed exits, and writes the comparison ledger. It does not reclassify
cells or constitute an independent certificate replay.

| Input and continuous envelope state | State-area bracket, decimal display | Gap target |
|---|---|---|
| Corrected, exact pure-old plus three mixed branches | [.010064875217462937,.010268420287157458] | Not met |
| Corrected, four soft branches | [.010205825738730889,.010412841339252006] | Not met |
| Far-face, either variant | [.01921122504438199,.019276964150105845] | Met |

Full exact rational endpoints, rather than these rounded displays, are in
the summary and certificates. Every state upper endpoint is below 1/8 and
80003/640000. These are fixed-input outer-state bounds, pending independent
theorem and certificate review; they are not a final general upper bound.
State lower endpoints are not actual halo lower bounds.

The corrected case retained 45,007 complete leaves from 60,008 created
nodes. Three of its four tiles reached the fixed node protection; its
case gaps are about .0002035451 and .0002070156, exceeding 1e-4. The far-face
case retained 23,905 leaves from 31,872 created nodes. Its case gap is about
.0000657391, within 1e-4, even though one tile did not meet the stronger
individual budget of 1/40000. The per-tile termination field and the
case-level target comparison are both retained; neither replaces the other.
Every unresolved leaf remains in its compressed certificate.

The inherited independently audited R78 actual Gamma areas are approximately
.007652366201647474 (corrected) and .016617159266386683 (far-face). The main
continuous envelope exceeds the actual area by a certified interval
approximately [.0024125090,.0026160541] and [.0025940658,.0026598049],
respectively. The corrected four-soft variant adds some further slack.
These differences compare two enclosures; no actual-area computation was
rerun and no state lower bound was substituted for actual area.

Artifacts are under `computation/results/round079-quad`:

- `production/summary.json`: actual version/PIDs/timing and both state brackets;
- `production/{corrected,farface}-certificate.json.gz`: full tile/leaf records;
- `process-receipt.json`: the unique quadtree batch receipt;
- `comparison-pending-audit.json`: exact comparison with frozen R78 areas and exits;
- `postprocess-receipt.json`: the read-only postprocess receipt.

The production source SHA256 is
`428c3c90516d9d51a062c602fefb85c8da722ae9db7abf331c1f0c5eb7645aee`;
parameters SHA256 is
`380edaf5c4fa12b9fa4cb58caab5fc3e8e5bc8d8fd0c7d7dc2fd42e0d5cfaf73`.
The source used in the run was not changed after execution. The original
R78 sources, results and manifests were not modified.
