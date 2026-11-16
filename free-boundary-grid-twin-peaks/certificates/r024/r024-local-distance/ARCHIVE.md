# r024 completed local tau=2 distance recovery

This package retains the completed n=16, P=131072, R=128, L=16 terminal
distance raw record, frozen terminal source, strict aggregator, sidecar,
reports and exit telemetry. `completed-n12-distance-and-geodesic-r024.json`
also records independent completed-batch checks for the earlier n=12 distance
and n=12/n=16 geodesic runs, but their r023 raw material is not duplicated
here. Trace output and unfinished r024 queues are excluded.

All reported values are finite SMC ratio-of-means diagnostics for the
complete \(K=2\) law reconstructed through the tau=2 L-lift; delta SE is not a
confidence interval or mixing certificate.

Portable replay from this package root (the result is compared byte-for-byte
with `reports/n16-distance-eps-1-2.json`):

```sh
python3 -B source/aggregate_tau2_distance.py \
  raw/n16-p131072-r128-l16.jsonl \
  --sidecar receipts/n16-distance.sidecar.json \
  --source source/tau2_distance_terminal_r023.swift --epsilon 1/2
```
