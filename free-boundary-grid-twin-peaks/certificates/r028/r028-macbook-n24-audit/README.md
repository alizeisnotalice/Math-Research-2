# n=24 returned-batch integrity and discrepancy audit

This is a finite audit of one completed 32-replicate n=24 batch.  The raw
JSONL is losslessly gzip-compressed; `raw/RAW_GZIP_SHA256.txt` records both
the compressed and uncompressed SHA-256 values and original size.  The stored
integer payload, IDs, seeds, frozen aggregation, and Fraction cone arithmetic
all replay exactly.

This is not a rank-by-rank trajectory replay: it audits only the saved
terminal group output and its frozen aggregate.

The result is deliberately negative about numerical interpretation: the
stored plug-in violates the necessary cone and has near-one weight
concentration.  Conditional on the exact cone, the included Fraction
certificate proves

`max(|B/A_hat - B/A_true|, TV(q_hat,q_true)) > 11/50`.

This is neither a confidence interval nor a mixing, implementation-bug, or
asymptotic conclusion.  It says this batch cannot simultaneously provide
accurate values for both displayed target quantities under that necessary
condition.
