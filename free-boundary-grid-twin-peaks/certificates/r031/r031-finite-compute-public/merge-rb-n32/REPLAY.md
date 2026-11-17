# Public finite replay: n=32 merge-RB production audit

This package contains a gzip-compressed production raw stream, the frozen
`f1593...` Python auditor, the production receipt/resource record, and the
recorded independent aggregate.  It excludes the historical binary, the large
trace, and telemetry flood.

## Replay command

Run in a writable isolated work directory; decompression creates a 5.3 MiB
raw file there.

```sh
python3 source/verify_n32_public.py \
  --raw-gz raw/n32-p262144-r32-l16.jsonl.gz \
  --auditor source/auditor_frozen_r026.py \
  --derived-source source/production_r026.swift \
  --recorded reports/n32-p262144-r32-l16.independent-aggregate.json \
  --receipt receipt/n32-p262144-r32-l16.receipt.json \
  --workdir replay/work \
  --output replay/public-audit.json
```

The command checks the compressed raw after decompression against the receipt,
runs the frozen auditor, and compares n/P/R/L, ESS, shares, weighted moments,
`A2_over_B_mRB`, `Q_mRB`, legacy Q, and log-normalizer fields with the
recorded aggregate.

## Scope and interpretation

This is a finite input/implementation-invariant audit.  It gives no cone
check, exact comparison, probability claim, mixing guarantee, or SMC-accuracy
certificate.  In particular, the recorded `ESS=1.0000000003349663` and
`maxshare=0.9999999998325169` show extreme finite weight imbalance.  The
recorded `Q_mRB=0.22509415843036934` is only a finite imbalance diagnostic,
not evidence about probability, mixing, or the original question.

## Public whitelist

* `raw/n32-p262144-r32-l16.jsonl.gz` — required production raw only.
* `receipt/` and selected `reports/` — receipt, resource, finite aggregate,
  and small-P structural/replay reports.
* `source/auditor_frozen_r026.py` — SHA-256
  `f1593cc93db21d778289966c80539606b079d702fa3fd5c197394c27d4737389`.
* `source/production_r026.swift` — SHA-256
  `3a9f3c0a8c10409346c25f16f2d82463f928bfaacf7307b4a88a770ec8f69ec6`;
  it is needed for the frozen auditor's source-fingerprint check but is not
  compiled by this replay.
* `source/verify_n32_public.py`, `PROVENANCE.json`, `SHA256SUMS`, and this
  file.

The historical Swift source/binary SHA-256 values are provenance fingerprints
only; neither binary nor a claim that public users can reproduce its build is
included.
