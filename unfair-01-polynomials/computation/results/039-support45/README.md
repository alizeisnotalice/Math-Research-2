# A30 total-degree-45 independent replay package

This package proves a finite support-layer statement only:

> Both supports have at most six points and their degrees sum to 45.
> The singleton-unit classification leaves no unresolved support pair.

It does not certify degree total 46, and it does not turn the rebooted remote
45/46 parent controller into a completed 45/46 run.

## Actual result

| raw supports | canonical representatives | excluded | fair | unknown |
| ---: | ---: | ---: | ---: | ---: |
| 606,304,193 | 303,168,708 | 301,344,396 | 1,824,312 | 0 |

The completed semantic audit was four local Python 3.13.18 worker processes,
with a restart-safe cache.  Its runner PID 20390 was absent when
`total45-exit-receipt.json` was written.  The actual elapsed time was
1668.603 seconds.

## Files and bindings

`total45-full-semantic-report.json.gz` is the gzip-compressed 33 MiB complete
per-chunk report.  It is 3,483,614 bytes with SHA-256
`a26427cf90b5b90f618e33ec71f7ed0b47a93514e8098362689120e56f25904f`.
Decompressing it gives SHA-256
`0c23d3d74f62271d3f8cd51f76870b4570a94c292ffed98b7ac81de01ad8f957`.

`total45-chunk-bindings.json.gz` is a deterministic 3,177,536-byte index of
all 35,442 chunks: `(m,n,I-index,I,counts,gzip SHA, uncompressed SHA, cache
SHA)`.  Its SHA-256 is
`ded359e2671a9349ddfbbacd64eeb23c2695982c146fd0f7278320b42b298571`.
Each compressed file is below the project 10 MiB archive limit.

`total45-degree-pairs.json` contains the 22 degree-pair aggregates.  It has
SHA-256 `9d5cc8668f2092b5d57ac3bd5322eacaaca0b5c0e8c2276a6b37aaa8d62b9d8a`.

The source contract is fixed by these hashes:

| source | SHA-256 |
| --- | --- |
| `sources/frozen_support.py` | `d68f1010c61ad1d911749742fef8c790c0a8064bf34512b402e7e65e4703e710` |
| `sources/remote.py` | `6fe5ea85657d6c9ed49b09b3410b915f62e14ede7d6b9864745b59d85ea106a6` |
| `sources/audited_chunk_verifier.py` | `a557bf120c6bed03adbd8b2ce153a88bf76e678ed94c6bc8815d14dc33edbe4b` |
| `sources/verify_total45.py` | `cb7cb83d095a0dc6a7a4143f6d12d09c0fe0e39e6a10695e4cc9c3896c80c784` |

## Reproduction

The full temporary streamed ledger remains at
`/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round038-macbook-support46/run`.
To reproduce the semantic audit, copy `sources/verify_total45.py` and
`sources/audited_chunk_verifier.py` into one fresh directory, then run:

```sh
/Users/zzh/.venvs/main/bin/python verify_total45.py \
  /Volumes/data/CodexTemp/T2/unfair-01-polynomials/round038-macbook-support46/run \
  --out total45-report.json --workers 4 --cache total45-cache
```

This enumerates every `J` for every total-45 `I` itself; it does not import
the searcher classification.  Its cache accepts a previous chunk only after
rechecking current gzip bytes as well as frozen/remote/events/auditor hashes.
For the adversarial regression suite, run `test_verify_total45.py`; it passed
the independent twenty-row smoke audit and rejected a deleted or duplicate
support index, a witness mutation after both hashes were rebound, and a stale
cache gzip binding.

## Parent full-content aggregate check

The parent independently decompressed all 1,170,786,246 compressed bytes (16,440,295,671 raw bytes), checked all gzip/raw hashes, 303,168,708 row counts and the exact combinatorial I/J coverage. The archived checker was actually run and exited 0 in 30.129 seconds; see parent-content-hash-and-coverage.json. This aggregate check supplements the full semantic replay, not replaces it.

```sh
python check_total45_coverage.py --ledger-root <full streamed run directory> --record <data temporary result path>
```

The compressed original events.jsonl is also included; its raw SHA is bound in both the semantic and aggregate receipts. Complete reproducible ledgers remain in CodexTemp rather than Git.
