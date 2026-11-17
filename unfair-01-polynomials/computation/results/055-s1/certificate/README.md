# Round 055, degree-8 S={1}: compact archival index

This directory freezes the sources, summaries, canonical records, source-byte ledgers and independent-comparison bindings used by the necessary-prefix/root-screen work.  `verify_archive_index.py` checks archive hashes, gzip payloads, cardinalities, stage links, root mother/H/code links and the copied independent comparison summaries.  It is **not** a new complete mathematical enumeration or a substitute for the parent portable full verifier.

The 1,140 stage-4 and stage-5--7 original per-core JSON outputs are intentionally represented by canonical records plus source-byte ledgers: copying them would add redundant small-file noise.  All 60,134 row-4 records, all 67 row-7 records and all 67 root mothers (with exact `H=P*Pbar` and failure code) are retained.  Every gzip member is below 10 MB.

Historical disclosure: an early, unarchived independent-root draft used an incorrect `N*x^2*F` term. It disagreed immediately and was replaced before any passing certificate; it is not in this package.  The archived independent root verifier uses `F^2 + T*x*F + N*x^2`.  The historical stage-4 source has an integer-key JSON-resume binding issue (JSON makes those keys strings); it is frozen byte-for-byte here, and later portable verification must use a fresh run rather than trust that cache.

Run only the archive/index check:

```sh
python3 verify_archive_index.py
```

This package uses no symlinks and does not depend on the original absolute paths for verification.  Source locations and source hashes in the ledgers are provenance information.
