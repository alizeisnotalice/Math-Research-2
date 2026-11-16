# Scope of the P=65536 address audit

The Metal terminal code forms prefix addresses such as `id*MAXN2+k` in 32-bit `uint` arithmetic, with `MAXN2=256`. The old entry guard `P<=UInt32.max` does not by itself prove that every such product is address-safe for every accepted P.

For the completed production run, `P=65536` and `N=144`. Therefore `id<=65535`, `k<=N-2=142` for a prefix access, and

Here `65535*256 = 16,776,960`, so the maximum is `16,777,102`, strictly below `2^32`. The allocated prefix extent is `P*256=16,777,216` entries, so this run's relevant `uint` product and its ordinary buffer index are safe.

The source also has the separate m2/m3 `UInt64` count bounds documented in `RANGE_AND_INDEX_AUDIT.md`. This file adds only the production-P address calculation. It does not repair or certify the old upper P range generally.
