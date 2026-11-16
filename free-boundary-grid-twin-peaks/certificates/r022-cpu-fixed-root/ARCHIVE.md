# r021 CPU fixed-root certificate: corrected n=11 staircase

Source snapshot: `a022-cpu-fixed-root-r021`.  The retained production data are
exactly the eleven roots `(0,0),(1,0),(1,1),(2,1),(2,2),(3,2),(3,3),
(4,3),(4,4),(5,4),(5,5)`, i.e. IDs `0,1,12,13,24,25,36,37,48,49,60`.
Each has P=65536, R=64, L=16 and its CSV plus metadata are included.

The wrongly selected off-stair IDs `23,34,35,45,46,56` are intentionally
excluded from numerical aggregation; `OFF_STAIR_NOTICE.md` preserves the
minimal correction record.  Results are finite SMC diagnostics, not exact
root counts or proof of either candidate inequality.  To reaggregate without
running SMC: `python3 -B summarize_n11_roots.py . > replay-summary.json`.
`fixed_root_smc.cpp` is retained, but its historical binary is not.
Final source SHA-256: `b6fc769e5f496146b9473b8dc66e2ece122cd4a18cd329f0d1bb845c6dbd84bd`; replayed summary SHA-256: `544988b5109b29ebb9921befa1b80db57f0717c50d48e35595ebf9c84aaec3e4`.

The earlier `n11-root-summary.stdout` was omitted during final integration: it predates the corrected coordinate-specific candidate fields in the authoritative JSON.
