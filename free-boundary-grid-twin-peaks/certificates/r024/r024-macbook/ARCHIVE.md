# r024 MacBook finite compute records

The remote host is deliberately identified only as a MacBook test machine;
no SSH configuration, host key, private path or credential is archived.

`smoke/` contains the n=3 P=257 R=2 cross-machine equality check. `n8/`
contains the completed n=8 P=16384 R=32 GPU run and strict aggregation.
`n12/` contains the independently audited completed n=12 P=131072 R=128
run. All retained sources are identified by source hashes in the sidecars and
receipts. The staged n12 audit receipt redacts only the original host label;
its raw/source hashes and audit fields are otherwise retained. The remote
n11 batch completed only after this round's cutoff; it, n16, and n24 are
excluded from this r024 package.
Finite SMC estimates are not asymptotic or certified-probability claims.
