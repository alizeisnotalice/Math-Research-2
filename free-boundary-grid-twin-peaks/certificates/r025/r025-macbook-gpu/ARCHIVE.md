# Completed MacBook terminal distance batches

The package contains two completed finite `tau=2` terminal SMC batches that
were already fully terminated before this round's cutoff.  The n=16 batch
has `P=262144,R=128,L=16,multi`; the n=24 batch has
`P=131072,R=64,L=16,multi`.  Each sidecar pins the frozen r023 source and
expected metadata.  Exit and resource records show normal completion.

The n=16 raw JSONL is 13.8 MB uncompressed, so only its existing 3.41 MB
gzip representation is archived.  Reproduce it on the data volume with:

```sh
gzip -dc n16-distance-p262144-r128-l16.jsonl.gz > /Volumes/data/CodexTemp/T2/replay-n16-distance.jsonl
```

The n=24 raw is below 10 MiB and is retained directly.  The supplied summary
files are finite ratio-of-means diagnostics, not confidence intervals,
mixing guarantees, or asymptotic evidence.  They must not be pooled with
the distinct local `strat` batch solely because they target the same
quantity.

For the n=24 batch, the normalizer--distance consistency audit in the
sibling directory finds that its plug-in balance normalizer and distance
histogram violate a necessary exact first-merge inequality.  This constrains
the interpretation of its Q/tail diagnostics; it does not by itself identify
an implementation or mixing failure.
