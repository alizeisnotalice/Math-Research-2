# r031 public certificate: n=7 short-tail D4 orbits

Status: passed archival verification for the finite `n=7` model. The package
combines the nine non-center representatives `0,7,8,14,15,16,21,22,23` with
the previously audited center `24`, using D4 multiplicities totalling 49.

The replay checks 36 raw JSON/exit/resource receipt triples, nine certified
r029 layer histograms, four Pratt-certified primes (49 verified nodes), and
both direct-cofactor and Garner CRT reconstructions. The exact global QE8
checks pass at `t=2,3,4`; their details are in
`n7-short-tail-orbits-crt-audit-public.json`.

The frozen C++ source SHA-256 is
`9d757fcfd65729242d2742cc94b9c34e3e88ba23d1052b985c3b70f0b8a02d13`.
The old executable SHA-256 is retained only as a historical receipt in
`PROVENANCE.json`; no executable is shipped or reverified.

This is neither a proof for all `t` nor an asymptotic result. In particular,
finite `n=7` exactness does not settle the concentration or limit question.
See [REPLAY.md](REPLAY.md) for the portable command and [PROVENANCE.json](PROVENANCE.json)
for source hashes.
