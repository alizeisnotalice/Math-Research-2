# Total-51 semantic cache evidence (round 48 staging)

This package covers the finite cap=6 total-degree 51 search only.  It does not
claim a general theorem.  `cache/cache-mXX.jsonl.gz` stores one lossless record
per original cache JSON: `relativepath`, its SHA-256, and its original UTF-8
bytes.  Recreate any cache file by writing `utf8.encode("utf-8")` to
`relativepath`; verify its SHA-256 first.

All listed archive files are below 10 MiB.  `input/` contains gzip-compressed
exact input bytes.  `code/` contains the independent coverage/cache audit and
the full raw-gzip byte audit reports.  The raw `t51` gzip corpus is deliberately
not included because it is reproducible and 3.7 GiB; `raw51-bytes-report.json`
records that it was fully read locally.

To validate staged cache archives, decompress each JSONL, verify every record
SHA after UTF-8 reconstruction, and compare all files to `package-manifest.json`.
The existing `code/REPRODUCE.md` gives the original audit commands and scope.
