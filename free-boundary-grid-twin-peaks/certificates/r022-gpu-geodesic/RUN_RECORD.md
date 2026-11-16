# Actual r022 run record

Source SHA-256: `97b0cb005a85d7ec1dc04800bf1bf0cd338ff1046e096bb4991c1ff7d763a766`.
The frozen input source was fc45 SHA-256
`fc45ad181f6de297001c835a77d9be4620b944169269ed3cdb0f9c9a9a5a8629`.
The compiled binary SHA-256 is
`d87b2d1f11cb8eea89eb2873aa6c1d2837f4f2388e41c6f511b077e4f9b6080a`.

Actual device/runtime: Apple M6 Metal GPU, macOS 27.2, Apple Swift 6.4.2
targeting arm64-apple-macosx27.2.0, and Python 3.15.0.  The shared Swift
module cache was `/Volumes/data/CodexTemp/cache/swift-modules`.

| run | GPU seconds | wall seconds | independent check |
| --- | ---: | ---: | --- |
| n=3, P=17, tau=10, L=16 | 0.0067522 | 0.030361 | full prefix/group trace |
| n=3, P=257, tau=10, L=16 | 0.0172817 | 0.036433 | two groups, full trace |
| n=8, P=17, tau=65, L=4 | 0.0710545 | 0.242516 | full prefix/group trace |
| n=8, P=257, tau=65, L=4 | 0.1155364 | 0.254460 | two groups, full trace |
| n=32, P=17, tau=1025, L=4 | 2.619782 | 6.322845 | full 32-word-state trace and prefix/group replay |
| n=8, P=257, tau=2, L=4 | 0.0901705 | 0.226914 | full two-peak prefix/group trace |
| n=32, P=17, tau=2, L=4 | 2.9188303 | 5.435942 | terminal output (no large trace) |

The n=3 exhaustive oracle and all trace receipts are in `results/`.  The
n=3/n=8 same-seed comparisons against a separately compiled fc45 binary are
`results/n3-legacy-fields-receipt.json` and
`results/n8-legacy-fields-receipt.json`; every non-timing pre-r022 replicate
field matched exactly.  The P=4,194,305 pre-allocation guard rejection is in
`results/invalid-particle-receipt.json`.

The independently replayed n=8, P=257, tau=2 trace has d=2 group counts
`[0,0]`, demonstrating that the kernel does not infer a geodesic merely from
an endpoint condition or apply a one-peak-only shortcut.

All invoked GPU processes completed before this record was written.  The
large P=17 n=32 trace is retained only in this temporary computation folder,
not in a production or archive directory.
