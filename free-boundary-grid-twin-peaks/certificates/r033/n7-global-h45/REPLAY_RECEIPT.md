# Actual audit receipt

The packaged Python auditor was actually run against this package's 120 raw
evidence files: 40 residues, 40 exit receipts, and 40 resource receipts. It
passed all asserted checks. It uses CPython's arbitrary-precision integers and
no random input, GPU, SSH, or production enumeration.

It checked each CRT coordinate both by cofactor CRT and incremental Garner
CRT. The product of the four certified primes exceeds `48!`, and every
reconstructed coordinate is nonnegative and at most `48!`. The source also
passed a separately compiled n=3 centre smoke with residues
`2688,1600,1024,640,384`.

The finite result has global R45 exactly zero and finite QE8 checks true for
`T2` through `T6`. Neither statement is an asymptotic or all-dimension claim.
