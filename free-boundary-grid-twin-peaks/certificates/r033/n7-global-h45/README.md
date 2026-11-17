# Exact n=7 global h1--h5 audit

This is a finite, all-49-hole D4 reconstruction from 40 saved modular
checkpoints: the ten D4 representatives `0,7,8,14,15,16,21,22,23,24`, four
Pratt-certified primes per representative. The centre p0 is the canonical
foreground retry; an earlier failed nohup attempt is not present. The v22
placeholder in a different source package is excluded: canonical v22/v23
records are the eight MacBook receipts included here.

The independent audit obtains

```
h4 = 3396771146185020250675732839307818411684105820449792
h5 = 2011907760585005003686570084046477536187007211969536
h3 - 3*h4 + 2*h5 = 0.
```

It also reproduces r31's global and per-orbit h1--h3 exactly. QE8 holds for
the finite global `T2,...,T6` checks. These facts are finite `n=7` statements,
not an asymptotic result or a proof of any all-n tail principle.

No historical executable is distributed. The `source/` C++ file has SHA-256
`08c31508f6d3cb0ae50f74973aa372a32e3bac8439d9d6bfd51d4bf6c0312fc8`.
