# Finite certification report

All 40 raw residues passed the portable exact audit. The four certified primes
have product greater than `72*48!`; `h` and each far counter are reconstructed
in `[0,48!]`, and `E` in `[0,72*48!]`. Direct CRT and Garner agree field by
field. The global h vector agrees with public r033, and exact first-five
contributions `sum(t h_t)/A2` and `sum(t E_t)/(49 A2)` are in the audit JSON.

No L>=6 statement is inferred: this package does not use a full Q/tail integer.
Far counters are checked only for source-implied nesting; R45 is recorded per
orbit rather than assumed zero.

Separated residual audit: the independently stored exact full squared-Euclidean
numerator is `2313395312707855007943501445611841392991226454185634816`.
Subtracting the exact first-five spatial numerator
`1598021776624057449634826922769523822439279480787325952` gives certified
L>=6 residual `715373536083797558308674522842317570551946973398308864`.
No full exact Manhattan-threshold numerator for `d1>7/2`, `d1>7`, or
`d1>21/2` was located, so this package deliberately does not report their
residual probabilities. The subsequent finite fibre audit found that all three
events are constant on every squared-Euclidean difference fibre and that the
r018 `distance_counts` source is explicitly Manhattan; exact far residuals are
now recorded in `evidence/n7-spatial-residuals-audit.json`.
