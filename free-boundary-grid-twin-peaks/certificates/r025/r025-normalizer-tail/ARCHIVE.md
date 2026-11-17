# Top-two normalizer--distance consistency check

For the true two-peak law, the exact first-merge identities imply

```
B/A2 <= sum_{d>=2} d/(2(d-1)) P2(rho=d),
P2(rho>r) <= 2r/(r-1) (1-B/A2),  r an integer >=2.
```

The proof is pointwise in the first merge, so it requires no independence
between branch balance and root distance.  Exhaustive 3 by 3 enumeration
checks `A2=143808`, `B=81984`, `B/A2=61/107`, the exact first-merge identity,
and every listed distance inequality.

When the two completed n=24 finite SMC outputs are inserted separately, both
violate the exact linear constraint: the local `strat` slack is about
`-0.23139` and the MacBook `multi` slack is about `-0.27946`.  At r=12,
their plug-in tails also exceed the corresponding necessary upper bounds.
This does not on its own identify a coding error, prove statistical
significance, or establish a mixing failure.  It rigorously says that each
collection of finite plug-in ratios cannot all be treated simultaneously as
exact values of the original K=2 law.  Their Q/tail outputs must therefore
not be used for trend claims.

Using exact rational normalization of the archived decimal summaries, the
certificate further proves that at least one of the absolute `B/A` error and
the total-variation error of the distance law exceeds 0.15 for local `strat`,
and exceeds 0.18 for MacBook `multi`.  These are deterministic feasibility
floors for the displayed plug-in pair, not confidence intervals or a bug
diagnosis.
