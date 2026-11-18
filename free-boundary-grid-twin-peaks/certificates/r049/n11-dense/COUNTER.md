# Ranked-layer protected-peak counting for size eleven

This package implements the protected-peak and top-two-label upper bound
from round 47 with two dense cardinality arrays. The counted events are
unchanged. Each original peak at a protected vertex is visible inside its
region, so labels with exactly two global peaks have at most two protected
peaks. If the two protected peaks already exhaust the total, the global
maximum must be protected and the second label cannot give an additional
unprotected peak. The merger uses the three regional histograms described
in the previous protected-peak proof and the exact shuffle coefficients.
It produces an upper integer for \(A_2(n)\).

## Dense cardinality indexing

Within a fixed cardinality, order independent sets lexicographically by
row mask, using increasing numeric order within each row. Precompute the
number of feasible suffix completions for each row, previous row mask and
remaining number of occupied vertices. Also precompute the rank term for
choosing a given row mask, which sums the completion counts for all smaller
compatible masks with that same previous row and remaining size. This is
a backward row transfer. Ranking a path consists of adding its terms.

Suppose \(S\) is enumerated at cardinality \(k\), and \(S-\{v\}\) deletes
a vertex in row \(t\). The rank of this predecessor at size \(k-1\) uses the
following terms:

- Rows before \(t\) have the original masks and preceding masks with one
  less remaining occupied vertex.
- Row \(t\) uses the decreased remaining size and the modified mask.
- Row \(t+1\), if present, has the modified preceding mask and the original
  remaining size.
- Rows after \(t+1\) use the original terms verbatim.

Compute the first and last groups as cumulative sums once per set. Each
predecessor rank is then two table lookups and additions, including the two
already computed cumulative sums. Lexicographic enumeration writes the
new residues consecutively. The `--verify-ranks` gate compares each
shortcut with a direct rank evaluation and each enumeration ordinal
with its direct rank. Only two residue layers are kept, with one
64-bit residue per state and no per-entry key.

The protected mask and the interface mask use separate coordinate systems.
A region of up to 66 vertices has up to 55 protected vertices, and at most
11 unprotected points, one in each row. Those two sets are indexed
independently. For a peak set, the union of its closed neighborhoods is
tracked in both coordinate systems. Its size is the sum of the two
bit counts. The same masks give the unsafe-outside and safe-outside
populations and the edges between them for the top-two histogram. There
is no shift by the number of vertices in the whole side.

## Five-prime coverage

The fifth prime is \(2305843009213693693=2^{61}-259\). The certificate
retains the first four primes and all their witness nodes. The five-prime
certificate has 50 independently verified Lucas/Pratt nodes. The product
has 305 bits and exceeds \(65!\), including the 273-bit \(60!\) and
279-bit \(61!\). It is smaller than \(66!\), so the driver refuses exact
production for a split containing a 66-vertex side. The C++ counter can
still plan such a side or compute an individual modular residue if
invoked explicitly with adequate resources.

For \(p=2^{61}-c\), the two-fold modular multiplication uses the congruence
\(2^{61}\equiv c\pmod p\). After the first fold the value is below
\((c+1)2^{61}\). After the second it is at most \(2^{61}-1+c^2<2p\), also
for the largest correction \(c=259\), so at most one subtraction is
needed. All divisors used in the side recurrence have size at most 66.
The CRT range concerns the final nonnegative counts, each at most the
side factorial. Intermediate signed sums are accumulated modulo the
primes and need not satisfy that factorial bound.

## Size-eleven resource comparison

| Side profile | Area | Protected independent sets | Dense live payload | Older open-addressed payload | Predecessor lookups per modulus |
|---|---:|---:|---:|---:|---:|
| `55555555555` | 55 | 178,150,221 | 557,950,936 B | 2 GiB | 1,913,264,884 |
| `66666666666` | 66 | 16,884,649,135 | 47,805,806,160 B | 128 GiB | 224,308,257,215 |
| `55555566666` | 60 | 1,498,356,771 | 4,389,077,248 B | 16 GiB | 17,909,118,760 |
| `66666655555` | 61 | 2,250,165,837 | 6,589,688,840 B | 24 GiB | 27,397,073,080 |

The balanced single-turn cut minimizes the larger live table requirement
and then the combined predecessor lookup count among all 462 cuts with
widths 5 and 6 and left area 60. The comparison uses 236 representatives
under row reflection. It is a resource comparison, and does not rank the
uncomputed upper bounds on \(A_2(11)\). The exact planning data are in
`profile-resource-scan.json` and `n11-plan/plan.json`.

The two balanced payloads are about 4.088 and 6.137 GiB per modulus. This
can be approached one worker per machine. These are simultaneously live
array sizes, with allocator retention, C++ metadata and actual RSS still
to be measured in a first production modulus. Timing also requires actual
measurement. Neither region was counted in this local package. The
55/66 alternative is substantially more expensive than the balanced cut,
and both its 66-point memory requirement and modulus product prevent
production under the planned resources.

## Build and replay

Put `protected_mod.cpp`, `count_protected.py` and `pratt-61bit.json` in the
same data-volume directory. Build with Apple clang and the shared module
cache. The Python driver targets Python 3.9 or later and macOS `/usr/bin/time
-l` for receipts. It defaults to planning and requires `--execute` to
start counting.

```sh
clang++ -std=c++20 -O3 -Wall -Wextra -Werror protected_mod.cpp -o protected_mod
python3 -B count_protected.py balanced-plan 55555566666
./protected_mod 55555566666 0 --max-table-mib 4480 --progress
./protected_mod 66666655555 0 --max-table-mib 6656 --progress
```

Those last commands are examples for separate guarded workers. Prime
indices run from 0 through 4. Each successful individual worker prints
one JSON residue record on stdout and progress on stderr. For local batch
reconstruction use `python3 -B count_protected.py <fresh-dir> 55555566666
--execute --workers 1 --max-table-mib 6656 --aggregate-table-mib 6656`.
The driver checks coverage and capacity before creating the worker pool.
It independently verifies old and new layer sizes by forward row transfer
and the final regional histogram totals by factorial identities.

The small gates used 20 five-modulus jobs on both sides of two cuts,
`112` and `444445555`. Their three CRT histograms and merged upper bounds
match the prior independent open-addressed round-47 ledger exactly.
Direct-rank mode also checked all 6,186,802 and 9,569,917 predecessor
indices for the two (n=9) profiles, plus small (n=3,4) gates. The sparse
size-eleven profile `11111111112` has one protected point of degree two.
Its exact three histograms from five moduli matched elementary factorial
counts. A 6,284 MiB budget correctly refused the larger balanced side
before allocation, and the 55/66 run refused insufficient CRT coverage.
All local runs were waited and checked to have exited.

To replay the gates, first run:

```sh
python3 -B count_protected.py small-replay 112 444445555 --execute --workers 2 \
  --max-table-mib 32 --aggregate-table-mib 64
python3 -B compare_profiles.py
python3 -B audit_gates.py <round47-top-two-ledger-summary> <round48-pratt-json>
```

`build_pratt_five.py` regenerates the certificate that the driver verifies
independently. All checks remain valid if the caller performs only the
small replay and local planning.
