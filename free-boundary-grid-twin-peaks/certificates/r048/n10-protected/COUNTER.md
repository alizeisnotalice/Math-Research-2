# Size-ten straight cut protected-peak counter

The counter extends the size-nine protected-peak method to 50 regional
vertices while retaining compact masks, pruned row transitions, and two
live cardinality tables. It returns three exact histograms by CRT. The
first tracks protected peaks, the second additionally requires that the
side maximum is unprotected, and the third requires a protected maximum
followed by an unprotected nonadjacent second label. The size-nine proofs
carry over unchanged. The merger counts all words with at most two
protected peaks, then excludes words with exactly two protected peaks
and an additional peak among the globally highest two labels. Each
protected vertex has all its full-grid neighbors on its own side. The
surviving count is a rigorous finite upper bound for \(A_2(n)\).

## Four moduli

The fourth prime is
\(2305843009213693723=2^{61}-229\). The three earlier primes and their
recursive certificates are unchanged. The four certificate roots have
40 distinct Lucas/Pratt nodes. The product has 244 bits, compared with
215 for \(50!\), and guarantees unique reconstruction of every frequency
in \([0,50!]\). The product of the first three alone has only 183 bits.

Set \(B=2^{61}\) and \(p=B-c\), with \(c=1,31,45,229\). A product of two
residues is below \(B^2\). The first fold replaces the high base-\(B\)
digit by that digit times \(c\), so the folded value is below
\((c+1)B\). After the second fold it is at most \(B-1+c^2<2p\). One
conditional subtraction yields the representative below \(p\). All
closed-neighborhood sizes used as denominators are at most 50 and smaller
than the primes. Final frequencies are reconstructed by CRT and checked
against factorial and marginal count identities. Intermediate signed
histogram contributions need only be stored modulo the primes.

## Resource plan and commands

For `5555555555`, each induced region has 50 vertices, of which 40 are
protected. The two sides are reflections and share one numerical region
profile. The independent sets of the protected region have 32,641,916
members in total and a largest cardinality layer of 6,905,186. Two adjacent
layers require at most 512 MiB of live table entries per modulus at the
implemented power-of-two capacities. Four parallel moduli would have 2 GiB
of live array payload. Allocator retention and process overhead can increase
RSS above those payload values. The straight count was run on MacBook in this round. The collected
raw residues and all five production cuts are in the round ledger.

In a fresh data-volume directory, put `protected_mod.cpp`,
`count_protected.py` and `pratt-61bit.json` side by side, and compile using
Apple clang with the shared module cache. The driver uses macOS
`/usr/bin/time -l`, Python 3.9 or later, a fresh output path and an explicit
`--execute` to start enumeration.

```sh
clang++ -std=c++20 -O3 -Wall -Wextra -Werror protected_mod.cpp -o protected_mod
python3 -B count_protected.py straight-plan 5555555555
python3 -B count_protected.py straight-count 5555555555 --execute \
  --workers 2 --max-table-mib 512 --aggregate-table-mib 1024
```

Without `--execute`, the driver verifies the C++ suffix-state counts by an
independent forward row-mask DP and writes only `plan.json`. With
`--execute` it runs the four modular jobs, saves original JSON and resource
logs, and writes the three final exact regional histograms and merged
upper count. The job budget is enforced before allocation and the aggregate
payload cap before starting the worker pool. The driver records actual RSS
and waited exit statuses. An independent numerator for the conditional lower bound is in the
quadrant package.

## Replayed gates

The full four-modular build replayed the two sides of the `112` and
`444445555` cuts, 16 jobs in all. Every histogram coordinate and both
merged upper integers matched the exact size-nine-method ledger. An 11-point
sparse size-ten side `1111111112` has a single protected vertex with two
neighbors. Its three analytic histograms also matched all four modular
calculations. The 511 MiB gate refuses the 512 MiB straight side before
allocation. Every local gate process was checked to have exited.

To repeat, first run:

```sh
python3 -B count_protected.py small-replay 112 444445555 --execute \
  --workers 2 --max-table-mib 32 --aggregate-table-mib 64
python3 -B audit_gates.py <prior-top-two-ledger-summary> <prior-pratt-certificate>
```

The prior inputs are `top-two-ledger/summary.json` and `pratt-61bit.json`
from the round-47 protected-peak package. `build_pratt_four.py` regenerates
the expanded certificate. The verifier inside `count_protected.py` checks
recursive Lucas witnesses independently of candidate generation.
