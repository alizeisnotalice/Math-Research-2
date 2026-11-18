# Two-layer modular root counts on free-boundary grid regions

For a profile \(h=(h_y)\) in the free-boundary square, the region consists of
vertices \((x,y)\) with \(x<h_y\). It uses only as many mask bits as there are
regional vertices, at most 41, including when the global grid has 81
vertices. A safe root is an outer-band vertex whose full-grid neighbors all
belong to the region. Row widths are positive and below the grid size.
The counter accepts sizes 3 through 9, including profiles whose row width
changes by more than one. The disjointness of a family of profiles is a
separate geometric condition.

For an independent set \(S\) in a region \(H\), the number \(C_H(S)\) of
orders making all vertices of \(S\) peaks satisfies
\[
 C_H(\varnothing)=|H|!,\qquad
 C_H(S)=\frac{\sum_{v\in S}C_H(S\setminus\{v\})}{|S\cup N_H(S)|}.
\]
The recurrence conditions on the largest label in the union of the closed
neighborhoods. At a safe root \(r\), the number with \(r\) as the unique peak
is the signed sum of \(C_H(S)\) over independent sets containing \(r\), with
sign \((-1)^{|S|-1}\). This is inclusion–exclusion on the other possible
peaks. A backward row-mask dynamic program counts the cardinality layers
and precomputes viable row transitions, keeping row pairs with no vertical
intersection. The solver then streams one cardinality at a time through
two open-addressed hash tables. Each transition reads a subset key of size
one smaller from the preceding table. Contributions to the safe-root
residues are accumulated immediately.

All three moduli are certified in `pratt-61bit.json`. The Python driver
checks all 30 certificate nodes before using CRT. Each modulus is larger
than 41, so every denominator is a unit. The C++ counter precomputes the
inverses and checks their products. For modular multiplication, write
\(B=2^{61}\) and \(p=B-c\) with \(c\in\{1,31,45\}\). The product of two
residues is below \(B^2\). Replacing its high base-\(B\) digit by that
digit times \(c\) twice gives a nonnegative value at most \(B-1+c^2<2p\).
One conditional subtraction gives the representative below \(p\). The
product of the three primes exceeds \(41!\), enough to recover each exact
root count, including zeros. The driver checks the recovered total against
the regional factorial, without bounding the much larger signed
intermediate sums by that factorial.

## Build and ordinary pilot

Use a fresh data-volume directory with `root_mod.cpp`, `run_crt.py` and
`pratt-61bit.json`. Compile using Apple clang and the existing shared
module cache. Python 3.9+ standard library is sufficient. The batch wrapper
uses macOS `/usr/bin/time -l` for per-worker RSS measurements.

```sh
clang++ -std=c++20 -O3 -Wall -Wextra -Werror root_mod.cpp -o root_mod
python3 -B run_crt.py pilot-plan 654343456 345656543 --plan-only
python3 -B run_crt.py pilot-count 654343456 345656543 --workers 3 \
  --max-table-mib 1024 --aggregate-table-mib 3072 --timeout 1800
```

For separate processes, `./root_mod PROFILE INDEX --progress` emits each
completed layer on stderr and a single JSON record on stdout, with its
root vector, verified processed layer populations, predicted table sizes
and layer timings. Index 0, 1 or 2 selects the corresponding prime. The
optional `--plan` mode omits enumeration and returns the same independent
set cardinality counts and planned capacities. The per-job cap is checked
before allocating either large table. The batch wrapper independently
checks the histogram with a forward row-mask recurrence, limits concurrent
processes by the aggregate table cap, saves per-job resource logs and
exit receipts, and reconstructs the root counts after all residues arrive.
All user-supplied output directories must be fresh.

For the pilot cut the two regions are complements up to reflection, with
40 and 41 vertices. They have 63,608,332 and 80,643,260 independent sets.
At the peak, the two live tables require 512 MiB and 1,024 MiB of array
payload per modulus, respectively. A slot holds one 64-bit key and one
64-bit residue. At most three jobs run concurrently under the default
budget of 3,072 MiB of live table payload. This budget concerns live arrays.
Actual peak RSS also depends on allocator reuse, metadata, and the process
runtime. The driver preserves that measurement separately. For a pair of
complementary regions it multiplies the two total safe-root counts by the
binomial number of interleavings and reports the exact individual cut
counts. It does not add mirrored events or infer an original-law
conditional probability without \(A_2\).

## Local tests completed

- All three moduli and CRT were checked using eight fresh exact U192 root
  vectors across sizes 3 through 8 and one sparse size-9 profile whose sole
  safe root is global vertex 72. The latter has 9 rooted orders and all
  other root counts zero. Its reference is the frozen independent subset
  completion gate. See `gate-audit.json`.
- Thirty separate balanced size-9 region plans matched every cardinality
  count and live table payload in the independently frozen forward
  planner. The 41-vertex pilot side was refused before allocation with a
  1,023 MiB job budget. Balanced size-9 counting was not run locally.
- A full size-8 residue for `65454345` ran 10,610,055 independent sets and
  matched every coordinate of the earlier U192 root vector, modulo the
  second prime. The forced maximum live table payload was 128 MiB.
  Measurements were 6.03 wall seconds, 2.91 user seconds and 225,869,824
  bytes maximum RSS. These are local size-8 timings, not a size-9 benchmark.
- The production driver completed and waited nine jobs in the mixed-size
  small batch. A separate one-second deadline on the size-8 full task
  exercised cleanup by process group. All processes started for these
  tests were waited or checked to have exited.

To repeat the local comparison, additionally build `u192_reference.cpp`
into `u192_reference`, prepare `small-driver-gate` with profiles `111111112`,
`34343`, and `2323`, and save an actual `root_mod 65454345 1` JSON record to
`n8-large-residue.json`. The three supplied `reference-*.json` files are the
frozen independent plan, small compact-mask and large root-vector evidence.
Then run `python3 -B audit_gates.py`. It uses the newly built binaries and
regenerates `gate-audit.json`.

The frozen coordination directory `r046-n9-gates` was checked byte for byte
against the copied prime certificate, sparse gate and n9 resource ledger.
Their original filenames are `n9-pratt-61bit.json`, `layered-crt-gate.json`
and `n9-layer-plan.json`.
