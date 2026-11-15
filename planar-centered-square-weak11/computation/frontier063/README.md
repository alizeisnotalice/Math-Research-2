# Exact finite-source structural frontier: prepared, not yet executed

This new implementation follows R61's finite source-union cut certificate.
It has no new mathematical result before actual execution and review. Frozen
R57/R61 files are unchanged. Root owns all computation launches and Git work.

`exact_source_frontier.py` enumerates geometric source closures by intersections
of all coordinate-interval masks, rather than by all source subsets. Any such
intersection is its tight bbox closure, and any source-subset tight bbox occurs
as such an intersection. The same exact strict radius enclosure, labeled outer
arrangement, all-cut zeta sum and nonproportional Hminus flow from R61 are then
used. Mass and four extrema arrays have signed 64-bit entries; source input
bounds are checked. All areas and comparisons use arbitrary precision integers.
Fraction objects occur only when final values are recorded. The full N=20 cut
loop still evaluates all 1,048,575 nonempty cuts; no geometric relaxation or
random cut selection is hidden in this optimization.

`build_structural_inputs.py` creates twelve distinct structural strata, six
cases each with N=16,17,18,19,20,20. It uses no random sampling. Long connected
rows, diagonal chains, unequal two-scale clusters joined by light bridges,
dyadic shells, frames and heavy-core/light-grid geometries test noncompact and
proper high-rank cuts. These are prescribed tests, not adversarial maxima or
general theorems. The second N=20 case also perturbs geometry, not merely scale.

Root should run these gates in order, with `HERE` set to the absolute directory
containing this file and `OLD` to the frozen author R61 computation directory:

```sh
python3 "$HERE/verify_against_r061.py" --reference "$OLD/finite_source_cut_certificate.py" --saved-inputs "$OLD/r061-bounded-inputs.json" --out /Volumes/data/CodexTemp/T2/frontier063/equivalence.json
python3 "$HERE/build_structural_inputs.py" --out "$HERE/structural-inputs.json"
/usr/bin/time -l python3 "$HERE/exact_source_frontier.py" --input "$HERE/structural-inputs.json" --stratum 0 --case 0 --out /Volumes/data/CodexTemp/T2/frontier063/pilot-s00-n16.json
```

No CPU or memory benchmark has yet been measured by this implementer. Each
frontier job is one Python process and uses one core. At N=20 the five packed
signed-64-bit arrays cost 40 MiB; weight/zeta pointer tables cost another 16 MiB
before integer objects and geometry. A provisional engineering estimate is
150–300 MiB per case, but it must be replaced by the pilot's actual RSS. Default
hard address-space budget is 768 MiB per process and wall budget is 900 seconds
per case. This is a safety ceiling, not a measured resource claim. The pilot
decides concurrency; the source does not launch workers itself.

For each independent stratum s in 0,...,11, root may run:

```sh
python3 "$HERE/exact_source_frontier.py" --input "$HERE/structural-inputs.json" --stratum "$s" --out "/Volumes/data/CodexTemp/T2/frontier063/stratum-$s.json"
```

Each completed case is atomically checkpointed. The process stops on failed
assertion, timeout, memory exhaustion, bound above B, or bound above 7/4. A
failure checkpoint records the input case index; it is not a successful
certificate. `--case` and `--limit` provide bounded pilots without changing
inputs. Every final output stores source/input hashes, versions, stage timings,
RSS in platform units, exact critical mask/rank/compactness, all cut counts,
Gamma/min mass, sparse cell weights and exact flow arcs.

The optional `--deletion-rank 1`, `2` or `3` additionally checks **every** proper
compact block through that rank, keeping all old complement mass. It encloses
the full deletion loss by

    |H_full^- \ H_remainder^+| <= |H_full \ H_remainder|
        <= |H_full^+ \ H_remainder^-|.

Both bounding areas use exact integer union differences. Each block's area and
mass-normalized lower/upper ratios are retained. This can test the proposed
small-compact-block deletion route; it does not establish the existence of a
cheap block for arbitrary inputs. Deletion mode is disabled by default because
all rank-three blocks can make it much more expensive than one all-cut run.
Root should first pilot the included four-source minimal input with
`--mode deletions --deletion-rank 2`, then prioritize structural deletion
probes with `--mode deletions --deletion-rank 1` or `2`. This mode skips the
all-cut arrays and flow, recomputes every genuine remaining-source halo, and
saves remaining inner/outer areas as well as the deletion-loss bounds.

Each input also records exact source-hypergraph components and whether all
components satisfy 33 d(C)^2 <= 35 m(C), already paid below 3.5 by R53. A
large deletion score on such a terminal is not a selection-conjecture
counterexample. Failure of every rank-three candidate only excludes the
bounded-rank strengthening; a cheap higher-rank block may still exist.

Open research remains a uniform compensation mechanism for noncompact source
cuts, or a general mass-paid deletion theorem. These finite certificates do not
change C2's general bounds.

`--max-mass-compact` adds every proper maximum-mass strictly eligible block
(d(S)^2 < m(S)), including all ties. It can be used alone in deletion mode,
or alongside bounded-rank candidates. Maximum-mass eligible subsets are bbox
closures: an omitted interior positive source would increase mass and preserve
strict eligibility. Thus the complete positive-rectangle closure catalog finds
all maximizing blocks without a 2^N source-subset geometry loop. The output
separates maximum-mass candidate minima from bounded-rank minima; if the full
input itself is a maximizing eligible block, its exclusion as improper is
explicit. No maximum-tie sampling is used.
