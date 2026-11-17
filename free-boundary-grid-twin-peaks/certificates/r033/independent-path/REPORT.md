# Independent audit: first-merge path-prefix code (r033)

## Verdict

For the top-two `B` slice actually supplied to this audit, the proposed path
move has a valid explicit inverse code. The two displayed code-space
identities are correct:

\[
\sum_{\ell=3}^{q+1}\binom{q-1}{\ell-2}(2^{\ell-2}-1)
=3^{q-1}-2^{q-1},
\]

and, when the low root `r` is fixed at original position two,

\[
\sum_{\ell=3}^{q+1}\binom{q-2}{\ell-3}2^{\ell-3}=3^{q-2}.
\]

The second count is the sharper one for top-two words and is the computation
audited here. It is an upper bound:
the general code space allows labels that are not valid branch assignments,
which does not harm injectivity.

## Explicit code and inverse

Fix the deterministic parent of each pre-merge nonroot as its earliest exposed
neighbour. At the first merge `w` (one-based rank `q+1`), select in each
branch the earliest exposed neighbour of `w`. Their parent chains form the
simple path `p -> ... -> w -> ... -> r`, of `ell` vertices. Move that path,
in this orientation, to the front and keep every nonpath vertex in original
order.

For a fixed one-peak output, retain:

1. `ell` and the `ell-2` original positions of path vertices other than `p,w`
   among positions `2,...,q` (with `w` fixed at `q+1`);
2. the position `j` of `w` along the output path; and
3. which of those retained positions belong to the `p` arm.

The inverse fills the selected `p`-arm positions in forward path order and
the remaining `r`-arm positions in reverse path order, places `p,w` at their
fixed ranks, and fills all holes with the preserved output tail. Hence the
same output and code cannot come from two input words. With `r` known to be
at original position two, it is forced to the `r` arm, leaving exactly the
`3^(q-2)` code capacity above.

The path begins at the global-high root `p`, so the moved prefix is connected.
Every retained tail vertex still has an old predecessor or a newly earlier
path vertex; the output is one-peak. Since all path vertices lie among the
first `q+1` original vertices, and the tail supplies the remaining members
of that prefix first, the output has the same first-`q+1` vertex set.

## Actual finite audit

`audit_path_prefix.py` used all `9! = 362880` permutations of the `3x3` grid
and encoded the 81,984 words in the top-two `B` slice. For every encoded word
it checked exact decode equality,
one-peak validity, preservation of the fixed prefix vertex set, and absence
of a duplicate `(output, code)` pair.

The final low-priority, single-process run used CPython 3.15.0 with niceness
10 and exited zero. Its detailed version, UTC time, output and source hashes
are in `audit.exit.json`; full output is in
`path-prefix-independent-r033.json`. In particular, observed preimage maxima
for each fixed output were bounded by `3^(q-2)` for every realized `q`.

Scope correction: this run is **not** an exhaustive encode/decode audit of
the full `K=2` class. The complete `3x3` `K=2` total is 143,808 words; its
separate full computation belongs to `primary merge_geometry`. Therefore the
`3^(q-1)-2^(q-1)` A2 expression is audited here only as an analytic code-space
identity, not as an independently enumerated all-A2 injection. No conclusion
about a full-A2 finite enumeration should be drawn from the 81,984 count.

This audit supports the first-merge rank counting lemma and its base-3
exponent. It does not by itself prove a peak-distance statement, nor does it
replace the separate probabilistic step needed to turn the count into the
claimed high-probability rank tail.

Reproduce:

```sh
cp audit_path_prefix.py "$OUT/"
(cd "$OUT" && nice -n 10 python3 -B audit_path_prefix.py)
```
