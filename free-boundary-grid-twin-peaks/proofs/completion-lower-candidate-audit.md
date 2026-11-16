# Two-step completion lower bound: audit r010

Let a finite connected graph have a revealed connected set `S`, with
`m=|V\S|>=2` and initial frontier `F=N(S)\S`, `D=|F|`.  A completion word
lists the remaining vertices, each adjacent to the previously revealed set.
Under the uniform law on these words, put

\[
p=\Pr\{\text{the second vertex is not in }F\}.
\]

The candidate statement examined here is

\[
p\ \ge\ \frac{m-D}{D(m-1)}. \tag{L}
\]

It is a useful auxiliary statement only.  Even a proof of (L) would not
by itself give a positive lower bound under the original two-peak
conditioning.

## Proven special case: rooted forests

If every non-frontier vertex belongs to one of `D` disjoint rooted trees,
whose root is a member of `F`, let their sizes be `L_1,...,L_D`.  Interleaving
the independently legal tree orders gives

\[
p=\frac{\sum_iL_i(L_i-1)}{m(m-1)}
 \ge \frac{m-D}{D(m-1)},
\]

where the inequality is Cauchy--Schwarz.  This includes arbitrary internal
tree shapes.  It does not extend by a naive edge-monotonicity argument.

## Strict failure of the naive forest reduction

For root `0`, take the five-vertex graph with adjacency bitmasks

```
[12, 28, 3, 3, 2]
```

so the root frontier is `{2,3}`.  Direct exact completion counts give
`p=2/3`.  Adding only the edge `{2,4}` inside the remaining graph (2 is initially available and 4 is not) gives `p=3/5`.
Thus adding an internal edge can decrease `p`; deletion to an arbitrary
rooted spanning forest cannot prove (L).  Both values still satisfy (L).

## Actual exact searches

* The historical `baseline-search.py` (temporary audit directory) is the inherited full labeled-graph enumeration:
  it tests all connected graphs through six vertices and all connected
  prefixes (944824 prefixes at six vertices).  It finds no counterexample.
* The historical `random_search.py` (temporary audit directory) uses integer completion recursion on 22000 deterministic
  random rooted graphs: 2500 each for 7 through 14 vertices and 500 each for
  15 through 18.  It finds no counterexample.  Its equality example at
  `n=7,D=2` has a triangle in one root component, so equality is not a
  tree-only phenomenon.
* `edge_monotonicity.py` is an exact search whose planned upper size was six vertices;
  it stopped at the five-vertex counterexample above and did not complete six vertices.

The scripts use `Fraction` and integer recurrence values; these are finite
checks, not a proof of (L).

## First unproved step

What remains is a direct injection or a valid extremal theorem saying that
the uniform rooted-graph completion law is bounded below by the balanced
root-component expression.  The apparent assertion that adding an internal
edge can only increase `p` is false.  Existing generic antimatroid references
on balanced pairs do not supply (L), so it must not be cited as an available
theorem.

Historical exploratory search outputs remain under /Volumes/data/CodexTemp/T2/a022-completion-lower-proof. Only the strict minimal monotonicity counterexample and its reproducer are required in this round archive; random search is not invoked as a proof.
