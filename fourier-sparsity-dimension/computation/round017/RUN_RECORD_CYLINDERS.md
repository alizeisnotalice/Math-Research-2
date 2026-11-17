# Round 017 independent cylinder audit

Source certificate: `nonparallel-cases.json` from the main worktree.  The
independent script `check_cylinders.py` imports no project code; it uses exact
GF(2) bit arithmetic and integer Walsh sums.  Run:

```sh
python3 check_cylinders.py \
  --input /Volumes/data/CodexProjects/T2-A28/fourier-sparsity-dimension/computation/round017/nonparallel-cases.json \
  --output cylinder-audit.json
```

The script also has `--input` and `--output` parameters whose defaults are
`nonparallel-cases.json` and `cylinder-audit.json` in its own directory.  The
recorded corrected run exited with status 0: script elapsed time was
12.059140708 s (wall `time -p`: 12.08 s), on Python 3.15.0.

It asserts input domain size 64, basis independence, and that blocks are a
disjoint cover before it reconstructs each affine block direction, restricts every phase in the
provided code, obtains local ANFs, freezes the combined non-star quadratic
map, and then constructs the two half-cylinders (7) for every double block
and (8) for every minimal block.  It exhaustively checks all pairs in
`V0 × F2^6`, so it checks cylinder affine-flatness, disjointness, covering,
and equality of each cylinder slice with the direct integer partial-Walsh
support.  It additionally asserts `r<=3b`, `dim(V0)=dim(V)-r`, normal rank
at most 7, `dim(K)+rank(projection)=dim(V0)`, every one/two-dimensional
block restriction bound, and the joint restriction bound.

Results in `cylinder-audit.json`:

| case | b | r | dim V0 | normal rank | projected rank | dim K | restriction-image rank |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| four nonparallel double blocks | 4 | 2 | 9 | 4 | 4 | 5 | 5 |
| three double + two minimal | 0 | 0 | 10 | 5 | 5 | 5 | 5 |

All eight cylinders in each case have codimension three in the joint space,
are pairwise disjoint, and cover it.  Every tested Walsh support agrees
exactly.  The per-block restriction ranks of K are respectively `(2,2,2,2)`
and `(1,1,1,1,1)`, both within the asserted one/two-dimensional local bounds;
the joint image rank is 5, below the bound `m1+2*m2=8`.

This is a finite exact audit of the two supplied codes, not a proof of the
general rank-budget theorem.
