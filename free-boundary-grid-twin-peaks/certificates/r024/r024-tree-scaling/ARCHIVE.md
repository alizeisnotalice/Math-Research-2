# r024 frozen exact tree scaling computation

The sole production source in this package is
`tree_peak_dp_r024_frozen.py`, SHA-256 `76eddd4a…`. `rerun_frozen.py` was
actually run against that exact source and produced `frozen-rerun.json`
(SHA-256 `f967ef59…`), with resource receipt `frozen-rerun.resource.txt`.
It exactly recomputes complete-binary trees through height 9 and two equal
complete-binary trees joined root-to-root through height 7.

For the joined family the finite ratios rise from 1.4211 at height 1 to
17.5825 at height 7. These are finite exact confirmations.  The separate,
reviewed all-height theorem in `../r024-global-normalizer/double-binary-b-upper.md`
now proves that no such universal \(C_3\) exists; this computation is not its
proof.

`CUT_FORMULA.md` contains the corrected ancestor factor \(G_e(r)\) in the
hook-length identity. The older formula without this factor is excluded.
`verify_cut_formula.py` (SHA-256 `88b6d164…`) independently checks the
corrected formula; its actual output and resource receipt are retained.
`exhaustive_cut_formula.py` (SHA-256 `537bf42b…`) is the original Prüfer
scan through six labelled vertices (3+16+125+1296=1440 trees), whose retained
output agrees with the subset DP.  Its historical absolute module path is
preserved in the source; `verify_cut_formula_archived.py` is a path-only
archive adapter for a portable rerun.
Neither tree calculation is a free-boundary grid counterexample.
