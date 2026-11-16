# r024 tree root-message audit

`verify_tree_messages.py` is the historical audit driver: it enumerated all
rank permutations on random small-tree instances and checked the tree
hook-length \(A_1\) calculation and root-message invariance. `receipt.json`
records that completed run against source SHA-256 `6615204c…`. The matching
source bytes are not present in this staging package (the live precursor has
since changed), so this is retained as a historical receipt, **not** a
portable rerun claim. Any later frozen/re-run tree-integral source is
represented separately in `../r024-tree-scaling/` with its own exact source
and receipt.
