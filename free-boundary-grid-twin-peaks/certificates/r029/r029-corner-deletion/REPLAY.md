# Portable replay

Run from this directory.  The first command audits the n=3/n=4 local map and
integer DP; the second is the independent n=5 inclusion-exclusion check.

```sh
python3 -B corner_deletion_audit.py
python3 -B audit_h5_independent_ie.py \
  --tail-summary ../r029-tail-full-n5/n5-tail-summary.json \
  > verify-n5-independent-ie.json
```

The n=5 output asserts `T_22 = 8 A1(H)` exactly.  The general corner-deletion
lemma is analytic; neither finite audit closes the twin-peaks problem.
