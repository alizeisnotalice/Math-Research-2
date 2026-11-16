# Reproduction

```sh
python3 domino_classification.py > domino-classification.json
python3 orientation_audit_n3.py > orientation-audit-n3.json
```

The program uses only exact finite set operations.  It audits every domino in
the complete free-boundary grids \(P_n\square P_n\), \(3\le n\le7\), and has
no random input.

The second program enumerates every \(3\times3\) word and audits the precise
two-to-one collapse caused by renaming the heavy domino side as the helper.
