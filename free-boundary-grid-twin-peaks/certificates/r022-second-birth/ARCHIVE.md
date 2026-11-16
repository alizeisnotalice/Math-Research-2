# r022 second-birth survival certificate

Source snapshot: `a022-second-birth-survival-r022`.  This package keeps the
three scope-limited proofs (second-birth survival, ideal-exact-arithmetic
tau=2 SMC first moment, and first-merge Rao--Blackwell identity), their n=3
exact JSON outputs, and their Python verifiers.  The finite computation is an
exact 3-by-3 audit; it is not a tail bound or a finite-SMC ratio theorem.
Primary verifier SHA-256: `4779df6cee5b6470e9b542c2ca9a01095c0350f4921e3d817f4900c2a24f0bd3`.
Relative-path archival adaptation SHA-256: `853d99f61c2d753adfc554336212169b9839f62cc8aec6dbf4d9c04d77539621`.

`independent-root-audit/` is a separately written direct 9! position-mask
enumerator.  Its actual result has A2=143808, D=633600, top2=81984 and
sum(L^2)=358080, and verifies every distance-by-birth cell against the
L-lift table.  Run it with `python3 -B independent-root-audit/audit.py`.

Run locally from this directory: `python3 -B verify_second_birth_survival_n3_relative.py`,
`python3 -B verify_first_merge_rb_n3.py`, and
`python3 -B verify_tau2_prefix_fk_n3.py`.  Historical absolute paths retained
inside the untouched historical verifier describe its original run only; the
`*_relative.py` adaptation uses the two retained minimal inputs and is the
archive reproduction entry point.
