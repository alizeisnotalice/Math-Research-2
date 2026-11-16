# R76 general-grid joint-null payment

The proof `../../proofs/r076-general-grid-joint-null-corner-payment.md` proves all n>=3 under an actual joint-null condition. It retains original internal positions in the joint-witness test and full-grid capacities; single-outside auxiliary radii are explicitly truncated at1/2 before the clipped annulus estimate. It does not assert that the collapsed measure satisfies the old caps or that a collapsed joint label is absent.

From the topic directory, reproduce the small exact scalar regression with:

```sh
python3 -B computation/round076/check_general_joint_null_constants.py computation/round076/general-joint-null-input.json /Volumes/data/CodexTemp/T2/r076-general-joint-null-replay.json
```

Python3.15.0/Fraction actually checked127 n values,8255 pure-prefix inequalities and635 truncated-radius squared identities. PID50576 exited0 and was confirmed absent. The saved result records all rows and input/source pins; the receipt records the actual command and timing. This is not a source-position search, a finite-input capacity certificate, a joint-prefix activity test or a halo-area integration. General quantifiers are proved analytically; the independent reviewer reads this scalar run without repeating it.

The independent n5 original-prefix fixture and n3 finite probe are separate components. They must not be counted as runs of this code. Frozen old-round files remain unchanged. No MATLAB, agent creation, commit or push was performed in this clone.
