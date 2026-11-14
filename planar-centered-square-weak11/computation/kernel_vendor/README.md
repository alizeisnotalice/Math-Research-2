# Unchanged kernel evaluator

`hl_core.py` and `LICENSE` are byte-identical Apache-2.0 files from the old
project's `kernel-rescaling/vendor/`, originally from
https://github.com/tejasvi/square-maximal-weak-type-v4 at commit
`3af09fc359ce262f4a7f63d34b0d97e6d4d95d1a`.

SHA256 of `hl_core.py`:
`c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf`.
SHA256 of `LICENSE`:
`cfc7749b96f63bd31c3c42b5c471bf756814053e847c10f3eb003417bc523d30`.

The new family-floor program uses only its closed-form kernel evaluator,
not the historical C1-C4 checker. A replay from `floor-neumann.json`
does not require the old candidate coefficients, old T directory, SciPy,
the floating LP or its full point menu. It requires python-flint0.9.0 and
the decoded archived preconditioner; see `proofs/kernel-family-floor.md`.
