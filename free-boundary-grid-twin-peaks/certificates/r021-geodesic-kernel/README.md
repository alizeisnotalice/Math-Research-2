# Geodesic-prefix average audit

This directory studies `p_d <= C^d/(d!)^2` for the uniform one-peak law. It remains open: no uniform constant is proved, and no construction here refutes every fixed constant.

`local-extension-audit.md` concerns a stronger pointwise estimate, `q(gamma) <= C/(j+1)^2`. On actual 5x5 data one local witness requires C about 20.67. Therefore the particular attempt with C=5 is rigorously false. This does not refute the global average bound or prove the required local C diverges with n.

`frontload-completion-comparison.md` gives an all-n completion lower bound. Its direction cannot establish the candidate, and it is kept to prevent that logical error.

The reference prefix numerators for n=4 and n=5 are in `reference/`. Reproduce the finite derived data without project-specific absolute paths:

```sh
python3 analyze_geodesic_candidate.py
python3 analyze_extension_kernel.py n3-kernel.json n3-local.tsv
python3 analyze_extension_kernel.py n5-kernel.json n5-local.tsv
python3 verify_frontload_lower_n3.py
```

The last command checks all 218 connected nonempty 3x3 subsets. The compiled extension executable and caches are omitted.
