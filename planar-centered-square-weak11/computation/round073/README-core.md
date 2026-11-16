# R73 core: exact-cap atomic repair

The mathematical source is `../../proofs/r073-cap-preserving-atomization-and-local-full-halo-obstruction.md`. It proves all-square caps for every n>=2 and0<delta<=1/(2n), exact preservation of maximum mass1, and weak approximation of the R72 diffuse family. Outside contraction makes the selected block unique. A formula-defined79056-source input is connected/nonterminal and has local full-halo area greater than7/2. This is not an actual deletion-area counterexample: the counted outside singleton halos are residual halo and must be subtracted.

`verify_cap_preserving_atomization.py` uses only Python standard-library Fraction. The parameter lists and explicit n9 fixture are in the source, and the result records every fixture parameter. Reproduce from this topic folder with:

    python3 computation/round073/verify_cap_preserving_atomization.py /Volumes/data/CodexTemp/T2/r073-core-replay.json

The original run used Python3.15.0, PID26306, exit0, confirmed absent. It checked3092 k,l pairs,9276 exact radius identity samples and42 zero-capture cases. The analytic proof, not these samples, handles all parameters and radii. Exact index counts give78975 outside atoms and28179 singleton squares in C outside the collar. No79056-source enumeration, full subset list or complete halo area computation was performed. The explicit lower bound is7832077922743/2230740000000; it is a lower bound on full halo intersect C, not on full minus residual.

`cap-repair-result.json` and `cap-repair-receipt.json` preserve actual output, software version, source hash, process status and executed command. Other components independently check the formulas and region bookkeeping; their runs are not counted as core runs. Root owns Wolfram and the separate29-source complete-area control. No MATLAB process was started by core in R73. No commit or push is performed in this workspace.
