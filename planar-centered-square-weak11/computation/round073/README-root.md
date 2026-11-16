# R73 root computations

All computations were actually executed on 2026-11-16. They support the
round's structural investigation; they do not establish an exact C2 or a
general 3.5 deletion bound. No MATLAB process was started in this round.

## Small D4 diagnostic

`probe_d4_atomic.py` used the archived R72 complete-closure implementation
read only. It considered n=3,5,9, four equal outside atoms placed on the
four axes or four diagonal corners, and eight rational displacements.
All 48 position configurations were checked with exact rational maximum
mass caps. Four caps were zero; the other 44 each used five mass ratios,
for 220 floating complete-area evaluations. Area and ranking results are
numerical evidence only. No general or family-wide optimum is claimed.

Python 3.13.18 used three worker processes, 25474/25475/25476. Parent
25472 exited zero and all four PIDs were absent after the pool joined.
Wall time was 3.751 seconds. The finite probe found no near violation of
the corner-only candidate; this is not evidence proving that candidate.

## A simpler exact finite control

The selected n=5 axial configuration has the 25 old grid atoms, plus
(-1/2,12/25), (73/50,12/25), (12/25,-1/2), (12/25,73/50), each with mass
(8/25)*(99999999/100000000). Complete 72-bit directed area calculation
with the archived R72 code verified 476 closures, 171 eligible closures,
a unique maximum S25 of mass1, one eligible component, and nonterminal
gap 733950007/515625000. Its exact deletion lower endpoint is greater
than 79/50; the decimal display is about 1.5826878638260273.

This simplifies the already established finite deletion-budget obstruction.
It is not a new C2 lower bound. PID26372 exited zero and was absent after
wait; wall time .111 seconds. An independent algorithm is recorded by
the R73 independent component, not by this root production receipt.

## Wolfram scalar checks

`r073_cap_repair.wls` checked nine scalar identities or real-domain signs
for the enlarged-hole atomization. All returned True in Wolfram15.0.0.
In particular the unit-square excess factors as
-(1-2*n*delta)*(n-k)*(n-l)/(n^3*(1+2*delta)); the general-radius remainder
has the required sign for R>=1 and 0<=delta<=1/(2*n).
The full geometric capture and cell-count arguments remain analytic proof
obligations. The scalar calculation does not enumerate all source subsets.
Launcher25839 and kernel25841 exited zero and were absent; wall1.514s.

## Portable reproduction

From the topic directory, send outputs to an existing directory under
`/Volumes/data/CodexTemp/T2`:

    python3 -B computation/round073/probe_d4_atomic.py --geometry computation/round072/actual_halo_geometry.py --output /Volumes/data/CodexTemp/T2/r073-replay/d4.json --workers 3
    python3 -B computation/round072/actual_halo_geometry.py --input computation/inputs/r073/s25-four-axis-unique-input.json --output /Volumes/data/CodexTemp/T2/r073-replay/unique29.json --bits 72
    wolframscript -file computation/round073/r073_cap_repair.wls /Volumes/data/CodexTemp/T2/r073-replay/cap-scalar.json

The exact input, final code, outputs, versions, parameters and process
receipts are in the root component manifest. Archived process IDs refer
to those actual executions, not to any future reproduction.
