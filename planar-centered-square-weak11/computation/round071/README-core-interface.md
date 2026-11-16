# R71 core exact-interface checks

Run from the topic directory with an existing Python 3 interpreter:

    python3 -B computation/round071/check_shared_mass_interface.py --input computation/round071/shared-mass-interface-collar-input.json --output /Volumes/data/CodexTemp/T2/r071-replay/collar.json
    python3 -B computation/round071/check_shared_mass_interface.py --input computation/round071/shared-mass-interface-input.json --output /Volumes/data/CodexTemp/T2/r071-replay/interface.json
    python3 -B computation/round071/check_receiver_cell_minimax.py --input computation/round071/receiver-cell-minimax-input.json --output /Volumes/data/CodexTemp/T2/r071-replay/cells.json

Create the destination directory first. No external mathematics package
is needed. All arithmetic is Fraction/integer; no random choices occur.
The inputs contain complete positions, mass menus, receiver cells and
radii. These checks do not optimize an LP or compute complete halo areas.

Actual final runs used Python 3.15.0. The two seven-position fixtures each
evaluated all 64 mass vectors and all 127 source-subset predicates per
vector (8128 actual row evaluations), 37 distinct closures and 512 residual
receiver comparisons. The collar fixture has only one admissible vector;
the exterior fixture has nine. Their full-radius-menu and deletion-point
checks number 8 and 72, respectively. Four strict/tie/zero-distance edge
fixtures and a matrix-only weak-duality diagnostic also passed; the latter
is not a geometric counterexample. Final PIDs 9710 and 9713 exited zero
and were absent after wait. Wall times were 0.110432 and 0.121929 seconds.

The receiver-cell kernel checked 279 exact subset/cell minimizers, 279
maximizers and 1953 open-receiver/closed-cell intersection comparisons.
Degenerate cells are included. PID11745 exited zero and was confirmed
absent; wall time was 0.051877 seconds.

Two earlier successful logical checks used a short-circuit generator,
so their field named all_source_subsets described a potential range,
not actual complete row execution. They are not used for the final actual
count. The exact correction and initial receipt identities are recorded
in interface-scope-correction.json. Both final interface fixtures were
rerun after replacing the generator by explicit complete row evaluation;
final outputs and receipts bind the corrected script.

The accompanying proof treats fixed positions, all residual radii,
strict shared witness margins, complete tied-distance menus, arbitrary
outside measures through pointwise geometric cuts, zero-cost closed
collar coverage, and exact whole-receiver-cell residual exclusion.
Independent proof review is a separate component. No new C2 bound or
general deletion-area upper bound follows from these finite checks.
