# R72 core proofs and exact scalar checks

The two core proofs are:

* uniform-grid Borel square caps, sharp constant density, critical-radius
  halo criterion, and a canonical atomic-transfer obstruction;
* outside-mass contraction, unique maximum-block selection and continuity
  of actual deletion area, giving the stated universal non-strict budget
  equivalence on perturbation-stable domains.

Run the scalar check from the topic directory after creating the output
folder:

    python3 -B computation/round072/verify_borel_family_scalar_bounds.py --input computation/round072/borel-scalar-input.json --output /Volumes/data/CodexTemp/T2/r072-replay/borel-scalars.json

Only standard-library Fraction arithmetic is used. There is no random
seed or tolerance. Actual Python was 3.15.0. The run checked 4552 cap-pair
values, 4460 exact differences, 1119 interval-overlap event/gap samples,
five sharp-density and limiting-upper-bound formulas, and six exact
canonical-grid cell counts. These are finite checks of analytically
proved formulas, not a replacement for their all-n proofs. The cell
counts certify residual singleton-halo area lower bounds only; neither
complete atomized halos nor actual maximum eligible masses were computed.

The process PID20772 exited zero and was confirmed absent. Wall time was
0.102582 seconds, with 0.071031 seconds inside the checker. Input, output,
source hashes and the actual command are retained in the receipt.

The scalar output also records the exact integral w/2 of the left-unit
strip cumulative capacity 1-k_n(a)/n, where k_n(a) counts internal columns
at coordinate at most 1-a. It is a necessary weighted outside-mass budget,
not a proved deletion-area inequality and not the main result of this
round. No new 3 or 3.5 area bound is inferred from it.

Independent internal-halo shape checks, actual finite-source area/cell
optimization and uniqueization, Borel adaptive-area certificates, and
root MATLAB/Wolfram work are separate components. The core did not rerun
those computations. In particular the new finite contraction example
uses the R72 actual input rather than relabeling an old R71 calculation.
