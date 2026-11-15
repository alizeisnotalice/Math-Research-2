# R64 one fixed-position mass-polytope pilot

`mass_polytope_variation.py` is a new module and does not change/import frozen
R63 production code. It consumes the core's exact 21-source positive seed,
with input SHA `46acda723e8590054912a18f71c67efc95b5a0bcfd449a9c60ce615413ac41dd`.
Root owns execution, resource monitoring and final integration.

The seed is symmetric and lies on endpoint-order walls. Central finite
differences along mass-S-preserving coordinate directions propose an LP
objective; neither those numbers nor SciPy's LP status are certified ordinary
or chamber derivatives. The LP maximizes this numerical proposal over all
complete bbox mass caps, the unique-maximum gap, positive floor, nonterminal
total-mass cap and prescribed connection inequalities, with a local trust box.
The finite line may cross walls. No fixed-chamber theorem is required or claimed.

The LP vertex is rounded to dyadic rational masses; one S mass is adjusted so
mass(S)=1 exactly. Each of eight default dyadic line steps is then independently
checked against **all** exact constraints. Infeasible steps are recorded and
cannot produce a success claim. Every feasible candidate is evaluated with
complete source bbox closures, strict inner/outer rational radius enclosures
and exact integer rectangle union differences for the full and residual halos.
Success requires

    candidate deletion lower bound > seed deletion upper bound.

The saved JSON retains all candidate masses, exact checks, area enclosures and
gain lower bounds, as well as software/source/input identities. A certified
improvement is for this fixed input family only, not a general deletion theorem
or a new C2 bound. A failed line does not prove local optimality.

Typical root command with an existing NumPy/SciPy interpreter:

```sh
python3 -B /absolute/path/mass_polytope_variation.py \
  --input /absolute/path/r064-fixed-position-positive21.json \
  --out /Volumes/data/CodexTemp/T2/r064/lp-fixed21-pilot.json
```

Defaults: one LP, trust radius 1/16, finite-difference step 1/262144,
dyadic mass rounding denominator 2^28, 8 line steps, exact endpoints at 72 bits.
The seed's exact area enclosure is checkpointed before numerical proposals.
Temporary numerical proposals need NumPy/SciPy; all geometry certificates use
only Python integers and Fraction. Source coordinates must be dyadic at the
chosen endpoint precision.

`continue_certified_line.py` consumes the actual pilot result and reuses its
strictly certified winner on the same positions. It pins the original module
hash, changes only the next input's initial masses, and preserves all original
constraints. It runs at most eight child LP calls, accepts a winner only if its
new lower area bound exceeds the immediately preceding upper bound, and
records every input/output/parent hash. Exact source-hypergraph connectivity,
positive support and actual nonterminal margin are additionally checked for
each accepted state. One no-gain call halves the trust radius, two consecutive
no-gain calls stop, and success restores radius 1/16. Exhausting this finite
budget proves neither local optimality nor an upper bound for the class.

Accepted-state records contain every strictly eligible bbox closure's mask,
rank and mass. The `competitor_overlap_max_mass` field is the largest such mass
for a closure intersecting S and different from S. It is a conservative
geometric score; the code does not claim that each closure is realized as a
positive-area maximal actual-source state.

Actual verification recorded on 2026-11-15: root ran the scaled 40-bit seed,
the 72-bit pilot (one LP/eight finite line points), and the continuation (eight
LP calls/64 finite line points). All 72 line candidates passed exact feasibility;
42 had strict improvement over their own current seed, and nine winners were
accepted successively. The final exact deletion lower bound is
27563576327882167757252602439668073627139913/
22300745198530623141535718272648361505980416, about 1.2359935097459571.
Continuation stopped at `FINITE_CALL_BUDGET_EXHAUSTED`, with eight accepted
continuation winners. This is neither convergence nor local optimality.

The saved-ledger summarizer verifies source/parent/output hashes and successive
strict interval comparisons; it does not rerun geometry. Independent geometry
review is supplied separately by the audit branch. The pilot has an actual
startup/exit/PID receipt. The continuation's retrospective observation receipt
records actual tool-session exit and all eight child exit codes; individual
PIDs were not saved, so no per-PID disappearance is asserted for that phase.
