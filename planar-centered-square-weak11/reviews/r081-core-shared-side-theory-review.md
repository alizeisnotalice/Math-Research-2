# R81 core adversarial review: where the shared-side budget is actually used

2026-11-16. This review concerns the fixed twelve-position alternating family in [the shared-side theorem](../proofs/r081-dual-shared-side-corner-payment.md). It is mathematical and source reading, with no additional numerical run. The separate author, producer and independent-audit components retain their own actual receipts.

## 1. Qualification is a continuum statement; area at two endpoints is not

The production source `computation/round081-compute/shared_side_family.py` constructs all closed coordinate-interval products from the twelve fixed positions. The441 possible products yield157 distinct nonempty full bbox captures. A captured set's tight bbox gives the same capture mask: shrinking to the extreme captured coordinates cannot add a source outside the original box and cannot remove a captured source. Thus the menu includes every possible square capture and every bbox closure.

For each capture, mass is affine in t while bbox and capacity are constant. The exact half-line intersection in the saved output is[0,1/2]. Positive masses restrict the physical family to0<t<1/2; the requested[1/4,499/1000] is wholly inside. Intersecting these affine constraints, rather than checking a parameter grid, certifies all scales and positions for every requested t. Since S has mass1 and span squared<1, the complete caps establish S as a maximum block.

The saved area engines run only at t=1/4 and499/1000. Both endpoints have precisely the five maximum masks15,99,390,1548,2073. The report should retain this observed endpoint scope rather than silently treat endpoint enumeration as a complete classification of all intermediate ties. The continuous area theorem needs S to be a maximum, not a unique maximum. Neither the calculation nor the proof contracts the outside masses, and no monotonicity under R72's uniqueness perturbation is invoked.

The source separately builds the outside-only closure menu, which is essential: a residual closure must not be formed against the old atoms. Its104 nonempty captures give the complete residual halo. Algebraic cancellation of pure-outside full labels precedes directed area bounds. Lower area uses full inner minus residual outer, with the reverse orientation for upper area. Internal old-halo endpoints are exact in this fixed four-corner geometry; the retained-old and net-mixed directed pieces really add to deletion at both endpoints.

The second production source `regional_credit_postprocess.py` performs new area operations on saved endpoint rectangles. It recomputes cropped full and residual union areas and the uncovered-full/residual deficit intervals; it is neither a new parameter input nor a mere hash/postprocessing summary. Its mathematical operations and its separate process receipt must be reported.

## 2. Local geometry: the special mass assignments matter

Every effective corner window contains exactly one old atom of mass m=1/4 and its own two outside atoms. Full caps imply all strict radii<1/2, so restriction to that window is exact. The common outside measure remains present in the global caps and complete residual. At each corner let a,b denote near/far mass and let

    v=1/4, u=√(m+a)/2, z=√a/2,
    U=√(m+a+b)/2, V=√(a+b)/2.

On the requested family a lies in[.001,.25], so z>e, z<=v<=u and u<d-e. These inequalities validate the author's pure-old thin strip and near-singleton overlap. The overlap is exactly(u+z-d)_+(z-e), and the near union formula does not discard the old strip. The resulting G(a) formula is exact for this geometry.

Even corners have a=m: the far mixed receiver and pair residual have the same radius and the same crop whenever nonempty. Odd corners have b=m: the far mixed receiver is empty since√(1/8)<D0/2. Both facts use the particular alternating assignment. For arbitrary a,b, the far mixed label cannot simply be omitted.

The pair residual lies above the whole near part, while the far residual singleton misses the crop. The joint receiver also lies above the near part because D0-U>u throughout this family. Thus G(a) and joint contributions may be paid separately here without overlapping residual bookkeeping.

The author's bound G(m)<.063 and the two cases proving G(a)<.1 are valid. In the active-overlap case, z²−zu<=0 and

    .605u−.356z >=sqrt(.605²−.356²)/4 >.12

follow from u²−z²=1/16 and the positive reverse-Cauchy square. This is a continuum inequality, not an interpolated numerical observation. The bound for four near contributions is2(.063)+2(.1), and uses exactly two corners with a=m.

## 3. The cross-corner step is genuine but restricted

The actual unit window on side i contains two old endpoints of total mass1/2, near_i, and far_(i−1). Thus full global caps force

    a_i+b_(i−1)<=1/2.

The family saturates these four constraints. Under its prescribed alternating assignment, this implies the actual local totals

    B_even=1/2+eta, B_odd=1/2−eta,
    0<=eta<=249/1000.

These are common source masses; four independently optimized local measures could not be substituted here. The active joint threshold is B0=D0²−m, and the derivative bound on the full joint receiver gives

    actual joint deletion<=kappa(B−B0)_+,
    kappa=1143/4000.

Using full joint area only on the right is safe, since the unchanged complete pair residual can reduce its contribution. Equality at B=B0 gives an empty strict receiver and is handled correctly. The function is continuous there. The derivative bound uses B<=.75, hence√(m+B)<=1.

The essential shared-budget calculation is

    sum_i(B_i−B0)_+
      =2[(delta+eta)+(delta−eta)_+]<=10719/20000,
    delta=1/2−B0.

One cannot maximize all four hinge terms independently. This hinge estimate, together with the two strong/two weak near estimates, is what proves the strict aggregate<.481. Adding the entire noncorner area gives complete deletion<3.480997 for every requested t, once the exact cap qualification is imported.

It would be an overstatement to say that arbitrary masses satisfying only sum B_i<=2 have now been paid. Both the paired relation and the special identities a_even=m,b_odd=m, as well as d,D0,w, enter the proof. The result is more than two endpoint computations, but remains a single fixed-position, one-parameter theorem. It does not cover general source positions, arbitrary old CDF profiles or unrestricted outside rank.

## 4. Limits and final scope

The uncapped affine-cost counterexample in the author note has mass1.01 in a unit window, so its old singleton fails the mass1 maximum-block qualification. Its algebraic inequality is correct, but it does not disprove a qualified affine cost statement. The successful family proof does not use that false uncapped shortcut.

The two upper-endpoint strong corners reproduce R80's actual local half-mass failure; the other two corners provide shared reserve. No pointwise D4 symmetrization or independent pasting is asserted. No general C2 bound or arbitrary maximum-block selection theorem follows. The useful transferable structure is the need for a common side-capacity state together with geometry-dependent local cost bounds; proving such costs for general positions and old mass profiles remains the missing step.
