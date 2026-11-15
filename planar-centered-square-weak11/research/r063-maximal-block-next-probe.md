# R63 directed probe and continuation

2026-11-15. Core proof: `proofs/r063-maximal-compact-block-deletion.md`. The actual general upper bound has not changed. The selected next mathematical target is its conjecture (6), not an assertion that the finite probes will settle all inputs.

## Decision rule

For every complete source input, first determine the eligible hypergraph components. A component with d²≤(35/33)m is already paid by R53 and is not an obstruction. On each remaining connected component A, enumerate all bbox closures and let L be the largest mass among strict eligible closures. Retain **every** maximizer S, not merely those with a bounded source count.

For each S, compute certified full residual halo differences as in (7). A single upper deletion ratio≤7/2 supplies a valid next deletion for that input. A claimed counterexample to the maximum-block conjecture requires a lower ratio>7/2 for every maximizer of one genuinely nonterminal connected input. Low-rank arbitrary compact blocks are an optional stronger menu and have separate minima; their failure is not a counterexample to the selected conjecture.

The elementary rigorous score (1+w/sqrt(L))(1+h/sqrt(L)) already pays S whenever one bbox span is ≤3sqrt(L)/4. Prioritize examples for which all maximum blocks are wide in **both** coordinates. Use exact checks 16w²>9L and 16h²>9L. A 3×3 equal-mass grid with suitably large spacing gives this hard shape but is only a control, not an adversarial maximizer of true birth.

## Where the useful geometry lies

Every S has an exact source-free collar and every outside T obeys m(T)≤d(S∪T)²−L. All actual witnesses have radius<sqrt(L)/2. The true deletion problem can therefore be normalized to L=1 and restricted to bbox(S)+[-1,1]² without changing the answer; this local input has mass≤9. These constraints permit a finite local mass relaxation or geometric deficit argument, while preserving true source overlap and the old residual halo.

The precise missing inequality is

    |I_S\H(A)| + |I_S∩H(A\S)| ≥ |I_S|−7L/2

for at least one maximum S whenever its right side is positive. A proxy area or a raw shell mass does not furnish either term. A proposed critical exterior background must obey the large-square caps, not just the empty collar; otherwise it can introduce eligible blocks of mass>L and invalidate the normalization.

## Bounded computation ownership

The selected Sol `next_compute` owns the new frontier implementation in its isolated `pcw-position-slab-impl` workspace. Root owns actual frontier executions and receipts. The initial menu is twelve structural strata with six 16–20-source inputs per stratum, plus earlier small-case equivalence controls. This design and the existence of source files are not a claim that all 72 cases have run. Actual completed counts must be copied from root's receipts and the saved results.

Strata should emphasize wide connected inputs, tight two-dimensional clusters with bridges, thin frames, unequal masses, and boundary additions that retain old source reserves. Do not spend an unbounded sequence of rounds merely enlarging low-cost random finite inputs. If the finite deletion scores remain far below 3.5, the next work is the geometric deficit inequality above or a rigorous local mass relaxation using the collar and square caps.

## Completed core sanity run

The core itself ran only the small independent check saved in `computation/results/round063/maximal-block-structure-exact.json`: Python 3.15.0, 32-bit strict root brackets, exact integer/Fraction decisions, nine inputs, 1,981 nonempty subsets, 25 maximal blocks, 211 square-cap checks, 96 collar point checks, 500 exterior-growth checks, and eight strict spatial restrictions. All checks passed in 0.1116 seconds, exit 0. Localized deletion intervals overlap their independently recomputed full-input intervals. This tests finite identities, not the uniform selection conjecture.

No new mathematical general bound, exact C2 value, new legal extremizing sequence, commit, push, or backup verification is claimed here. Parent integration remains responsible for the full R63 report and independent review.
