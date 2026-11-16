# R79 independent review: continuous states, saved quadtrees and general-block direction

2026-11-16. Final mathematical review PASS for the sources below. This review separates the single new independent execution from analytic/source review; it does not certify an unrestricted C2 bound.

## Mathematical checks

The quadratic envelope proof retains one actual common radius in (R_J,1/2). Its origin-bearing reflected bbox gives Rlo >= max(XL,YL) and Rhi >= max(XU,YU), so coordinate comparisons of F are used only on their legal radius domain. The minimum over a closed radius interval is equivalent to the strict open-radius test only when R_J<1/2 and the mass inequality is strict. Equality, B_J=0 and closed crop cells are handled correctly; cell conclusions concern D_Gamma. Both state variants subtract the same complete actual residual. Neither a cell inside the outer state nor its state-area lower endpoint certifies actual full capture. The optional clipped-mass construction in the proof is not an implemented production variant.

The general-block note correctly distinguishes open positive-width rectangle queries from closed queries, including degenerate closed rectangles and atomic ties. Full-witness necessity uses an OLD UPPER bound, while a safe capacity relaxation uses an OLD LOWER bound. The four CDF interval coordinates need not be jointly realizable; this extra relaxation is safe only in the stated upper-bound direction. Shared-cut error has the stated coefficient-weighted adverse signs. The fixed five-source block is uniquely maximal, connected and nonterminal; its real whole-cell deletion lower bound exceeds the equal corner quota. This rules out that quota for general blocks, not the full 7/2 budget.

The thin-prefix note pays the union of all 2K-1 thin joint labels by one increasing pair-residual band. Its clipped-square derivative estimate and odd/even parity formulas are valid. Adding the higher-prefix far-face upper bound is a union upper bound, not an additive area identity. The conditional odd-n four-corner result requires the SAME globally qualified measure in all corner frames. The even-n coarse estimate leaves a positive deficit. No arbitrary-corner or general-block theorem follows.

## One actual independent replay

`computation/review_r079_quadtree.py` ran once under Python 3.15.0: PID 19134, exit 0, confirmed absent; wall 13.8686 seconds. Its persisted result and receipt are in `computation/results/round079-independent/`.

Starting from the two exact 402-source inputs, it rebuilt the actual grid, outside subsets, reflected bboxes and all six residual root enclosures. It then reconstructed each complete base-four prefix trie, requiring all four children at every split and exact root-tile coverage. It recomputed every saved leaf classification using rational derivative-sign minima, legal-radius checks and a horizontal slab union algorithm for residual subtraction, independently of the producer's inclusion-exclusion area routine. No adaptive tree was generated.

The replay checked all 68,912 leaves (45,007 corrected and 23,905 farface), 1,447,152 leaf fields, 274,920 branch tests, 283,049 quadratic minima and 2,412 nonempty horizontal residual slabs. Both variants' complete rational endpoint sums, unknown-area accounting, node counts, stop reasons and case-versus-tile targets agree. Corrected's main state interval is approximately [0.01006487521746,0.01026842028716], with case gap above 10^-4 after the node cap. Farface's interval is approximately [0.01921122504438,0.01927696415011], with case gap below 10^-4 despite one tile missing its allocated target. Unknown leaves remain included in every upper endpoint. Production used five recorded task-worker PIDs; eight was its configured limit.

The exact endpoints are in the result. All lower endpoints are OUTER-STATE area lower bounds, not true halo/deletion lower bounds. R78 qualification and actual-area results were reused by their frozen input/result hashes; no old caps, actual halo union or R78 run was repeated.

## Source and author-run scope

Reviewed final proof pins:

- Quadratic envelope: `b9561e32f5403d6294a7eb4e801df4c63fa76e0a21ee432fcb3d9d8f632e2cb7`.
- General profiles/obstruction: `7b8859301ca4b035c95678ebcc2a6c0f4a6f79ea05bc735112dfcfb9fb2b9ede`.
- Shared thin-joint band: `59f4f8c529ec4c9c791f5ef92a189b5536e35d9fd8cd70f7d09ecfe491956f38`.

The author's quadratic 34-fixture kernel, core 31-subset/14-bbox five-source kernel, 200 CDF queries and eight parity rows were read with their saved outputs/receipts; they were not independently rerun. Root's SymPy source/output/receipt were also read: one identity run, no area or convergence check; it was not rerun. Full proof quantifiers are carried by the analytic arguments, not those finite fixtures. The independent replay source SHA is `1f4898560e64dfa301a14a55be2dadc0c369bfad6ac9ab803c858282bb2d77aa`; result SHA is `805ecaaaf63748fc70504d6cf436c68e77f11d0e145c97fba6fe2fe78d1efa2d`.

Final report SHA: `614119ecf4e2f6a5a6729a1763a14e5899ba892c0d05f46a01325ff7ce7cb9af` (Python version and postprocess wrapper-time corrections applied); report scope reviewed PASS. The component manifest binds final report and component versions. Root's unified integration performs archive byte/link validation. General profile construction, controlled shared-cut weights and the universal continuum area payment remain open; mass error alone gives no area convergence. No R80 statement is included in this review.
