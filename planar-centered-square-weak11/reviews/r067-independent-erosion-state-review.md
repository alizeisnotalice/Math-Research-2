# R67 independent proof and implementation review

2026-11-15. Reviewer workspace: pcw-sol-cover-structure. The support-rank reduction, the full erosion-state containments and the prescribed positive-padding obstruction pass independent review. The W and W′ selection statements remain candidates; neither two passing inputs nor rank reduction proves the original Δ selection or a general C2 bound. R66 files remain frozen. No producer source was edited, no MATLAB was launched, and no commit/push was made here.

## General analytic interface

Every actual strict centered-square witness has captured mass at most L and half-side less than sqrt(L)/2. A receiver lost on deletion has no surviving witness and therefore has a witness meeting S. This places it in the union Z1 of equal-radius source squares over S. Distinct positive support positions are the relevant count: rank≤3 gives Δ≤3L and also |W|≤3L. This is a sufficient reduction, not a claim of a sharp support bound. Coincident positions must first be merged as the note assumes.

For every eligible closed label T, its full-mass receiver rectangle lies in the actual full halo. Its same-full-bbox core at mass b=m(T\S) lies in the actual residual halo. The single-label ring area is at most a=m(T∩S), by the reviewed R66 locally absolutely continuous area function with derivative between0 and1. The b=0, point-bbox and critical/empty-core cases are covered. Taking the **union** of rings is essential; the R66 repeated-capacity obstruction still prohibits summing a over labels.

The actual residual subset B=T\S has a tighter bbox. Its receiver core contains the same-full-bbox core and still lies in the residual halo. It may have another center, and may protrude outside the full receiver rectangle. Hence W′ must use an actual rectangle difference; pretending it is a concentric erosion would be wrong. All complete source states needed for a true deletion receiver occur among the included eligible closed labels. This proves

    D⊂W′⊂W⊂H(full)∩Z1.

Because D=H(full)\R, both |W|=Δ+|W∩R| and its W′ counterpart follow as exact measure identities. Their surplus is old halo, rather than artificial points outside the full halo. The whole remaining halo is not part of either state's construction. It is used separately to diagnose Δ and surplus. Closed labels not individually reachable as complete N-states may still overpay; neither identity establishes a selection theorem.

The rank and fat/connected/nonterminal domain restrictions in the candidate statement have the correct quantifiers. A refutation must cover every maximum-block tie in a qualifying component and certify a strict **lower** bound above3.5L. A single failed block or an outer upper enclosure above3.5L is insufficient. The unconditional filled C1 area bound remains below4L, which does not close the3.5 target.

## Positive padding obstruction

All nine support-grid masses are strictly positive, sum to1, and the heavy corner has mass4/5−7·2^-20>3/4. Every old proper subset has mass below1. A crossing subset omitting the origin has mass below1; one containing both origin and outside atom has span squared above total input mass. Thus S is the unique eligible maximum, rank9, fat in both axes. S and the heavy-corner/outside pair give connectivity. The nonterminal margin6982433/34603008 is positive.

The prescribed source-neighborhood radius is sqrt(θL)/2, with θ=3/4. The Cartesian support grid makes its neighborhood union exactly the old filled bbox expansion: adjacent one-dimensional intervals overlap, then their product fills C. The heavy singleton's exterior area is φ(h)=3h/4−3/16−sqrt(3h)/4. On h≤t≤4/5, 0≤φ′(t)≤3/4, so the R66 strict rational lower certificate loses at most21·2^-20/4. The resulting strict excess is982963/524288000>0. The residual square lies inside C, so no omitted residual credit changes this sign.

This disproves the proposed θ=3/4 support-neighborhood repair with that common radius, even after the low-rank reduction. It does not disprove mass-sensitive W/W′, arbitrary other radii or thresholds, or true deletion. The final proof wording is bound by the final proof pins in the accompanying manifest.

For this particular input, all same-full-bbox residual cores vanish: any crossing label spans at least1/8 while sqrt(b)=1/32; internal labels have b=0. The residual singleton square is contained in the heavy internal receiver. Therefore W=H(full). The same internal receiver also covers the only parts removed by the tight residual cores, so W′=W here too. Its surplus is exactly the old singleton halo ε, not a low-container hole. This is an example-specific conclusion.

## Independent algorithm and actual results

I read the frozen `whole_erosion_state.py` SHAab0e7dba2993869d69871384e9f95f0f095563b453dbde897561097d06251ccf, the two-control runner, and the positive-padding builder. Label enumeration includes **every** strictly eligible full-input closed bbox meeting S; it has no θ filter. b and the actual B bbox are recomputed per label. The producer's difference calls use full-inner minus core-outer for the lower ring and full-outer minus core-inner for the upper. Full residual halo computations occur only after the un-discounted state lists have been assembled. Its per-label capacity and exact inclusion checks have the stated scope. A passing input is not universal source completeness for another input.

My new `review_r067_erosion_states.py` imports only the frozen independent R66 grid helper interface and its independent R64 geometry helper. It does not import any producer module or call strip unions or rectangle-subtraction routines. It independently enumerates coordinate-tight bboxes, uses exact integer binary root brackets and places twelve layers on a compressed2D difference grid: full inner/outer, residual inner/outer, per-label full receivers, full-bbox cores, actual-subset-bbox cores, and rank-neighborhood inner/outer. Distinct label bits compute each union of differences pointwise, preserving overlap and nonconcentric cores. All endpoints and area comparisons are integer/Fraction operations.

Both actual prescribed inputs were replayed from a six-file minimal archive-layout fixture in dataTemp. Its files were copied bytewise from the compute final manifest destinations and the two already published independent helpers. The entry checks all six SHA256 pins before importing the core or opening case data, checks input/output/source chains, and has no old-workspace fallback. Only the portable spec and entry themselves are supplied outside that fixture. Python3.13.18, two workers,0.323156seconds, exit0. The owned entry PID5874 was checked absent. In total71 label rows and76,026 compressed cells were checked.

| Input | Full closed bboxes | Eligible labels meeting S | W upper / L | W′ upper / L | Δ upper / L |
| --- | ---: | ---: | ---: | ---: | ---: |
| Positive9+1 padding |46|17|1.00066538634|1.00066538634|0.99968882384|
| Frozen R64 final |274|54|1.29922593162|1.29922593162|1.23599350975|

Both inputs have unique S511, rank9, both fat axes, one connected component and the strict nonterminal condition. Every saved Fraction endpoint matched: W, W′, full halo, residual halo, Δ, both old-halo overpayments, refinement interval and every label's mass/span/core/ring/capacity record. All six direct2D inclusion violations were exactly zero, including outer W outside the outer rank neighborhood. Both states are certified below3.5L on both inputs. On R64 the two state *intervals* agree; this alone does not prove equality of their exact regions or areas. Its old-halo surplus is approximately0.063232421875L. Tiny positive refinement upper endpoints reflect directed enclosure slack.

The padding check independently enumerated all1023 nonempty subsets, checked positive masses and unique maximum, and checked Cartesian interval overlap and the rational982963/524288000 sign. The original fixed-θ padding budget execution was source/receipt-bound, not independently rerun here; its analytic failure is established above. Author and producer padding inputs have different serialization/order/metadata, but the same exact position-mass triples and deleted set. The recorded hashes are kept distinct.

Root canceled a third N32 fat-family halo probe after the prior C1 area bound settled its outcome. Only the rational identity (1+217/256)(1+7/8)=7095/2048<7/2 is checked by this entry. No N32 halo, extra input family, θ search, optimization or MATLAB run is claimed.

## Portable reproduction and delivery scope

In a clean published checkout, run:

```sh
python3 -B planar-centered-square-weak11/computation/review_r067_erosion_states.py \
  --output /Volumes/data/CodexTemp/T2/r067-replay-new/result.json --workers 2
```

`--output` is mandatory and must name a new JSON file below dataTemp; a sibling receipt is saved. `--archive-root` may name the repository root or its topic directory; default is inferred from the installed entry. The default pinned `r067-independent-replay-inputs.json` lives beside the entry; `--spec` can select an explicit spec. Python3.10+ standard library is sufficient. Python `-O` is rejected because the frozen mathematical checks require assertions. All case paths and dependency paths in the default spec are archive-relative.

The supplementary manifest records the exact entry/spec/proof/source pins, actual independent result/runtime/process receipt, and fixture-copy provenance. Saved producer receipts and author arithmetic are hash-bound; they are external saved executions, not independent reruns by this reviewer. The two proofs and source scope were audited directly. No blocking mathematical or implementation issue remains within this specified interface. W/W′ selection and the original general Δ/C2 goals remain open.
