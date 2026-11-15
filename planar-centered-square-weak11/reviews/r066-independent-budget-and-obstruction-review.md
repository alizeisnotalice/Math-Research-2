# R66 independent review: exact budgets, erosion, and the fat obstruction

2026-11-15. Reviewer workspace: `pcw-sol-cover-structure`. Status: PASS for the four final analytic notes and the specified frozen implementations, with the fixed θ=3/4 envelope selection explicitly refuted. The original Δ-based cheap-deletion question and general C2 bound remain open. Author sources were read only; R61–R65 files were preserved. No MATLAB was launched and no commit/push was made by this reviewer.

## Decisive new result

The three-source input with masses `(1/5,4/5,1/1024)` at `(0,0)`, `(1023/1024,1023/1024)`, `(1151/1024,1023/1024)` has unique maximum eligible block S consisting of the first two atoms, L=1. Both S widths are 1023/1024. The input is one eligible-connected, nonterminal component. Independent enumeration of all seven nonempty subsets and independent full-halo/residual-halo geometry confirm all these domain conditions.

The complete high closure list is masks 2,3,6. With C the θ=3/4 filled bbox expansion and R the residual singleton halo, R⊂C and its area is ε=1/1024. The heavy singleton square Q gives

    U ≥ |Q\C| = 33/80 − sqrt(3/20).

Thus

    B ≥ α²+93/80−ε+α sqrt(3)−sqrt(3/20),  α=1023/1024.

The signs are certified with rational arithmetic, independently of geometry intervals: sqrt(3)>433/250 and sqrt(3/20)<97/250 follow by squaring positive rationals. Substitution gives

    B−7/2 > 246397/131072000 > 0.

The stronger shape-free compensated condition also fails, by at least 707/128000 using the same rational substitutions. Since S is unique, neither failed condition can be repaired by choosing another maximum block. This disproves the fixed θ=3/4 **actual envelope budget selection**, even in the all-fat branch.

The independent ten-layer 2D difference-grid calculation matches every saved Fraction endpoint and obtains B_lower≈3.50309040474105>3.5 while Δ_upper≈0.999688668942791<1. The proof's separate analytic Δ bound is valid: the pair receiver rectangle, heavy singleton square and residual square lie in the enlarged square Q_b about the heavy atom. Together with the light singleton receiver and E_S, this yields Δ≤1+(1−α)²=1+1/1048576. Therefore the obstruction is to filling C, not to the original Δ-based selection statement.

The final main proof directly rewrites section3 to state the refutation; section5 and section6 consistently retain the original Δ question as open. I reread these final changes and the three-source note's final execution paragraph. The three-source adapted producer input and the author's minimal input are distinct serialized artifacts; their different hashes are intentional, not an input mismatch.

## General proof audit

The complete source set N(x) uses one actual strict common-radius witness, as in the reviewed R64 interface. For a finite positive atomic input, witness capture sets at fixed x are nested by radius; the finite union is attained by one of those witnesses. Bbox completion preserves coordinate spans, captures every atom inside the bbox, and can only increase mass. Consequently all actual source states are among the eligible closed bboxes used in Pθ, although not every eligible closure need be a separately reachable N-state.

For x in H\R, its full source set meets S. If its mass is at most θL, its strict witness radius puts x in Cθ; otherwise x belongs to Pθ. Since Pθ⊂H, the exact exterior identity D\Cθ=Pθ\(Cθ∪R) follows. Subtraction then gives

    Bθ=Δ+|Cθ\H|=|(Cθ∪Pθ)\R|.

Hence Bθ is nondecreasing in θ. Searching θ for this exact filled-container budget favors θ=0; it does not prove a stronger selection theorem. At θ=0 a degenerate open bbox may be empty. At θ=1 there are no high closures since L is maximal; the proof still holds. The producer intentionally accepts only 0<θ<1 and includes S in the high list. Omitting S in general is invalid; the singleton example in the main proof correctly diagnoses that omission.

The directed intervals in the main proof have the correct subtraction directions. High membership is the exact strict comparison m(T)>θL, with equality low. Full and residual halos must each use their own masses and closures. The direct interval for |(C∪P)\R| has width at most the sum of the C, P and R uncertainty bands. The producer instead adds separate directed area/credit/high endpoints; this is valid but potentially wider. Its compatibility interval for B−Δ−hole is a finite consistency check, not the proof of the identity.

For erosion, E_T(b) uses the same full-T bbox with b=m(T\S), not an invented residual core. Every point of this rectangle is covered by the actual residual halo. With F(t)=(sqrt(t)−u)_+(sqrt(t)−v)_+, F is zero below max(u,v)² and has derivative 1−(u+v)/(2sqrt(t)) in [0,1] above it. Continuity and local absolute continuity include the threshold, b=0 and u=v=0. Thus a single unpaid rectangle has area at most a=m(T∩S). The high **union** of the erosion rings is a safe upper bound; summing a across labels has no such global source budget.

The rim calculation and nearly critical receiver-area factorization are correct. Combining worst-case container and rim bounds only returns 4L; the half-L saving is relative to these worst-case bounds, not a requirement for literal old-halo credit on every input.

The two-row fat family has the claimed exhaustive closed-box classification. A closed box using the external source has a one-row suffix, a two-row suffix, or just that source. The formula d_x²−m=(3/4)z+(49/64)z² excludes mass≥1 crossing competitors, while the selected negative-z tails are strictly eligible. The unique maximum is S, it is fat, the full input is connected/nonterminal, and the selected capacity sum is exactly (41/64)(7N/32−1), unbounded in N. This refutes a sum-of-intersections budget, not the actual union budget or deletion. The general argument is analytic; the finite root census does not replace it.

The separate five-source example correctly has unique maximum S of widths (3/4,0), high masks15 and30, K=1/64, G−(7/4−sqrt(3))=(428sqrt(3)−741)/512>0 because 3·428²−741²=471. Its actual budget is (539+108sqrt(3))/512<91/64. It refutes only the unqualified shape-free condition. Its stated scope remains correct when read as the conclusion of that example, separately from the new fat three-source result.

## Independent implementation and actual runs

The frozen producer budget code is SHA211d33e21e47e6e14290e4ddbc3b3a36e247131de1bc53114091b7766d0fe8b6; frozen erosion code is SHA4dc4fb494736c40334359f9950f052ee2e4a53b503883ffbc3f9bc0da69b1361. Both rely on the reviewed R64 module SHA4ef61b815b9c4a9a8972cdaf31716ae7917b59243a29fb8697c8042137788246. I read both R66 programs and checked strict eligibility, inclusion of S, full/residual mass recomputation, exact high selection and directed set subtraction. Their union/strip algorithms were not imported by my independent replay.

`review_r066_grid_budgets.py` imports only my R64 independent helper. It enumerates coordinate-tight bboxes independently, brackets roots by exact integer binary search, and places full, residual, low-container and per-label high/core rectangles on a ten-layer compressed 2D difference grid. Distinct label bits preserve unions of ring differences without calling the producer's strip/merge/subtract-rectangle routines. Every resulting metric uses integer/Fraction arithmetic. The saved exact producer endpoints are compared for equality; I did not merely compare decimal summaries.

The actual complete replay used Python3.13.18, twelve local worker processes, and took 2.73224 seconds. It checked 74 input records and 201 proper maximum blocks: all 72 R63 cases/199 blocks, the R64 final block, and the old thin five-source control. It checked every saved maximum tie, complete full/residual closure geometry, low old credit, high unpaid union, raw high exterior area, low hole, Δ, B, identity residual, high compensation scores, and all 1103 individual high-label receiver/erosion/outside/unpaid records. It processed 4,295,571 compressed grid cells in total, with 99,225 at most in one block. All comparisons passed.

Within the original R63 sample, all 199 blocks have B≤3.5L. The stronger condition has 181 passes and 18 strict failures. These are finite sample facts. In particular, the isolated fat failure s05-c05-n20/S557056 is an entire independently terminal component of a disconnected input, not an obstruction in the intended fat connected/nonterminal branch. The new deterministic three-source control is recorded separately and supersedes any temptation to generalize the sampled fat passes.

All 1103 erosion labels matched. In 27 labels, the outward **upper enclosure** of ring area slightly exceeds a; the exact mathematical ring bound remains valid and all lower endpoints are ≤a. Such enclosure slack is not a counterexample to the single-rectangle theorem. Equality between ring-union and actual-unpaid interval values on some inputs is also not a general identity.

`review_r066_critical_controls.py` separately reran the five-source and new three-source controls with the same independent grid algorithm, enumerated all31 and7 nonempty subsets, checked the radical signs using rational inequalities, checked 375 exact rational-square erosion boundary cases, and performed a generic-coordinate bbox census for the smallest N=32 fat repetition input. That census found1681 closed bboxes,1012 eligible closures, unique maximum mass1, six selected tails and capacity sum123/32. Actual runtime0.13265 seconds. The old201-block output is unchanged; the new fat example is a separate supplement.

Saved root budget/erosion receipts, Wolfram v2 receipt/output and the twelve-case fat-family census were inspected and source-bound. Those are saved external executions, not independent reexecutions by this reviewer. Wolfram v1's export failure is excluded. The root fat-family source was read, but I did not repeat its entire twelve-N census. The final binding output records hash-chain/receipt checks and the manifest binds the exact final four proof pins.

Reproduce in this workspace with `/Users/zzh/.venvs/main/bin/python3 -B computation/review_r066_grid_budgets.py --workers 12 --out computation/results/round066-independent-grid-replay.json`, then `/Users/zzh/.venvs/main/bin/python3 -B computation/review_r066_critical_controls.py`. The source paths are currently absolute local provenance paths; the integration owner must preserve or remap them when publishing.

No blocking issue remains in the final four notes. What is now disproved is the fixed θ=3/4 filled-container envelope selection route, including its stronger shape-free condition in the fat branch. A true Δ-based maximum-block deletion theorem or a general C2 upper bound has not been proved or refuted here.
