# R70: mass-sensitive necessary mixed-receiver envelope

This module computes an envelope M depending only on the mass and support of a selected maximum block S, with L=m(S). It is a necessary outer carrier for mixed receivers, not the erosion state W or Wprime. The three fixed pilots are the R69 balanced eight-source construction (S four corners, spacing 9/10), a four-source S-only square of spacing 99/100, and the historical s00-c00-n16 maximum block. The uniform control has no ambient A or residual halo. No new ambient configuration, random search, or LP was constructed.

For every nonempty proper S-bbox-closed U, let a=m(U), r(x)=max over p in U of ||x-p||infinity, and c(x)=max over p in S of ||x-p||infinity. All bbox closures are included, including those which are not internally eligible. A square capturing exactly U from S has this closure property. The relaxed carrier allows a candidate total mass m in (a,L], r<sqrt(m)/2 and c>sqrt(L+m-a)-sqrt(m)/2. It does not require U to be exactly captured by that candidate square, which is one source of relaxation.

There is no mass grid. The function f_a(m)=sqrt(L+m-a)-sqrt(m)/2 decreases up to (L-a)/3 and then increases. On the feasible interval its infimum is attained in the closure at

q=max(a,4r²,(L-a)/3), provided q<=L.

Strict feasibility needs r<sqrt(L)/2, a<L, and the strict c inequality; q=L is retained only as an outer boundary. For rational c,r the comparison contains no numerical root: put z=L+3q/4-a-c². Then c>f_a(q) iff z<0 or c²q>z². This uses c>=0 and q>0. In integer units, with C2=2c, R2=2r, q3=max(3a,3R2²,L-a), and z4=4(L-a)+q3-C2², the last test is z4<0 or 4 C2² q3>3 z4². The lower branch automatically gives the piecewise g(a) minimum stated in the theory interface.

Each spatial cell has exact rational lower and upper bounds for r and c from its coordinate intervals and the two support bboxes. The minimized threshold is nondecreasing in r. A cell is certified inside if one label passes using c_min,r_max with r_max<sqrt(L)/2; it is excluded if every label fails using c_max,r_min (or is beyond the maximum receiver radius). Every other cell remains undecided and belongs to the outer area only. The bounded domain is bbox(S) enlarged by an integer outward bound for sqrt(L)/2. Every nonempty candidate U has its receiver in that domain.

The initial and final budget was 8000 classified nodes, maximum depth 10, for each of three inputs. Pending or depth-limited cells are preserved as undecided. There was one actual batch and no refinement run. Each gzip certificate includes the whole covering leaf partition, including excluded and undecided leaves, integer scale, support bbox, domain, all label bboxes/masses, and classification legend. The integer cell coordinates use a rational denominator clearing both source coordinate and mass denominators, followed by 40 dyadic bits. The certificates have deterministic gzip timestamps. A verifier can reconstruct the midpoint subdivision tree from domain/depth/cells, confirm complete disjoint coverage, and test each classified leaf against every label without running a mass search. Unresolved area is reported explicitly, not discarded as a rounding error.

I=H(S) uses the pinned R64 root/rectangle library. Inner and outer rectangles have exact integer endpoints. Reported costs include |M|, |I|, |M\I|, and the S-only budget |I|+|M\I|, with all overlap treated as a union/difference. Where an ambient A is supplied, |M\(I union R)| is a slack diagnostic and does not redefine the S-only envelope. Historical true P and Delta are reused byte-bound from the final R69 result. Balanced-eight true P and Delta were evaluated in this batch. No old W/Wprime geometry was recomputed.

From the repository root, reproduce one selected case with:

```sh
python3 -B planar-centered-square-weak11/computation/round070/mass_sensitive_receiver.py \
  --archive-root . \
  --spec planar-centered-square-weak11/computation/round070/three-pilot-spec.json \
  --job uniform4-t099 --out-dir /Volumes/data/CodexTemp/T2/r070-one-replay
```

The exact original fixed batch, including per-child PID/exit/absence receipts, is reproducible with:

```sh
python3 -B planar-centered-square-weak11/computation/round070/run_three_pilot.py \
  --archive-root . \
  --spec planar-centered-square-weak11/computation/round070/three-pilot-spec.json \
  --producer planar-centered-square-weak11/computation/round070/mass_sensitive_receiver.py \
  --out-dir /Volumes/data/CodexTemp/T2/r070-three-replay
```

The runner uses three Python subprocesses, a sampled 1 GiB RSS threshold, a 180-second wall threshold, and verifies each owned PID is absent after waiting. Original absolute commands are preserved in receipts. Actual execution used Python 3.13.18 at /Users/zzh/.venvs/main/bin/python3, 0.991167167 seconds batch wall, and peak sampled per-child RSS 28032 KiB. All three exited 0 and were absent. No MATLAB was launched.

Uniform spacing 99/100 gives a certified S-only budget lower bound exceeding 3.59729L, so this envelope cannot yield a universal 3.5L budget. This is an obstruction to the envelope, not a counterexample for true deletion or W/Wprime. The other two selected S-only budgets are below 3.312L, while their ambient mixed-unpaid envelope areas are much larger than true P. Moderate cell enclosures already decide both the obstruction and the slack; no finer geometry was needed. The core report supplies the separate analytic argument and interpretation.
