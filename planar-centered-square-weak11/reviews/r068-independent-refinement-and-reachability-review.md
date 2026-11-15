# R68 independent review: strict erosion refinement and source reachability

2026-11-15. The equality W′=W is false. Moreover the old-halo surplus of W′ cannot be removed merely by requiring the full source state to be reachable, even at the receiver being charged. Both conclusions hold for positive inputs in the R67 connected, nonterminal, unique-fat-rank-four-maximum domain. They do not disprove either3.5 state-selection candidate or the original Δ selection problem.

The R67 containments remain valid:

    D⊂W′⊂W⊂H(full)∩Z1,   |W′|−Δ=|W′∩H(remaining)|.

A full label T is eligible and closed, but need not equal the complete source union N_full(x) at the particular receiver. Requiring T to be reachable somewhere differs from restricting its charged region to the cell N_full(x)=T. Neither requirement identifies the actual residual source union. In general,

    N_remaining(x) ⊆ N_full(x)\S

can be strict. Full witnesses that use deleted mass may capture residual atoms that no surviving witness can reach. The actual-bbox core of N_full(x)\S can therefore miss a point already covered by another surviving witness. This is the precise reason that full-state reachability alone cannot turn W′ into D.

## 1. Five sources prove strict W′<W in the candidate domain

Let η=1/4096. Use four deleted sources

    A=(0,0),          m_A=1/4−2η=511/2048,
    B=(9/10,9/10),    m_B=3/4,
    C=(0,3/10),       m_C=η,
    E=(3/10,9/10),    m_E=η,

and one residual source q=(−1/2,0), mass1/2. Here S={A,B,C,E}, m(S)=1, both widths9/10>3/4 and rank4. S is the unique maximum eligible block. Indeed a proper subset of S has mass below1. A crossing subset with B has horizontal span7/5, whose square49/25 exceeds total input mass3/2. A crossing subset without B has mass at most3/4<1. S itself is strictly eligible. S and the eligible pair{q,A} connect the full input, and

    d(full)²−(35/33)m(full)=49/25−35/22=203/550>0.

Consider the open rectangle

    P=(−13/32,−25/64)×(−1/128,1/128),  |P|=1/4096.

At every point of P the complete full source union is exactly N={q,A,C}. The single radius r=27/64 captures these three sources throughout P, and

    4r²=729/1024 < m(N)=3071/4096.

Neither B nor E can be captured by any actual witness at P: their minimum distance throughout P already has square-area cost greater than the mass of the entire input. Thus N is an actual receiverwise complete source state, not a hypothetical eligible closure.

The full-N bbox residual core at b=1/2 misses P, since its x=0 face is farther than sqrt(b)/2. The true residual q square contains all of P. Consequently P lies in the N-label full-bbox erosion ring. Every crossing label has the same residual subset{q}, so all its W′ pieces avoid P. No internal S-only receiver reaches P: one omitting B has total mass at most1/4 and all its source x coordinates are nonnegative; one containing B is excluded by the preceding distance test. Hence

    P⊂W\W′,   |W\W′|≥1/4096>0.

Strict old-halo surplus also occurs in this same example. The rectangle

    Q=(−7/32,−3/16)×(−1/64,1/64),  |Q|=1/1024

lies in both the A singleton receiver and the q singleton residual halo. Its endpoint inequalities are 4(7/32)²<m_A and 4(5/16)²<1/2. The internal A label has no residual core. Therefore

    Q⊂W′∩H(remaining),   |W′|−Δ≥1/1024>0.

The support order A,C,E,B is a monotone chain. The new R68 chain theorem pays its support union exactly at131/50<3; this mechanism example is therefore not in the domain left after the new tree/chain reduction.

This latter surplus uses an internal label; it is not evidence that the complete state at Q equals{A}. In fact the A singleton is a genuine complete source state near A, but its receiver rectangle extends into cells whose complete source set is larger. The six-source construction below addresses the stronger receiverwise completeness question directly.

## 2. Six sources retain surplus even on the actual full-state cell

The author's stronger construction uses S of total mass1:

    p0=(3/10,3/10),   m0=1/10,
    p1=(6/5,6/5),     m1=449/500,
    p2=(3/10,6/5),    m2=1/1000,
    p3=(6/5,3/10),    m3=1/1000,

and residual q=(−1/50,3/10), mass1/10, and z=(7/50,7/50), mass1/1000. The whole input is positive, eligible-connected and nonterminal, with unique maximum S of rank4 and both widths9/10. Independent complete63-subset enumeration verifies the maximum and all ties; the strict nonterminal margin is17637/55000.

Let

    P6=(1399/10000,1401/10000)²,  |P6|=1/25000000.

The four S positions are the corners of a side9/10 square. Their minimum support-tree cost is27/10, with tree bound37/10; their uniform support union has area361/100. Thus the final example is not paid by the new tree or filled-support sufficient tests. These support lower bounds are not lower bounds for W′.

The actual complete full state throughout P6 is N={p0,q,z}, mask49 in the displayed order. The common radius161/1000 captures these three atoms and has area25921/250000<201/1000=m(N). Any other S atom costs at least112338801/25000000>m(full) to capture anywhere on P6, so no full witness can add it.

The residual subset N\S={q,z} has mass101/1000. Capturing q from P6 costs at least2556801/25000000>101/1000. Thus its actual tight-bbox residual core excludes the whole patch. The z singleton residual halo nevertheless covers P6. More precisely the actual residual complete source state there is just{z}, mask32, a proper subset of mask48=N\S. Hence

    P6⊂[E_N(m(N))\E_(N\S)(m(N\S))]∩H(remaining).

The patch also avoids H(S): the p0 singleton has mass1/10 below the same capture-cost lower bound; every other S atom is already excluded even with all S mass available. Therefore this surplus cannot be removed by discarding internal H(S) overlap. It comes from a genuinely actual full source state and a different, smaller actual surviving source state. No full halo area or3.5 failure is claimed for this six-source input.

## 3. Independent computation and exact evidence

`review_r068_state_geometry.py` uses only the frozen independent R67 twelve-layer grid function and the independent R66/R64 math helpers. It imports no producer. It independently reconstructs all coordinate-tight closures, all eligible labels meeting S, both residual core geometries and the full/residual halos. Per-label bitsets calculate unions of actual rectangle differences. The five-source input has nondyadic rational coordinates, so the denominator is D=lcm(2^bits, all coordinate denominators); integer root brackets then give directed rational receiver rectangles. No decimal geometry is used.

The actual five-source run used Python3.13.18,72 bits and a bytewise minimal archive-layout fixture in dataTemp. It checked all31 nonempty subsets,15 closed bboxes,10 eligible closures and9 labels meeting S, on2401 compressed cells. Runtime0.028235seconds, exit0. All six grid inclusion violations were exactly zero. The saved rational intervals give the following display values with L=1:

| Quantity | Approximate interval value |
| --- | ---: |
| W |1.1479961034993444|
| W′ |1.0836349027370908|
| Δ |1.0319937836244355|
| W\W′ |0.06436120076225359|
| W′∩H(remaining) |0.05164111911265539|

The lower endpoints strictly exceed the two analytic patch-area certificates1/4096 and1/1024. Both W and W′ are far below3.5L here; strict refinement is not a failure of either selection candidate.

`review_r068_reachable_source_patch.py` separately executed the six-source63-subset and whole-patch rational check. It found24 closed bboxes and13 eligible closures, unique maximum mask15, and verified every rational radius/capture/exclusion bound above. Runtime0.001471seconds, exit0. This second execution is a patch and source-domain certificate, not a full halo union calculation. The weaker four-source thin toy considered during analysis was canceled before any run or archive and is not part of this delivery.

The first construction originated from the compute agent's fat three-source seed, independently padded and proved here; the six-source construction is the core author's stronger example, independently checked here. Source arithmetic and actual runs are kept distinct. The producer N128 family result and the four final author proof notes are reviewed in the accompanying final proof/implementation addendum; this file records the two independent mathematical obligations above.

## 4. Portable reproduction

From an integrated archive checkout:

```sh
python3 -B planar-centered-square-weak11/computation/review_r068_state_geometry.py \
  --output /Volumes/data/CodexTemp/T2/r068-replay-new/five-source.json
python3 -B planar-centered-square-weak11/computation/review_r068_reachable_source_patch.py \
  --input planar-centered-square-weak11/computation/inputs/r068-six-source-reachable-state-patch-review-input.json \
  --output /Volumes/data/CodexTemp/T2/r068-replay-new/six-source-patch.json
```

The first entry infers the repository root; `--archive-root` may name the repository or topic directory. Its adjacent default spec pins all archive-relative dependencies, including its case input, and there is no old-workspace fallback. `--spec` can select an explicit pinned spec. The second entry takes an explicit input path and records its SHA. Both use only the Python standard library and require a new JSON output below mounted dataTemp. Python `-O` is rejected. Historical R67/R66 files and pins are unchanged.

The review proves and verifies strict state refinement and the persistence of old-halo overpayment. Controlling the size of the full spatial union, or selecting a maximum block with state cost≤3.5L, remains unresolved by these examples.
