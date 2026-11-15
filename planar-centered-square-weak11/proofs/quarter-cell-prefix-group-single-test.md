# One symmetric four-side input: coarse quarter-cell grouping stays below 7/2

Date: 2026-11-15 (Asia/Shanghai). Root-authorized bounded single-input test for R49. New files only; the frozen R34/R40 and preceding R49 materials are not modified. No commit or push by this agent.

## Precisely which proxy is tested

Partition each coordinate into the intervals [0,1/4], [1/4,3/4], [3/4,1]. Assign mass 1/8 to each of the eight perimeter cells and zero mass to the center cell. Total mass is 1, each outer side band has mass 3/8>1/4, the total corner-cell mass is 1/2 and the center mass is zero.

For a block insertion order π, let B_j be the complete bounding rectangle of cell j and let k be that cell's insertion rank. The tested coarse container is

    R_j=B_j+[-sqrt(k/8)/2,sqrt(k/8)/2]²,
    G_π=union_j R_j.

Only cumulative prefix mass and whole-cell bounding boxes enter this definition. Actual common-capture eligibility is ignored. This is distinct from the original effective insertion proxy U_π, and G_π is not a formula for the true maximal-function halo. A point-insertion square centered inside cell j and having area at most k/8 is contained in R_j, which explains the conservative block cover.

## Exploratory enumeration

Index cells by increasing x-band and then y-band, omitting the center. D4 symmetry reduces the final cell to corner index 0 or edge index 1. The production script enumerated all 5040 remaining-cell orders in each case, 10080 representative orders total. The floating exploratory minima were

    corner last: 3.2254455140552167,
    edge last:   3.1954138576374205.

These numbers are exploratory only. No exact all-order minimum is asserted. No further mass vector or grid geometry was tested.

## One specified order has a rigorous upper bound below 7/2

Use the reproducible order

    π=(5,6,7,3,0,4,2,1),

whose insertion ranks by cell index are (5,8,7,4,6,1,2,3). A separate check gives the exact strict fact

    |G_π| ≤ 58944981652268129571/18446744073709551616
          ≈3.1954138582253653 < 7/2.                    (1)

Thus min_π|G_π| is below 7/2 for this single symmetric four-side-concentrated input. This does not establish a general grouping bound for all inputs with corner mass above center mass, all side-band concentration patterns, or the original effective proxy.

The certificate uses b=32. At rank k the expansion radius has square k/32. Set

    j=floor(sqrt(k·2^(2b−5))),
    J=j if j²=k·2^(2b−5), else j+1.

The radius is at most J/2^b. With common denominator D=4·2^b, all cell endpoints and outward expanded-rectangle endpoints are integers. Integer rectangular-cell coverage computes the union of these eight outward rectangles exactly. That integer union contains G_π, so its area divided by D² is a rigorous upper bound. The numerator before reducing the fraction is 943119706436290073136, with denominator D²=295147905179352825856. Its strict margin below 7/2 is

    5618622605715301085/18446744073709551616.

Only the specified order received this exact upper certificate; no all-order exact lower bound or optimality certificate was computed.

## Actual runs and files

Python 3.15.0, standard library only, one process, no randomness. The exploratory enumeration completed in 0.111126 seconds, exit code 0. The specified-order isqrt and integer-area check completed in 0.008845 seconds, exit code 0. Floating union areas selected the candidate order; they are not used in (1).

Run from the topic directory:

```sh
python3 -B computation/quarter_cell_prefix_group_proxy.py
python3 -B computation/quarter_cell_specified_order_upper.py
```

The two scripts and their actual output files are listed with SHA-256 in `computation/results/quarter-cell-prefix-group-proxy/HANDOFF-MANIFEST.json`. The first output and its CSV ledger retain all exploratory ordering results. The second output is the strict single-order certificate. The general four-side concentration ledger remains open.
