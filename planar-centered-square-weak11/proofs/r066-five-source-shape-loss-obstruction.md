# R66: a five-source obstruction to discarding the actual low-container shape

2026-11-15. This exact example refutes the shape-free compensated-high selection condition when it is asserted for **all** connected nonterminal inputs. Its unique maximum block has zero height and is already paid by R63. It does not refute the same condition on the unresolved all-fat branch, the actual shape-sensitive budget, the original maximum-block deletion conjecture, or any general C2 bound.

## 1. Input and complete high states

Let S consist of four sources at(0,0),(1/4,0),(1/2,0),(3/4,0), each of mass1/4. Add one source at(67/64,0) of mass1/64. Thus

    m(S)=L=1,  m(A)=65/64,  d(A)²=4489/4096,
    d(A)²−(35/33)m(A)=2537/135168>0.

S is the unique maximum eligible block. Any proper old subset has mass≤3/4, any set using the outside source and at most three old sources has mass≤49/64<1, and the entire five-source set has mass65/64 but span squared4489/4096>65/64. S itself has span squared9/16<1. The full S is an eligible edge and the three right old sources together with the outside source form another eligible edge, so A is connected.

At theta3/4 the only eligible closed high sets touching S are S and

    T={(1/4,0),(1/2,0),(3/4,0),(67/64,0)},
    m(T)=49/64, sqrt(m(T))=7/8.

Indeed three old sources have mass exactly3/4 and therefore are in the low layer. A different three-old-source group with the outside point has minimum coordinate0 and is not eligible. The complete receiver rectangles are

    E_S=(1/4,1/2)×(−1/2,1/2),
    E_T=(39/64,11/16)×(−7/16,7/16).

Their horizontal projections are disjoint. The low container is

    C=(−sqrt3/4,3/4+sqrt3/4)×(−sqrt3/4,sqrt3/4).

The residual halo is the singleton square

    R=(63/64,71/64)×(−1/16,1/16)⊂C,

so its low old-area credit is K=1/64. The inclusion follows, for example, from sqrt3>23/16 and the displayed endpoints. This uses actual residual source mass and an actual halo.

## 2. Exact area ledger

Both high rectangles have their entire horizontal projections in C. Their parts outside C are disjoint top/bottom strips, and already disjoint from R. Therefore

    U=|E_S\C|+|E_T\C|
      =(1/4)(1−sqrt3/2)+(5/64)(7/8−sqrt3/2)
      =163/512−21sqrt3/128.

The compensated high term is

    G=U−K=155/512−21sqrt3/128.

Compared with the shape-free margin delta=7/4−sqrt3,

    G−delta=(428sqrt3−741)/512>0,                         (1)

because 3·428²−741²=471>0 and both compared numbers are positive. There is only one maximum block, so no choice of another maximum block repairs this particular shape-free condition.

The actual shape budget remains small:

    |C|=3/4+3sqrt3/8,
    B=|C|−K+U=539/512+27sqrt3/128
      <91/64<3/2<7/2.                                    (2)

Here sqrt3<7/4 follows from3<49/16. The strict failure in(1) is caused by replacing the actual low-container area by its much larger worst-case value. The maximum block has widths(3/4,0), so the thin-axis theorem already applies. A valid main research reduction must perform that branch first and reserve the unresolved selection assertion for components whose every maximum block is fat.

## 3. Actual checks and scope

`computation/build_r066_five_source_shape_loss.py` actually checked all31 nonempty source subsets with exact integer/Fraction arithmetic, including strict eligibility, unique maximum mask15, high closed masks15 and30, actual connectivity, the nonterminal margin and the positive integer471. Python3.15.0, runtime0.00034seconds, exit0. The exact minimal input is `computation/inputs/r066-five-source-shape-loss.json`, SHA `b9e9888e8d9f2f12525b41dc85cfeda07f5554527fae838dc59ab4706fdd1940`.

The frozen production budget and erosion programs separately ran this input and obtained directed full-geometry intervals compatible with all displayed closed forms. Root's independent Wolfram v2 run additionally enumerated31 source subsets and used a two-dimensional algebraic coordinate-cell union calculation to verify the high area, old credit and budget, together with full deletion and hole formulas. Its final source is `r066-five-source-independent-v2.wl`, SHA `1085a949bf4acdd0abef07857ec33e3aa4196bdcf79f5e2aec051cce89b51da0`; Mathematica15.0.0 ARM,0.569151seconds kernel time/1.3312seconds wall, exit0 and owned PID verified gone. The earlier v1 AssociationMap export error is not counted as a successful run.

The elementary geometry above is the proof. The exact runs provide independent formula/implementation checks. No statement about a general fat-branch counterexample or a value of C2 follows from this example.
