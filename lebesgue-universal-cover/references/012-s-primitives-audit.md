# Gibbs 2018 Definitions 4--7: computable first arc cuts

Sources reviewed:

- Gibbs 2018, pp. 6--8, Definitions 4--7 and Propositions 3--6;
- Baez--Bagdasaryan--Gibbs 2015, p. 7, which fixes the dodecagon label order
  `A3,A2,B3,B2,C3,C2,D3,D2,E3,E2,F3,F2`.

This note reconstructs the three **Sprague-like** cuts only.  It does not
cover the reflection-dependent cuts `AH` and `EH`.

## Label map and coordinates

Use the base hexagon `H`, its normal vectors
\(n_j=(\cos(j\pi/3),\sin(j\pi/3))\), and the rotated-hexagon normals
\(m_j=(\cos(\pi/6+\sigma+j\pi/3),\sin(\pi/6+\sigma+j\pi/3))\).
Every supporting line has equation `normal . x = 1/2`.

If `I(a,b)` denotes the intersection of the two supporting lines, the labels
needed by Definitions 4--6 are

\[
\begin{array}{c|cccccccccccc}
 &A_3&A_2&B_3&B_2&C_3&C_2&D_3&D_2&E_3&E_2&F_3&F_2\\ \hline
 &I(n_2,m_1)&I(n_1,m_1)&I(n_1,m_0)&I(n_0,m_0)&
 I(n_0,m_5)&I(n_5,m_5)&I(n_5,m_4)&I(n_4,m_4)&
 I(n_4,m_3)&I(n_3,m_3)&I(n_3,m_2)&I(n_2,m_2).
\end{array}
\]

For independent normals `a,b`,
\[
I(a,b)=\frac{1}{2(a_1b_2-a_2b_1)}(b_2-a_2,\ a_1-b_1).
\]

The map agrees with the textual incidences: `C3C2` and `E3E2` are rotated
hexagon sides, while `D3C2`, `E2F3`, `GD1`, and `D1E3` are original-hexagon
side subsegments.

## Definitions 4--6 as finite operations

Let \(\operatorname{proj}_a z=z-(a\cdot z-1/2)a\), orthogonal projection
onto the support line with normal \(a\).

| Region | endpoints and circular arcs |
| --- | --- |
| \(C_S\) | \(K=\operatorname{proj}_{m_5}(F_3)\); \(G\) is the unique point of \([D_3,C_2]\) with \(\lvert G-F_3\rvert=1\). Its boundary is the minor `F3`-centred arc \(K\to G\), followed by segments \(G C_2\) and \(C_2K\). |
| \(E_S\) | \(I=\operatorname{proj}_{m_3}(B_3)\), \(H=\operatorname{proj}_{n_3}(C_3)\); \(J\) is the common point of the unit circles centred at \(B_3,C_3\) nearest \(E_2\). The boundary uses arcs \(I\to J\) centred at \(B_3\), \(J\to H\) centred at \(C_3\), and segments \(HE_2,E_2I\). |
| \(A_S\) | \(P=\operatorname{proj}_{n_2}(G)\), \(Q=\operatorname{proj}_{n_1}(E_3)\); \(X\) is the common point of the unit circles centred at \(G,E_3\) nearest \(A_1\). The boundary uses arcs \(Q\to X\) centred at \(E_3\), \(X\to P\) centred at \(G\), and segments \(PA_1,A_1Q\). |

The “nearest” rules choose the root pictured by Gibbs.  The other circle-circle
root is outside the base hexagon for the checked parameters.

## Area formula and actual checks

For an oriented unit-radius circular arc from `a` to `b`, centred at `c`, let
\[
\Theta=\operatorname{atan2}((a-c)\times(b-c),(a-c)\cdot(b-c)).
\]
Its contribution to \(\int(x\,dy-y\,dx)\) is
\[
c\times(b-a)+\Theta.
\]
Each line segment `a -> b` contributes \(a\times b\).  Thus one half of the
sum of these contributions gives the area of each displayed curvilinear
polygon.  This is an executable exact geometric formula up to the elementary
trigonometric evaluations and branch rules above.

`reconstruct_s_primitives.py` checked every stated unit distance and segment
membership at \(\sigma=0,1.3,1.5494\) degrees.  The corresponding cut areas
were:

| \(\sigma\) | \(|C_S|\) | \(|E_S|\) | \(|A_S|\) |
| ---: | ---: | ---: | ---: |
| \(0^\circ\) | numerical zero | numerical zero | 0.001161753185550909 |
| \(1.3^\circ\) | 0.000000304251975370 | 0.000000597803623060 | 0.001170265313731364 |
| \(1.5494^\circ\) | 0.000000516575992459 | 0.000001011646495025 | 0.001173792300040971 |

The zero limits of `CS` and `ES` agree with Gibbs's p. 8 statement; `AS` is
the nonzero Sprague reduction at zero slant.  These are numerical evaluations
of explicitly specified geometry, not interval certificates and not a new
universal-cover proof.
