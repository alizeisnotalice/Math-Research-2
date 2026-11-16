# R74: exact corner localization, compatibility and a Borel half-mass obstruction

2026-11-16, Asia/Shanghai. Threshold1, centered open axis-parallel squares, one common radius. Write H(rho) for the strict halo and D=H(sigma_n+nu)\H(nu) for actual deletion. This proof keeps finite atomic outside inputs separate from the explicitly named Borel example. Neither the general four-corner budget nor the unrestricted finite-atomic7/2 deletion bound is settled here.

## 1. Four disjoint effective source windows

Let n>=3, e=1/n², h=(n+1)/n², w=1−e, and let sigma_n give mass1/n² to each point of S_n={0,h,...,(n−1)h}². Assume finite positive outside nu keeps S_n a maximum eligible block of mass1. Equivalently, by the archived R64 theorem, every closed square Q of side R obeys

    (sigma_n+nu)(Q)<=max(1,R²).                        (1)

All strict full and residual witnesses have radius r<1/2: for R=2r>=1, both measures' masses are at most R². Put

    C=[−1/2,w+1/2]², J=[−e,1]²,
    Gamma_SW=[−1/2,−e]²,
    P_SW=(−1,1/2−e)².

The other corners Gamma and windows P are obtained by the two coordinate reflections s->w−s. Including the receiver boundaries changes no area. Every square centered at x in Gamma_SW with r<1/2 is strictly contained in P_SW, since its lower endpoints are greater than−1 and its upper endpoints are smaller than1/2−e. Thus source atoms exactly on the boundary of P_SW never participate in a valid corner witness. The two one-dimensional window intervals after reflection are

    (−1,1/2−e),   (1/2,2−e).

They have a separating gap of width e. The four two-dimensional P windows are pairwise disjoint.

Let sigma_P and nu_P be restrictions to one such P. Pointwise on its corner,

    D intersect Gamma
      =[H(sigma_P+nu_P)\H(nu_P)] intersect Gamma.      (2)

To check all radii in this identity, note that the restricted measures are dominated by sigma_n+nu, so they too have no strict witnesses at r>=1/2. For r<1/2 every witness square is contained in P, hence captures exactly the restricted masses. Both full and residual membership are therefore equal before taking their difference.

For feasibility, however, retain the entire internal measure: sigma_n+nu_P is dominated by sigma_n+nu and contains the eligible mass1 block S_n. It still satisfies all caps and has maximum eligible mass exactly1. The smaller sigma_P is only the captured internal measure in(2); it need not be a maximum block of the local input. In particular its mass is

    m_P=floor(n/2)²/n².                               (3)

Indeed the internal coordinates strictly below1/2−e are exactly those with indices0,...,floor(n/2)−1. For odd n the middle grid coordinate w/2 lies strictly above that cutoff. The local internal support remains in the northeast quadrant relative to its corner vertex0 (after reflection for other corners).

## 2. A genuine structural class that pays the complete7/2 budget

The separate proof `r074-single-effective-outside-corner-payment.md` proves the following all-scale statement, with no cap assumption: if sigma is any finite positive internal measure supported in the northeast quadrant of0, of total mass m, and the outside measure consists of one positive atom, then

    |[H(sigma+nu)\H(nu)] intersect {x<=0,y<=0}|<=m/2. (4)

The proof first moves all internal mass to0, which only increases full capture at such centers, and then estimates the exact two-atom full-minus-residual geometry. It is essential that the residual measure has one positive position. Caps need not survive the internal collapse because the two-atom bound itself needs none.

Combining(2)–(4), if each of the four effective windows P contains at most one positive outside position, with arbitrary additional outside sources elsewhere and with(1) still imposed globally, then

    sum_corners |D intersect Gamma|<=2floor(n/2)²/n²<=1/2.

The noncorner part of C has area

    |C|−4(1/2−e)²=3−3e².

Since D is contained in C, the complete actual deletion therefore satisfies

    Delta<=3−3/n⁴+2floor(n/2)²/n²<7/2.                 (5)

This closes a general structural subclass of the uniform-grid question for every n>=3. It is not a statement that a general outside measure can be replaced by one atom in each window. Multiple masses at the same location may first be combined; the condition concerns distinct positive positions.

## 3. Disjoint windows do not permit independent pasting

There is a very small finite counterexample to such a pasting step. Take n=5 and consider separately

    nu_SW=(8/25) delta_(-12/25,6/25),
    nu_NW=(8/25) delta_(-12/25,18/25).

Each point lies in its named effective window. Each individual outside measure, together with the FULL S_5, is an admissible unique-maximum, connected, nonterminal input.

For the first point, any mixed subset whose largest internal column index is j<=2 has mass at most(j+1)/5+8/25<=23/25<1. If j=3 or4, its horizontal span is at least(j+2)6/25, and

    (j+1)/5+8/25 <= [(j+2)6/25]².

Such a subset is ineligible. Subsets with no internal source have mass8/25; proper internal subsets have mass below1; S_5 is eligible with mass1. Hence S_5 is uniquely maximal. Its eligible pair with(0,6/25) connects the outside point to S_5. The second configuration is its reflection, with the same conclusions. Each has actual span36/25, total mass33/25, and nonterminal gap

    (36/25)²−(35/33)(33/25)=421/625>0.

But if both outside measures are pasted together, the closed unit square

    Q=[−12/25,13/25] x [0,1]

captures3×5 internal sources and both outside atoms. Its mass is31/25>1. The captured label even has bbox span24/25, so it is eligible with mass greater than1. Thus two individually admissible corner assignments violate global feasibility when combined. The windows' geometric disjointness is exact, while the square-cap constraints still couple their masses.

Mass-preserving D4 averaging is a different operation: because sigma_n is invariant, averaging rotated admissible outside measures does preserve the linear caps(1). There is nevertheless no justified Jensen principle or pointwise preservation of D under this averaging. The independent R74 review supplies a finite strict receiving-cell counterexample to pointwise preservation; it does not claim that the total corner objective is nonconcave or that averaging always lowers its area.

## 4. For outside densities<=1, unit caps imply all square caps

A useful general reduction applies to a different, explicitly Borel class. Let sigma be supported in a box whose two widths are at most1, and let nu have Lebesgue density between0 and1. If every closed unit square has(sigma+nu)-mass at most1, then every closed square Q of side R has mass at most max(1,R²).

For R>=1, if Q captures internal mass, choose a closed unit square Q_0 contained in Q that contains all the internal support captured by Q. This is possible coordinatewise: if[a,b] contains the captured old coordinates and[ell,ell+R] is the coordinate interval of Q, an admissible unit left endpoint lies in

    [max(ell,b−1), min(ell+R−1,a)],

which is nonempty because b−a<=1 and[a,b] lies in the larger interval. There is no internal mass in Q\Q_0. Thus

    (sigma+nu)(Q)
       <=(sigma+nu)(Q_0)+|Q\Q_0|<=1+(R²−1)=R².

If Q captures no internal mass, its density bound alone suffices. For R<=1, contain Q in any unit square to obtain mass<=1. This argument uses the outside density bound; it does not apply to arbitrary atomic nu.

## 5. An admissible single-corner Borel density

Set n=5 and, with P=P_SW and J as above, consider

    nu_c=c 1_(P\J) dx,   c=100/111.                   (6)

This is a finite Borel outside measure of density strictly below1, so H(nu_c) is empty at every radius. It is not a finite atomic outside input.

We prove every unit cap and then use Section4. Parameterize a closed unit square by its lower-left corner(u,v). If either coordinate is outside[−1,w], the square captures no internal mass and the outside density bound gives the cap. On[−1,w]², define

    F(u)=|[u,u+1] intersect [−1,23/50]|,
    G(u)=|[u,u+1] intersect [−1/25,23/50]|,
    K(u)=#{i=0,...,4: u<=6i/25<=u+1}.

The full unit-square mass is exactly

    K(u)K(v)/25+c[F(u)F(v)−G(u)G(v)].                  (7)

The complete list of one-dimensional breakpoints is

    −1, −19/25, −27/50, −13/25, −7/25, −1/25,
     0,   6/25,   23/50,   12/25,   18/25, 24/25.

On each open product cell, K is fixed in both coordinates and the outside area in(7) is bilinear. A bilinear function on a rectangle takes its maximum at a corner. The outside area is continuous at each endpoint, and the closed unit-square internal count at that endpoint is at least the neighboring open-cell count. Therefore checking the144 closed breakpoint vertices controls every product cell and every boundary stratum. No continuity of the atomic internal mass is assumed.

The accompanying Fraction certificate evaluates all144 values of(7), all at most1. The binding square has(u,v)=(−13/25,−13/25): it captures9 internal atoms, and its outside support area is

    (49/50)²−(1/2)²=444/625.

Hence its mass is9/25+(100/111)(444/625)=1. Any larger constant density on this same P\J violates this square. Thus c=100/111 is the exact largest density in this fixed Borel family. The full S_5 is retained in this cap check, including the middle row and column outside P; localizing sigma for this feasibility calculation would be wrong.

## 6. A strict single-corner area enclosure

For any x=(x_1,x_2) in Gamma_SW and r<1/2 capturing internal sources, the full witness square lies inside P and its intersection with J has area

    (x_1+r+e)(x_2+r+e).

The factors are positive once an internal source is captured. Thus a complete internal prefix of k×l sources yields strict full reception exactly when

    kl/25−4(1−c)r²−c(x_1+r+e)(x_2+r+e)>0.             (8)

All strict witnesses have r<1/2 by the already proved caps. Every captured internal prefix has1<=k,l<=2. Internal capture changes only at the complete source-distance events; between events the cost in(8) increases with radius. Closed capture with positive excess can be enlarged slightly to a genuine open witness.

Partition the closed corner square into64×64 equal cells of side23/3200. For a cell C_ij=[x_0,x_1]×[y_0,y_1], lower certification uses a radius r<1/2 at a source event for the lower-left receiver. Let k,l count the internal coordinates guaranteed to be captured throughout the cell. The strict sufficient margin is

    kl/25−4(1−c)r²−c(x_1+r+e)(y_1+r+e)>0.             (9)

The top-right receiver gives the largest J overlap, so(9) certifies every point of the closed cell. The square remains inside P uniformly. Increasing r by a sufficiently small positive amount produces an open common-radius witness for the entire cell, preserving the positive margin. The smallest saved margin is2209/284160000; a derivative bound of4 for the cost in(9), together with the finite positive distance of all saved r from1/2, makes this enlargement quantitative if desired.

For an upper exclusion, consider every possible prefix k,l and put

    r_min=max((k−1)h−x_1,(l−1)h−y_1).

Any witness capturing that prefix somewhere in the cell has r>=r_min. Prefixes with r_min>=1/2 are impossible. For the others a necessary positive upper excess is

    kl/25−4(1−c)r_min²
      −c(x_0+r_min+e)_+(y_0+r_min+e)_+ >0.            (10)

This follows by separately lowering the radius and the two receiver coordinates in the nondecreasing cost. If every such bound is nonpositive, the cell contains no full halo. All possible k,l are considered, so this is not an upper estimate based only on selected witness radii.

The final exact run certifies1655 closed cells by(9), excludes2161 cells by(10), and leaves280 additional cells undecided. Their interiors partition the whole corner. Consequently the ACTUAL Borel deletion, equal to the full halo here because H(nu_c) is empty, satisfies

    175099/2048000 <= |D intersect Gamma_SW|
                   <=204723/2048000.                  (11)

The lower bound is strictly greater than2/25=m_P/2, with excess11259/2048000. The upper bound is strictly below631/5000=1/8+3/(4·5⁴). Thus(6) is a strict Borel counterexample to the auxiliary half-local-mass claim, but is not a counterexample to the main single-corner sufficient threshold. The strict finite-atomic half-mass question with multiple outside positions remains open. R72/R73 show why a Borel deletion lower bound cannot simply be carried to finite atoms.

## 7. Actual execution and review scope

An initial32×32 exact midpoint membership probe used the single-window candidate densities100/111,567/605,3757/3897 for n=5,9,17. Its area sums were only midpoint estimates; only n=5 was then certified. The all-cap and area certificate uses Python3.15.0 standard-library Fraction. The final run checks144 cap vertices, executes4861 lower-event tests and3733 upper-prefix tests on the64×64 cells, with exact stored witnesses and exclusion bounds. It exited0 and its PID84403 was confirmed absent. The earlier lower-only run PID77691 also exited0; its144 cap vertices and1655 lower-cell records match the final run exactly. Its source snapshot and receipt retain that superseded execution without adding a second mathematical claim.

The independent reviewer additionally reconstructs all529 point/open-interval unit-origin strata and1156 corner limits, checks the1655 lower cells from the25 literal old atoms, and checks all2161 excluded cells with every feasible internal prefix. Its3283 upper tests differ from the author's short-circuit count. The same rational enclosure(11) and both threshold comparisons were independently recovered. Neither run enumerates a finite outside approximation or transfers the Borel area to one.

The main positive R74 conclusion is the genuine one-effective-outside-position payment(5). The general four-corner bound with arbitrary finite outside multiplicity, and the resulting uniform-grid7/2 target, remain unresolved.
