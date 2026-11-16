# R75: certified all-parameter orthogonal corner payment for S3

2026-11-16. This theorem uses the unchanged analytic interface `r075-orthogonal-two-outside-full-grid-caps.md` and the finite exact parameter certificate cited below. It concerns actual full-minus-residual deletion, not a full-halo surrogate.

## 1. Continuous two-parameter theorem

Let sigma_3 give mass1/9 to each point of {0,4/9,8/9}², let

    nu=b delta_(-a,0)+b delta_(0,-a),
    0<a<=1, 0<=b<=B_3(a),
    Gamma=[−1/2,−1/9]²,
    D=H(sigma_3+nu)\H(nu).

Here B_3 is the exact all-square capacity function in the frozen interface. Then

    |D intersect Gamma| < 1/18.                       (1)

The interval0<a<1/9 forces b=0, and the exact corner area there is1/324. For the remaining domain, the parameter-tree certificate starts with the entire closed rectangle

    [1/9,1] x [0,4/9].

At a node [aL,aU]x[bL,bU], monotonicity of B_3 discards only the infeasible strip b>B_3(aU); its boundary b=B_3(aU) is retained. Every split gives two closed children whose union is the retained parent box. Over each retained parameter box, Section4 of the analytic interface encloses the four possible full contributing labels from outside and encloses the three true residual labels from inside. This gives an upper bound on the actual deletion uniformly in every parameter in that box.

The saved certificate has71 nodes,35 splits,36 certified terminal boxes,7 removed cap strips and maximum depth9. There is no pending node or uncovered parameter. Its72-bit dyadic square-root bounds are checked by exact integer squaring. All rectangle union-difference areas, cap cuts, parent/child coverage tests and comparisons with1/18 are rational. Floating heap priority only chooses processing order.

The largest terminal uniform upper bound is exactly

    M=1563026333441017474428935886295754233898431
       /28224380641890319913506143438820582531006464,

and

    1/18−M=4994813330666965210294304749833684490817
       /28224380641890319913506143438820582531006464>0. (2)

Every feasible parameter belongs to a certified terminal box, so(2) proves(1). This establishes a strict payment for the complete continuous family, including a and b boundary values. It does not identify the family's optimal area or infer the theorem from sampled points.

The certificate is `computation/results/round075-orthogonal-boxes/parameter-tree.json.gz`, SHA256
`36bd0b4bad73c7047601a66f96c6cf76a3ac77e1794ecdcee7c1b4e819cfa99c`.
The producer summary and source record the successful single construction run. The independent replay reconstructs all71 node areas with horizontal slices, checks284 square-root bounds and all coverage/cap relations, and executes572 rational slice calculations. It recovers the same exact M. This is a full certificate replay, not a second search for a partition.

## 2. A genuine four-corner structural class

Now let nu be any finite positive outside measure such that the COMPLETE sigma_3+nu obeys every square capacity and S3 has maximum eligible mass1. In each of the four disjoint effective corner windows P from R74, assume the restricted outside measure is either

- supported at at most one positive position, with arbitrary position and mass; or
- a D4 image, relative to that internal grid corner, of an equal-mass orthogonal pair b delta_(-a,0)+b delta_(0,-a), with0<a<1.

The parameters may differ between corners. Outside the union of the four P windows, arbitrary additional positive outside sources are allowed. These assumptions describe one already feasible common input; they do not authorize pasting separately feasible corner constructions.

Restricting nu to one P while retaining the full sigma_3 only decreases every square mass. Thus that restricted complete-grid input still satisfies all capacities. In the pair case the exact capacity theorem gives b<=B_3(a), and R74 all-radius localization makes(1) the correct actual corner deletion bound. In the single-position case the R74 uncapped single-effective-outside theorem gives the non-strict bound1/18. Hence

    sum_corners |D intersect Gamma_corner|<=2/9.

All deletion lies in C=[−1/2,25/18]². Its noncorner portion has area

    |C|−4(1/2−1/9)²=3−3/3⁴=80/27.

Paying that portion by its full area yields the complete genuine deletion bound

    Delta<=80/27+2/9=86/27<7/2.                       (3)

This enlarges the R74 structural class for n=3 to permit an orthogonal equal-mass pair in any subset of the corners, while retaining full global feasibility and the complete residual halo. It is not an arbitrary two-outside-position theorem, and the unequal-mass capacity polytope alone does not prove the same area bound for unequal pairs.

For the pure family in Section1, only the southwest P contains outside mass; the other three corners each have exact deletion1/324 from their single internal corner atom. The sharper coarse full-area consequence is

    Delta<=80/27+1/18+3/324=109/36<7/2.              (4)

At a=1 the outside atoms are on excluded P boundaries, so the same conclusion holds. If one extends the pure orthogonal family to a>1 while separately retaining full global capacities, none of its outside atoms lies in any effective P; all four corners are pure internal and(4) remains valid. No all-a>1 mass-polytope formula or area search is needed for that observation.

## 3. Scope and evidence

The prior961-point screen used exact global capacity checks but floating local areas; it suggested that the specified parameter family was worth certifying. It does not prove(1). The complete finite covering tree and its independent exact replay prove(1) through the analytic envelope theorem.

The no-cap three-atom example from the companion R75 proof has southwest-quadrant deletion greater than half the internal mass. It violates the complete-grid capacity condition after the proposed S3 embedding, and its whole quadrant differs from Gamma. It therefore does not contradict(1)–(4).

The all-n exact capacity characterization and its unequal-mass extension are analytic results. The area certificate here is restricted to n=3 and equal orthogonal pairs. General finite outside measures, arbitrary two positions, unequal pair areas and the sharp C2 remain unresolved.
