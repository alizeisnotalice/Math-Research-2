# R57: a certified kink in the full periodic halo area

Let A(t,z) be the area of the complete R=10 halo union in the given fundamental domain for the R43-to-R48 affine mass family and the two-label x motion. The exact physical position is z0=P8/10^6. The total source mass M is constant and the level is lambda=1. The strict full-scale tail bound is pinned in the accepted endpoint inputs, so the finite-R halo union equals the full-scale halo union for this family.

The source-pinned R52 endpoint-product lemma supplies the complete common 26,327-bbox family throughout t in [0,1] at z0. The R57 prerequisite verifier checks the identical endpoint active keys, source geometry and coordinate order; verifies every endpoint expression against the midpoint affine-mass representation; and checks every candidate squared radius remains strictly between 0 and R^2 throughout [0,1]. It also pins the endpoint omission/source radii and strict tails. This uses the prior complete enumerations; it does not rerun their 45,633,320-box scans.

There is a unique crossing t* at z0 defined by

    sqrt(qA) + sqrt(qBm + qBprime*(t* - 1/2)) = C,

with all exact rational coefficients in the paired certificate. qA>0 is constant, qBprime<0, and C-sqrt(qA)>0. Thus the equality has the explicitly isolated root recorded there, approximately 0.43720851363210126589898486932909244. The endpoint A is constant and endpoint B strictly decreases in t.

For each side choose the exact rational midpoint of the wall enclosure plus or minus 10^-12. At that rational t every distinct normalized x/y endpoint is strictly ordered by directed radical intervals. Full endpoint identity retains the mass derivative q1; only a constant rational-square radical may be folded into an offset. There are 25,232 distinct x endpoints and 10,821 y endpoints. Over each closed segment from the side parameter to the exact wall, all consecutive endpoint relations remain strictly separated by directed interval bounds except the single known A/B pair; no y relation is ambiguous. The exact crossing equation and strict monotonicity prove its one-sided order. Consequently the full side orders extend to the corresponding one-sided limit. No unexamined simultaneous event is omitted.

Clip all rectangles by domain ranks. For each x gap, a dynamic coverage tree computes the exact rank formula for the y union length Y and its t derivative Y_t. A covered node contributes its upper-minus-lower endpoint and the corresponding derivative difference; an uncovered node contributes the sum of its children. Hence on either side

    A = sum gap*Y,
    A_t = sum (gap_t*Y + gap*Y_t),
    beta=A_z = sum (v_right-v_left)*Y.

The zero-width limiting gap at the collision still contributes gap_t*Y. Thus the calculation includes all events and possible cancellation in the complete area. Endpoint values and derivatives are evaluated at an outward enclosure of the exact root, not at nearby floating-point samples. All bounds use exact dyadic integers, with products at doubled precision.

The final 256-bit paired bounds imply strictly

    A_t(t*+,z0) - A_t(t*-,z0) < 0,
    beta(t*+,z0) - beta(t*-,z0) > 0.

Decimal displays are approximately -5.1936017367457632415007786512377074e-8 and +0.069460014800942633913127104246681787 respectively. The exact rational enclosing endpoints are in full-area-pair-certificate.json. Both one-sided A_t remain positive. The two area-limit enclosures overlap, as continuity requires. Since J=lambda*A/M and M_t=0, the normalized derivative jump is strictly negative, approximately -2.4612669543873546075851796487442997e-9.

Validation: both production sides exit 0; final 256-bit enclosures lie inside their initial 192-bit counterparts. A separate 5,000-update direct interval-merge comparison checks both tree length and derivative (seed 570021). The source prerequisite checker exits 0 after all 26,327 comparisons. Its first invocation found a schema-field mismatch (tail lambda lives under tail, not the slab root); after correction the recorded final run passed. No area-sweep arithmetic was changed by that correction. Independent four-state wall-jump review is maintained by the lower-audit worker.

This proves a nonsmooth point of the complete family area and excludes a single analytic raw-order cell crossing this wall. It does not prove local optimality (both mass-side slopes are positive), change the weak-type lower bound, or determine the exact weak-type constant.
