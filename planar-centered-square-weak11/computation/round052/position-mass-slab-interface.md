# R52 position/mass candidate and omission certificate

Scope: a local two-parameter candidate-bbox certificate, conditional on the pinned accepted R39 one-axis source catalogue, coordinate order and geometric prefilter stability. This stage does not certify endpoint order or the area union formula.

`build_position_mass_slab(center_data, mass_velocity, x_velocity, source_slab, T="1/1000", t0="1/2")` returns `position-mass-slab-v1`. The center data has masses at t0, x coordinates before applying the source slab s0. Each original atom has a strictly positive center mass w and an affine mass w+mu*dt. The source velocities are 0/1 and copy geometry is inherited without altering the frozen enumerator.

For the returned rational rho, |dz| <= rho H and |dt| <= rho T. Each candidate raw endpoint `[a,v,sign,q0,q1]` means

    a + v*dz + sign*sqrt(q0 + q1*dt).

Here q0=M0/(4 lambda), q1=MU/(4 lambda). The nonnegative square root is used. Each included radius branch stays fixed, and candidate radicands stay strictly positive. Rational-square center values with q1 nonzero remain radicals. Root/cap equality with q1 zero is persistent and is represented by the cap branch.

All eligible bboxes are scanned with two integer prefix sums, giving exact M0 and MU. Geometrically empty bboxes with M0=0 necessarily have MU=0 because center atom masses are strictly positive. All other rejected bboxes retain a sufficient failed condition: x/y span or one of four domain reach conditions. Write that condition d >= r, where d >= 0 is affine in dz and r^2=min(R^2,(M0+MU*dt)/(4 lambda)). At the center q=d0^2-r0^2. For |dz| <= rho H and |dt| <= rho T,

    |(d^2-r^2) - q| <= rho L,
    L=|d_velocity|*(2|d0|+|d_velocity|H)*H + |MU|T/(4 lambda).

The bound for r^2 follows from the 1-Lipschitz map u -> min(R^2,u), so omitted boxes do not require their own fixed cap branch. Taking rho <= q/(2L) gives d^2-r^2 >= q/2 > 0. Clip reasons additionally require rho <= d0/(2|d_velocity|H), keeping d >= d0/2. All comparisons use scaled exact integers. The global omission radius is the minimum, across all omitted bboxes, of a sufficient chosen reason; within one bbox the reason with the largest certified radius is selected.

When q=0, the failed strict condition can be omitted only if d_velocity=0 and either MU=0 or the mass radius stays capped on the full requested T interval. Otherwise that bbox remains in the candidate family, tagged `empty_equality_wall`, even though its center area vanishes. It may become active as mass changes.

Candidate root/cap branch separation is shrunk by half its affine squared-radius gap, and each atom mass is kept at least half its center value. The full requested T interval also has a strict position-independent total-mass R-tail bound, checked independently of the local rho. The final rho is the minimum of omission, cap-branch, atom-positivity and requested bounds.

`box_records` correspond one-to-one and in order to `boxes`/`candidate_boxes`. The `bbox` keys refer to the pinned x/y source coordinate classes. `mu`, `q0`, `q1`, `radius_branch`, `q_full_domain_lower/upper` and `cap_gap_lower` are exact rational strings. `center_active_keys_match_source_slab` verifies all strictly active center bboxes match the prior one-axis scan. Additional center equality candidates are separately listed. The counts ledger accounts for every eligible bbox.

A support-aware endpoint sweep must now establish each relevant raw endpoint order and clipping relation throughout a positive subrectangle. Equality at the center between different five-component expressions is not presumed persistent. Only after that step may the union area and its position slope/mass derivative be certified. Neither this stage nor sampled slices establish a stable cell covering all t in [0,1].
