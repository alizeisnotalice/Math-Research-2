# R67 positive padding and mass-sensitive erosion states

One prescribed 10-source input tests the proposed replacement of the old bbox
container by a union of uniform-radius source neighborhoods. All nine grid
masses are positive, the grid spacing is smaller than sqrt(theta L), and the
neighborhood union is exactly the old container. Complete directed geometry
still gives a budget lower bound about 3.50308723161359, while the true deletion
upper bound is about 0.999688823837595. This excludes that support-only repair.

The executed adapter input and the author's proof input have different metadata:

* executed adapter SHA: `26ae7ed968b11ab932c52b513bdcc2569647a6345d76e72fc4b782b2f3b4c164`;
* author input SHA: `68cfb1045685d6c51de1a3d8d6f320154454c34e8021012b56ce81f15700202b`.

Their exact position/mass triples and S agree. The executed builder enumerated
all 1023 nonempty source subsets and deduplicated their bbox closures; its
46-closure/18-eligible catalog agrees with the frozen coordinate-interval
algorithm. It also checked both directed uniform-neighborhood unions against
the container by exact zero-area differences. The author file was not silently
substituted for the actual executed input.

`whole_erosion_state.py` constructs two new costs from every strictly eligible
full-source bbox closure T intersecting S. Write b=m(T\S).

* W uses E_T(m(T)) minus the receiver rectangle at mass b on the **full T bbox**.
* Wprime uses E_T(m(T)) minus the receiver rectangle at mass b on the **actual
  tight bbox of T\S**. A zero b gives an empty kernel.

All label rings are decomposed and then aggregated as a true union. Neither
state uses theta or a complete-remaining-halo discount. The remaining halo is
computed separately only to enclose true deletion and overpayment. Directed
inclusion checks verify the computed compatibility of
Delta subset Wprime subset W subset H(full). Individual a=m(T intersect S)
caps are not summed as a common budget.

Exactly two state inputs were authorized and executed:

| Input | Eligible labels intersecting S | W/Wprime upper display | Deletion upper display | Overpayment upper display |
| --- | ---: | ---: | ---: | ---: |
| Positive 9+1 padding | 17 | 1.0006653863375947 | 0.9996888238375947 | 0.0009765625 |
| R64 final 21-source winner | 54 | 1.2992259316209571 | 1.2359935097459571 | 0.063232421875 |

Both have L=1. W and Wprime have identical computed interval endpoints on
these two inputs; this is not a claim of general equality. Full halo area,
remaining halo area, label masses/spans, both state intervals, refinement and
overpayment intervals are retained in the actual JSON outputs. The original
wrappers preserve absolute execution commands and PID receipts. All three
mathematical child processes (one padding budget, two state jobs) exited zero
and were checked absent. No additional inputs, LP, theta scan or MATLAB run was
performed. Independent geometry review is a separate artifact.

To reproduce from the repository root, using the archived dependencies and a
fresh data-volume output directory:

```sh
python3 -B planar-centered-square-weak11/computation/round067/run_positive_grid_padding_control.py \
  --budget planar-centered-square-weak11/computation/round066/theta_competitor_budget.py \
  --module planar-centered-square-weak11/computation/round064/mass_polytope_variation.py \
  --out-dir /Volumes/data/CodexTemp/T2/r067-replay/padding-build

python3 -B planar-centered-square-weak11/computation/round067/whole_erosion_state.py \
  --module planar-centered-square-weak11/computation/round064/mass_polytope_variation.py \
  --erosion-source planar-centered-square-weak11/computation/round066/erosion_ring_diagnostic.py \
  --input planar-centered-square-weak11/computation/inputs/r067-positive-grid-padding-actual-control.json \
  --out /Volumes/data/CodexTemp/T2/r067-replay/padding-state.json

python3 -B planar-centered-square-weak11/computation/round067/whole_erosion_state.py \
  --module planar-centered-square-weak11/computation/round064/mass_polytope_variation.py \
  --erosion-source planar-centered-square-weak11/computation/round066/erosion_ring_diagnostic.py \
  --input planar-centered-square-weak11/computation/inputs/r066-r064-final-input.json \
  --out /Volumes/data/CodexTemp/T2/r067-replay/r64-final-state.json
```

The actual runs used Python 3.13.18 with standard-library integers/Fraction and
72-bit dyadic endpoint enclosures. Source/hash/domain bindings are in
`domain-and-execution-binding.json`. These input-specific certificates do not
prove the new general ring-selection conjecture or change the C2 bounds.
