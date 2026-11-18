# Combining EH depth 0.00055 with the AH cap of Round 105

Keep the Gibbs angle \(\sigma=1.5497^\circ\) and Round 105's AH halfplane, defined by the hinge \(W=Z+0.38(T_a-Z)\) and the counterclockwise turn \(0.00033\). The EH line joins \(T_{\rm EH}\) to \(U_1=C_3+R_{-0.00055}(U-C_3)\). The intersection of the single-EH cover with the AH halfplane is a convex universal cover.

First apply the complete single-EH theorem at depth \(0.00055\) to the input. It returns a constant-width-one copy \(K\subset P\) avoiding the old AH and the deeper EH interiors. If it also avoids the new AH interior, retain this copy. Otherwise the new AH hit satisfies the Round 105 chart and trigger estimates, which depend on the same \(\sigma\) and AH parameters. Reflect by \(Q_A\). The reflected copy lies in \(P\) by the AH trigger. Since the input copy avoided the positive \(C\) cap, its reflection avoids \(F'\). An interior hit of the deeper EH would imply a hit of \(F'\), so this reflection also avoids EH.

If both copies hit the new AH interior, apply Round 105's five-point double-hit relocation. The local contact errors and the two angle branches are unchanged at the same AH data. Repeat all support comparisons against the deeper EH line \(U_1T_{\rm EH}\) together with the same two hat supports, old AH inequalities and new AH line. The directly tested EH row is `U1-Te-0`. Across 1340 parameter cells and 100 angle cells per bracket, the minimum perturbed compatibility-sign slack exceeds \(8.08818\times10^{-5}\) and the minimum placement slack after deducting the actual support error exceeds \(0.0002064112\). The five protected Sprague arc conditions hold at the same pose by the exact contact argument.

This composition applies after all branches of the deeper EH theorem. It uses the theorem's guaranteed output in the old AH cover, followed by the AH-only classification above. Avoidance of ear interiors gives placement in the closed candidate by the earlier density argument. Constant-width completion then covers all inputs of diameter at most one.

The AH cap is still the same region of the original Gibbs body. Its height is above \(0.518\), while the deeper EH cap is below the EH center height \(s_y<0\). Thus the two area cuts are disjoint, and the Round 106 AH Green integral and the deeper-EH Green integral can be subtracted from the original fixed-angle area, each using outward interval endpoints.

## Replay

From `computation`, run the standalone EH replay and the AH combination with Python 3.13.18 and python-flint 0.9.0

```sh
python -B certify_eh_joint_switch.py 1.5497 0.00055 36 "$OUT" \
  --eh-relocation --joint-eh --old-area eh-joint-1p5497-old-area.json \
  --replay eh-cap-switch-joint-550.json
python -B certify_deep_ah_eh.py --combined-deep
```

The output is `ah-eh-combined-550.json`. The triangle bound in this receipt also certifies an area upper endpoint directly. The sharpened area comparison uses the full, unchanged AH cap integral from Round 106.
