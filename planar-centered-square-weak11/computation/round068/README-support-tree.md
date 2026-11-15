# R68 bounded support-tree/chain supplement

From the topic directory, replay the three fixed rational fixtures using an
existing Python with SymPy (actual run: Python3.13.18,SymPy1.14.0):

```sh
python3 -B computation/round068/check_support_tree_payment.py \
  --input computation/inputs/r068-support-tree-payment.json \
  --out /Volumes/data/CodexTemp/T2/r068/support-tree-replay.json
```

The script enumerates31,31,511 subsets, checks positivity/unique maximum/
connectedness/nonterminal/thick-axis domain, uses exact Kruskal weights,
and integrates the support-square union by independent rational cell
compression. Its two SymPy identities check the edge-area and telescoping
algebra. No randomseed is used; all arithmetic outputs are exact Fractions.

The general proof is `proofs/r068-support-tree-payment-and-chain-union.md`.
This code does not compute W or Wprime, perform a floating LP, launch MATLAB
or certify a general3.5 selection proposition. The two positive cases are
paid by the theorem; the two-row case only obstructs making the support-
only/tree sufficient test universal. Saved receipt identifies the final
actual output and confirms the launched Python process exited.
