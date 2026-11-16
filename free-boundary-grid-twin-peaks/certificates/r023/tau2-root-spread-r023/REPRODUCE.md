# Reproduction

```sh
python3 tau2_pair_kernel.py > tau2-pair-kernel.json
```

This uses the exact independent-set inclusion--exclusion formula for the
top-two completion kernel.  It enumerates independent sets, not labellings;
all sums and divisions use Python integers for \(n=3,4,5\).
