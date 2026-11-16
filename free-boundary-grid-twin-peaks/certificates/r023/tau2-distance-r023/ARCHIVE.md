# r023 tau=2 terminal distance

Contains final terminal source, small n=3/n=8/n=32 receipts and raw outputs,
the strict Fraction aggregator, and the root-audit receipt. Large traces and
binaries are excluded.

The portable high-value test is

```sh
python3 -B root-audit/test_hand_fixtures_relative.py
```

It submits an independently specified unequal two-replicate rational fixture
to the archived aggregator, checks means and variance squares for n=3/n=7,
checks both strict-threshold sides at `1e-90`, and confirms malformed bin/L²
inputs are rejected. It is the source of `root-audit/hand-fixture-receipt.json`.

`aggregate/test_r023_relative.py` and `aggregate/test_tau2_distance.py` are
useful legacy regression checks of retained n=3 records, strict inequality and
underflow handling. They do **not** inject the hand fixture into the SUT, so
they are not a substitute for the root-audit test above. The absolute-path
legacy scripts are retained only as historical provenance and are not the
portable archive entry points.
