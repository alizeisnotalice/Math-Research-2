# r024 log-scaled distance aggregator

This is an independently derived aggregator with source hash `d112b9d9…`.
It adds a common \(\max\log Z\) shift, log-normalizer output, and raw
root-scale tuple guards. It does **not** repair a claimed r023 underflow bug:
the frozen r023 Decimal aggregator already reconstructs weights from `logZ`
when `Double Z` is zero. The new test explicitly verifies that both
aggregators preserve ratios and delta diagnostics after a global `logZ-1000`
shift.

Run `python3 -B test_logscaled_sut.py` from this directory for the retained
real-CLI hand fixtures. The frozen r023 aggregator and terminal-source copy
are included solely for this comparison; neither changes the terminal kernel.

The retained n=16 and remote-n12 summaries are finite outputs. The n=32
range record documents why incomplete random-tau data are outside this
terminal-distance protocol; it is not an estimate.
