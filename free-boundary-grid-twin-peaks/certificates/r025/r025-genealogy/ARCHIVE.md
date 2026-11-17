# Cross-replicate normalizer diagnostic

This read-only diagnostic compares the completed local `strat` and MacBook
`multi` n=24 records separately.  It reconstructs weights with Decimal
arithmetic from `logZ` after a common shift; it does not infer unrecorded
genealogies.

In both batches, concentration is already present in `Z`, and multiplying
by the recorded `L` factor modestly *increases* ESS.  In the local batch,
Z ESS is 1.138 with maximum share 0.9352, whereas ZL ESS is 1.450 with
maximum share 0.8080.  Thus these records do not support attributing the
concentration to the L-lift.

The two proposal methods are reported separately in this study.  We make no
combined higher-precision claim because both finite batches have severe or
material normalizer concentration and neither supplies a mixing certificate;
this is a methodological limitation of the present evidence, not a theorem
forbidding pooling across methods.  Delta-SE remains a diagnostic, not a
confidence interval.  Parent IDs, offspring counts, per-rank ESS, and
ancestor hashes were not written by the terminal kernel, so genealogy cannot
be recovered retrospectively.

In addition, the sibling normalizer--distance audit shows that the two
complete-law plug-in ratios fail an exact structural consistency check.  This
is stronger than the ESS warning for interpretation: no trend claim should
use their Q/tail values as if they were exact target-law quantities.
