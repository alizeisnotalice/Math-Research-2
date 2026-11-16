# Exact helper timing averages, r022

For a fixed K=3 peak set and an unoriented domino cover whose heavy-side
centre is the helper \(a\), the program counts all words with exactly that
peak set.  If its two same-side peaks are \(u,v\), it accumulates
\[
g=|\operatorname{pos}(u)-\operatorname{pos}(v)|,\qquad
h=N-\max(\operatorname{pos}(u),\operatorname{pos}(v))-1,
\]
and reports the exact sum of \(g/h\).

The backward subset DP \(H(S)\) counts valid completions after prefix \(S\);
the forward DP stores the prefix count and the first moments of the two
peak positions.  At the transition which inserts the later one at rank
\(k\), it adds \((k\,C-\sum\operatorname{pos})H(S\cup\{w\})/(N-k-1)\).
Thus every word is counted exactly once.  A vertex is admissible precisely
when it has no earlier neighbour iff it belongs to the fixed peak set.

Actual exact results:

| n | (P, unoriented-cover) records | covered-word weight | sum g/h | mean g/h |
|---|---:|---:|---:|---:|
|3|12|129024|51072|19/48|
|4|72|1194718147840|975017025920/3|217637729/800034474|

An earlier implementation mistakenly retained only the case in which the
smaller numbered endpoint was the heavy-side helper; its half-cover outputs
are retained under buggy-n*-half-covers.json and are not certificates.
The corrected n=3 result exactly matches r021's 129024 and 51072 totals.
An independent 9! enumerator constructs every unordered grid edge and
tests both its 2+1 and 1+2 heavy-side cases; it reproduces all twelve
individual records.  Thus the agreement is now a same-domain check, not
an inferred direction factor.

The n=5 single-DP memory layout consists of four U128 arrays over
\(2^{25}\) masks, exactly 2 GiB before allocator and output overhead,
within but close to the 2.5 GiB cap.  It was not launched in this bounded
round because the all-cover executable would repeat that allocation for
many peak pairs; it needs a fixed-record mode before a safe one-record run.
