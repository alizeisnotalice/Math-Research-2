# Round 005 matching computation record

## Scope

The seed is \(a_{(u,v)}=(-1)^{u\cdot v}/8\), \(u,v\in\mathbb F_2^3\).
For every nonzero \(\delta\), the program forms the 32 unordered pairs
\(\{x,x+\delta\}\), separates their products into 16 positive and 16
negative pairs, and joins a positive pair to a negative pair.  A join emits
its four-vertex incidence row over \(\mathbb F_2\).

The structured start labels a pair by the difference of its evaluations on

\[
 V=\operatorname{span}\{1_{v=i}:0\leq i<8\}
 \mathbin{+}\operatorname{span}\{u_j1_{v\bmod2^j=p}:0\leq j<3,\ 0\leq p<2^j\}.
\]

Positive and negative pairs having the same label are matched.  This is a
15-dimensional linear space.  The recorded matching has rank \(R=49\), so
the common sign-code has linear dimension 15 and its graph quotient has
affine dimension **14**.  Adding an independent character gives Fourier span
dimension 15; these three quantities are deliberately not conflated.
The original bent seed has full Fourier support \(\mathbb F_2^6\), hence
bare Fourier-span rank 6 before the graph lift.

## Actual runs

Apple clang 21.0.0 and Python 3.15.0 were used locally.

| seed/profile | moves | elapsed | initial/best `R` | result |
| --- | ---: | ---: | ---: | --- |
| `2026111711`, base | 500,000 | 58.5372 s | 49 / 49 | no strict improvement |
| `2026111712`, nonlinear Maiorana--McFarland | 500,000 | 59.4449 s | 57 / 57 | no strict improvement |

Moves were 70% two-pair swaps, 25% three-pair cycles, and 5% whole-delta
rematchings, accepted by a linearly cooled rank annealing rule.  These are
finite searches only; their failure to improve is not an upper-bound proof.
The MM comparison seed used the linear Gray-type permutation
`(0,1,3,2,6,7,5,4)` and a quadratic `h(v)`; it is not presented as a cubic
seed.  Its canonical matching does not preserve the specified `V15` labels.

## Certificate verification

[`structured-v15-certificate.json`](structured-v15-certificate.json) stores
all 1008 incidence rows, all 63 matching bijections, and a 15-vector kernel
basis.  The independent verifier checked:

- for every \(\delta\), exactly 16 positive and 16 negative pairs, each used once;
- `rank = 49`, kernel-basis rank 15, and annihilation of every incidence row;
- every one of the \(2^{15}=32768\) kernel combinations by exact integer
  64-point Walsh transform (all values are \(\pm8\));
- the seven-dimensional constant-plus-coordinate character subspace \(C\)
  lies in the kernel, and its 64 frequency-evaluation columns are distinct.

The cancellation mechanism itself proves Booleanity for all kernel elements:
each matched positive/negative pair has equal lifted pair parity, so its two
autocorrelation contributions cancel independently for every nonzero
\(\delta\).  The full enumeration is an additional finite cross-check.

Reproduce from this directory (put binaries and transient outputs under
`/Volumes/data/CodexTemp/T2/a28/r005`):

```sh
CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/clang clang++ -O3 -std=c++20 matching_search.cpp -o /Volumes/data/CodexTemp/T2/a28/r005/matching_search
/Volumes/data/CodexTemp/T2/a28/r005/matching_search 2026111711 500000 1.0 /Volumes/data/CodexTemp/T2/a28/r005/search-base.json base
python3 verify_matching_result.py /Volumes/data/CodexTemp/T2/a28/r005/search-base.json --out /Volumes/data/CodexTemp/T2/a28/r005/search-base.verify.json
```

SHA-256:

- `matching_search.cpp`: `ccc16fdfb8dfa44934f18a3e29772a05c7a81247c4127db099e27b1741378529`
- `verify_matching_result.py`: `342563a86407551a09aa2786e74dd77712a20e3b7afce3e7b13d1f0c86c4881b`
- structured certificate: `de16215621481fdc38ec76c29330ccd2a37f70eb9ab392db0df6e4f7d320a0b0`

All recorded local and remote search executions used the immediately prior
source hash `0fbac6b55e586696d2287f18b75d779c26b2074c633633afd32147edf951428e`.
The current source hash above differs only because the MM comparison comment
was corrected from “non-affine/non-linear” to “linear Gray-type/quadratic”; no
executable code, search move, or certificate serialization changed.

## Remote attempt

`macbook-math` briefly responded to `hostname`/`uptime`, but reported load
averages about 128.  Two subsequent low-impact SCP attempts (20 s and 30 s
banner timeouts) exited 255 before transfer; no remote source, process, or
certificate was created.  No other process or SSH service was modified.

After those transient failures, a single SSH stdin bootstrap did complete one
low-priority remote base search in
`/Users/zhengzhihao/CodexRemoteMath/a28-selector/r005/remote-search` without
SCP.  With seed `2026111713`, 500,000 moves took 26.9651 seconds and retained
`R=49`; the remote fast verifier independently checked all 32,768 masks in
0.788007791 seconds.  Its returned certificate and verifier record are
[`remote-search-2026111713-certificate.json`](remote-search-2026111713-certificate.json)
and [`remote-search-2026111713.verify.json`](remote-search-2026111713.verify.json),
with SHA-256 `ac1ecb7ab9708fdea549cab7d014116bcd1ffc2c2e2ed5ad10bdb076bc64db5a`
and `3b79246c0154a35ae1a271b448379d16791949e10eefe42d9f7b64ca52ef8291`.
The latter is the SHA-256 of the returned **verification receipt JSON**, not
the verifier source.  The remote receipt reports fast-verifier source SHA-256
`918de384903923aa2de96693d8f76cd7e089d2cdf0289be43eeb2aaad468c4aa`.

The archived remote bootstrap is a corrected reproduction command; the first actual validation invocation omitted --out, and a subsequent validation of the saved result succeeded. The search itself was not repeated. The MM result and receipt are archived as search-mm-2026111712.json and search-mm-2026111712.verify.json.
