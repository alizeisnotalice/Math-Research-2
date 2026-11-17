# r031 finite certificate replay

This is a finite-certification assembly, not a proof of a limiting statement.
All commands below write only to an external, initially empty output directory.
They require Python 3; the optional A source smoke additionally requires a
C++20 compiler.  Run from the directory containing this file.

```sh
set -eu
R31="$PWD"
OUT="${OUT:?set OUT to a writable directory outside r031}"
CERTIFICATES="${CERTIFICATES:?set this to the certificate root containing r029 and r030}"
mkdir -p "$OUT"
test "$OUT" != "$R31"

(cd "$R31/r031-finite-compute-public" && shasum -a 256 -c SHA256SUMS)
(cd "$R31/r031-n7-orbit-public" && shasum -a 256 -c SHA256SUMS)
```

## A. Saved h1--h5 finite checkpoints

This reads the saved n=3,4,5 residue checkpoints and the independent complete
r029/r030 L-histogram references. It does not launch n=7.

```sh
python3 -B "$R31/r031-finite-compute-public/short-gap-h45/source/verify_h1_h5_portable.py" \
  --source "$R31/r031-finite-compute-public/short-gap-h45/source/short_gap_h45_mod.cpp" \
  --checkpoints "$R31/r031-finite-compute-public/short-gap-h45/raw/checkpoints" \
  --r029-tail "$CERTIFICATES/r029/r029-tail-full-n5" \
  --r030-geometry "$CERTIFICATES/r030/r030-gap-geometry" \
  --output "$OUT/h1-h5-finite-crosscheck.json"
```

An optional finite source smoke is:

```sh
env CLANG_MODULE_CACHE_PATH="${CLANG_MODULE_CACHE_PATH:-/Volumes/data/CodexTemp/cache/clang}" \
  clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pthread \
  "$R31/r031-finite-compute-public/short-gap-h45/source/short_gap_h45_mod.cpp" \
  -o "$OUT/short_gap_h45_mod"
"$OUT/short_gap_h45_mod" 3 4 0 1 "$OUT/n3-v4-p0.json"
```

## B. n=32 merge-RB raw audit

```sh
python3 -B "$R31/r031-finite-compute-public/merge-rb-n32/source/verify_n32_public.py" \
  --raw-gz "$R31/r031-finite-compute-public/merge-rb-n32/raw/n32-p262144-r32-l16.jsonl.gz" \
  --auditor "$R31/r031-finite-compute-public/merge-rb-n32/source/auditor_frozen_r026.py" \
  --derived-source "$R31/r031-finite-compute-public/merge-rb-n32/source/production_r026.swift" \
  --recorded "$R31/r031-finite-compute-public/merge-rb-n32/reports/n32-p262144-r32-l16.independent-aggregate.json" \
  --receipt "$R31/r031-finite-compute-public/merge-rb-n32/receipt/n32-p262144-r32-l16.receipt.json" \
  --workdir "$OUT/merge-rb-work" --output "$OUT/merge-rb-public-audit.json"
```

This validates finite raw/input invariants only. The reported severe ESS/maxshare
imbalance and `Q_mRB` are not probability, mixing, cone, or accuracy results.

## C. n=7 orbit archival audit and parent full-count identity audit

```sh
python3 -B "$R31/r031-n7-orbit-public/audit_orbits.py" \
  --raw "$R31/r031-n7-orbit-public/raw" \
  --r29-verified "$R31/r031-n7-orbit-public/r29-verified" \
  --center "$R31/r031-n7-orbit-public/center/short-tail-n7-center-crt-audit-public.json" \
  --pratt "$R31/r031-n7-orbit-public/pratt/pratt-certificates.json" \
  --historical-audit "$R31/r031-n7-orbit-public/historical/n7-short-tail-orbits-crt-audit.json" \
  --historical-center "$R31/r031-n7-orbit-public/historical/short-tail-n7-center-crt-audit-original.json" \
  --output "$OUT/n7-orbit-audit.json"

python3 -B "$R31/r031-short-gap-parent/audit_n7_identities.py" \
  --short-gaps "$OUT/n7-orbit-audit.json" \
  --full-count "$R31/r031-short-gap-parent/n7-exact-r005.json" \
  --output "$OUT/n7-identities.json"
```

The checked-in parent result and its hash provenance are in
`r031-short-gap-parent/n7-identities.json` and
`r031-short-gap-parent/REPLAY_PROVENANCE.json`.

## D. Seven small finite proof/spatial checks

```sh
python3 -B "$R31/replay_proofs.py" \
  --certificates "$CERTIFICATES" --output "$OUT/proof-replay"
```

The included `r031-parent-proof-replay-receipt.json` records one already
completed replay: all seven scripts were run in isolated directories and their
complete JSON results matched. It is finite integer evidence, not an
asymptotic proof. The command above is the repeatable replay, not a reliance on
that receipt.

The computational replays above certify only their stated finite inputs and
outputs. Separately, the analytical documents in `r031-short-gap-identity/`
state and prove full-size local short-gap identities under their explicit
hypotheses; do not mislabel those local analytic statements as merely finite
experiments. Neither the finite computations nor those local identities settle
the original global concentration/asymptotic question.
