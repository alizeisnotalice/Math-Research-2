# Replay

The integer audit never runs an n=7 production enumeration. It reads the 40
archived JSON/exit/resource triples, verifies the 49-node Pratt certificate,
performs direct cofactor and Garner CRT, and cross-checks h1--h3 against r31.

```sh
set -eu
PKG="$PWD"
OUT="${OUT:?set a writable directory outside the package}"
mkdir -p "$OUT"
python3 -B "$PKG/source/audit_n7_global_h45.py" \
  --raw "$PKG/raw" --source "$PKG/source/short_gap_h45_mod.cpp" \
  --pratt "$PKG/references/pratt-certificates.json" \
  --r31 "$PKG/references/r31-orbit-h1-h3-audit.json" \
  --output "$OUT/n7-global-h45-crt-audit.json"
```

For an optional small build smoke only (not n=7), use a shared compiler cache:

```sh
env CLANG_MODULE_CACHE_PATH="${CLANG_MODULE_CACHE_PATH:?set a shared writable cache}" \
  clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pthread \
  "$PKG/source/short_gap_h45_mod.cpp" -o "$OUT/short_gap_h45_mod"
"$OUT/short_gap_h45_mod" 3 4 0 1 "$OUT/n3-v4-p0.json"
```
