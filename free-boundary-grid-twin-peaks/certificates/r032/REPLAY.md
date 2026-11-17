# Portable replay commands

Run from this directory. Set `OUT` to a writable directory outside this
package. Set `SYM_PYTHON` to an interpreter with SymPy installed; the successful
spatial replay used Python 3.13.18 and SymPy 1.14.0. Check the chosen runtime
explicitly before running the symbolic audit.

```sh
set -eu
R32="$PWD"
OUT="${OUT:?set an external writable output directory}"
SYM_PYTHON="${SYM_PYTHON:?set an interpreter that imports sympy}"
mkdir -p "$OUT"
```

For the analytic scripts, copy their directory before execution because their
JSON outputs are intentionally written beside the script:

```sh
cp -R "$R32/analytic-r45" "$OUT/analytic-r45"
for s in generic_fibre_nullspace.py exhaust_bipartite7_r45.py \
  search_bipartite_counterexamples.py check_bipartite9_r27.py; do
  python3 -B "$OUT/analytic-r45/$s"
done
cmp "$OUT/analytic-r45/exhaust-bipartite7-r45.json" \
  "$R32/analytic-r45/exhaust-bipartite7-r45.json"

"$SYM_PYTHON" -c 'import sys, sympy; print(sys.version); print(sympy.__version__)'
"$SYM_PYTHON" "$R32/spatial-r45/audit_interleaving_identity.py"
cp "$R32/independent-prefix-dp/verify_r45_prefix_dp.py" "$OUT/"
python3 -B "$OUT/verify_r45_prefix_dp.py"
```

The CRT auditors read saved raw checkpoints; they do not start n=6/n=7
production enumeration:

```sh
python3 -B "$R32/n6-h45/source/audit_n6_h45_portable.py" \
  --runs "$R32/n6-h45/raw/checkpoints" \
  --source "$R32/n6-h45/source/short_gap_h45_mod.cpp" \
  --pratt "$R32/common/pratt-certificates.json" \
  --punctured "$R32/n6-h45/references/punctured-table-n6-assert.json" \
  --top2 "$R32/n6-h45/references/top2-n6.json" --output "$OUT/n6-audit.json"

python3 -B "$R32/n7-center-h45/source/audit_n7_center_h45_portable.py" \
  --runs "$R32/n7-center-h45/raw/canonical" \
  --source "$R32/n7-center-h45/source/short_gap_h45_mod.cpp" \
  --pratt "$R32/common/pratt-certificates.json" \
  --old-center "$R32/n7-center-h45/references/old-center-h1-h3-public.json" \
  --normalizer "$R32/n7-center-h45/references/normalizer-v24.json" \
  --output "$OUT/n7-centre-audit.json"

python3 -B "$R32/short-gap-consequences/audit_short_gap_consequences.py" \
  --n6 "$OUT/n6-audit.json" --n7-center "$OUT/n7-centre-audit.json" \
  --output "$OUT/short-gap-consequences.json"
```

An optional source-build smoke, not an n=7 computation, is:

```sh
env CLANG_MODULE_CACHE_PATH="${CLANG_MODULE_CACHE_PATH:?set a shared writable compiler-module cache}" \
  clang++ -std=c++20 -O3 -Wall -Wextra -Werror -pthread \
  "$R32/n6-h45/source/short_gap_h45_mod.cpp" -o "$OUT/short_gap_h45_mod"
"$OUT/short_gap_h45_mod" 3 4 0 1 "$OUT/n3-v4-p0.json"
```
