#!/bin/sh
# Minimal source-only reproduction for commit 582cdcc4a17a84531a5bda6ee5ac8d3fc0ff7766.
# It intentionally does not attempt an unavailable data/proof-corpus replay.
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/ZeroOnePolynomials" && pwd)
OUT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)/repro_build
cd "$ROOT"
git rev-parse HEAD
shasum -a 256 LICENSE solver/ZeroOneSolver.cpp solver/ZeroOneSolver.hpp \
  solver/PositivstellensatzCertificates.jl solver/DistributedVerifier.jl \
  solver/data.sha512

mkdir -p "$OUT"
cd solver
c++ -std=c++20 -O2 -DZERO_ONE_SOLVER_M=5 -DZERO_ONE_SOLVER_N=7 \
  ZeroOneSolver.cpp -o "$OUT/tree_5_7"
"$OUT/tree_5_7" >"$OUT/tree_5_7.stdout" 2>"$OUT/tree_5_7.stderr"
c++ -std=c++20 -O2 -DZERO_ONE_SOLVER_M=7 -DZERO_ONE_SOLVER_N=13 \
  ZeroOneSolver.cpp -o "$OUT/tree_7_13"
"$OUT/tree_7_13" >"$OUT/tree_7_13.stdout" 2>"$OUT/tree_7_13.stderr"
wc -l -c "$OUT/tree_5_7.stdout" "$OUT/tree_7_13.stdout"
shasum -a 256 "$OUT/tree_5_7" "$OUT/tree_7_13"

if [ ! -d data ] || [ ! -d proofs ]; then
  echo "data/proofs absent: full certificate verification intentionally unavailable" >&2
fi
