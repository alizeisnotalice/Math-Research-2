#!/bin/zsh
# Reproduction of the round005 MacBook search, with corrected --out argument.
# Preconditions: the remote fast verifier has already been installed by the
# coordinating agent at $ROOT/verify_matching_fast.py.
set -euo pipefail
ROOT=/Users/zhengzhihao/CodexRemoteMath/a28-selector/r005
DEST="$ROOT/remote-search"
HERE=${0:A:h}

base64 -i "$HERE/matching_search.cpp" | \
  ssh -o BatchMode=yes -o StrictHostKeyChecking=yes -o ConnectTimeout=30 macbook-math '
    set -e
    ROOT=/Users/zhengzhihao/CodexRemoteMath/a28-selector/r005
    DEST=$ROOT/remote-search
    mkdir -p "$DEST"
    base64 -D > "$DEST/matching_search.cpp"
    cd "$DEST"
    nice -n 10 clang++ -O3 -std=c++20 matching_search.cpp -o matching_search
    shasum -a 256 matching_search.cpp matching_search
    nice -n 10 ./matching_search 2026111713 500000 1.0 result.json base
    nice -n 10 python3 "$ROOT/verify_matching_fast.py" result.json --out verify.json
    cat verify.json
    shasum -a 256 result.json verify.json
    gzip -c result.json | base64
  '
