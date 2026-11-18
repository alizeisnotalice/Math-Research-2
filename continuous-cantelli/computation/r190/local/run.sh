#!/bin/zsh
set -eu
mount | rg -q ' on /Volumes/data '
cd /Volumes/data/CodexTemp/T2/cc184-cantor-metal-local
export TMPDIR="$PWD" TMP="$PWD" TEMP="$PWD"
export XDG_CACHE_HOME=/Volumes/data/CodexTemp/cache
export CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-clang
mode=$1
xcrun swiftc -O -module-cache-path "$CLANG_MODULE_CACHE_PATH" cantor.swift -o cantor
./cantor "$mode" > "metal-$mode.log" 2>&1 &
job=$!
code=0
wait "$job" || code=$?
printf 'exit=%s pid=%s\n' "$code" "$job" >> "metal-$mode.log"
tail -3 "metal-$mode.log"
exit "$code"
