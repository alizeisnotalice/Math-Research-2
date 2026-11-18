#!/bin/bash
set -eu
mkdir -p temp
export TMPDIR="$PWD/temp" TMP="$PWD/temp" TEMP="$PWD/temp" XDG_CACHE_HOME="$PWD/temp"
export CLANG_MODULE_CACHE_PATH="$HOME/CodexRemoteMath/a28-selector/cache/clang"
printf 'Started UTC: '
/bin/date -u +'%Y-%m-%dT%H:%M:%SZ'
/usr/bin/swiftc --version
/usr/bin/swiftc -O -module-cache-path "$HOME/CodexRemoteMath/cache/swift-module-cache" R46MixedParityBlockerMetal.swift -o r46_mixed_parity_blocker_metal
/usr/bin/time -p ./r46_mixed_parity_blocker_metal macbook_r46_mixed_parity_result.json
/usr/bin/shasum -a 256 R46MixedParityBlockerMetal.swift run_r46.sh macbook_r46_mixed_parity_result.json
if /bin/ps -Ao comm= | /usr/bin/awk '$1 ~ /r46_mixed_parity_blocker_metal$/ {found=1} END{exit !found}'; then exit 2; fi
printf 'Verifier process absent after exit\n'
