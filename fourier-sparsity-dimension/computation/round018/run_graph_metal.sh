#!/bin/zsh
set -euo pipefail
root="$HOME/CodexRemoteMath/a28-selector/r018/graph-metal"; cache="$HOME/CodexRemoteMath/cache";mkdir -p "$root" "$cache/swift-module-cache" "$cache/tmp";export TMPDIR="$cache/tmp";cd "$root"
{ /usr/bin/swiftc --version;/usr/bin/sw_vers; } > environment.txt
ih=$(shasum -a 256 graph-cases.json|awk '{print $1}');sh=$(shasum -a 256 GraphMetal.swift|awk '{print $1}')
/usr/bin/swiftc -O -module-cache-path "$cache/swift-module-cache" -framework Foundation -framework Metal GraphMetal.swift -o verifier
./verifier graph-cases.json report.json "$ih" "$sh"
shasum -a 256 GraphMetal.swift graph-cases.json verifier report.json > manifest-sha256.txt
