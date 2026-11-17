#!/bin/zsh
set -euo pipefail
root="$HOME/CodexRemoteMath/a28-selector/r017/nonparallel-metal"
cache="$HOME/CodexRemoteMath/cache"
mkdir -p "$root" "$cache/swift-module-cache" "$cache/tmp"
export TMPDIR="$cache/tmp"
export PYTORCH_ENABLE_MPS_FALLBACK=0
cd "$root"
{ /usr/bin/swiftc --version; /usr/bin/sw_vers; echo "PYTORCH_ENABLE_MPS_FALLBACK=$PYTORCH_ENABLE_MPS_FALLBACK"; } > environment.txt
input_sha=$(shasum -a 256 nonparallel-cases.json | awk '{print $1}')
source_sha=$(shasum -a 256 NonparallelMetal.swift | awk '{print $1}')
/usr/bin/swiftc -O -module-cache-path "$cache/swift-module-cache" -framework Foundation -framework Metal NonparallelMetal.swift -o nonparallel_metal_verifier
./nonparallel_metal_verifier nonparallel-cases.json nonparallel-metal-report.json "$input_sha" "$source_sha"
shasum -a 256 NonparallelMetal.swift nonparallel-cases.json nonparallel_metal_verifier nonparallel-metal-report.json > manifest-sha256.txt
