#!/bin/zsh
set -euo pipefail
root="$HOME/CodexRemoteMath/a28-selector/r015/q4-metal"
cache="$HOME/CodexRemoteMath/cache"
mkdir -p "$root" "$cache/swift-module-cache" "$cache/tmp"
export TMPDIR="$cache/tmp"
export PYTORCH_ENABLE_MPS_FALLBACK=0
cd "$root"
{ /usr/bin/swiftc --version; /usr/bin/sw_vers; echo "PYTORCH_ENABLE_MPS_FALLBACK=$PYTORCH_ENABLE_MPS_FALLBACK"; } > environment.txt
input_sha=$(shasum -a 256 q4-input.json | awk '{print $1}')
source_sha=$(shasum -a 256 MetalWalshVerifier_r015.swift | awk '{print $1}')
/usr/bin/swiftc -O -module-cache-path "$cache/swift-module-cache" -framework Foundation -framework Metal MetalWalshVerifier_r015.swift -o q4_metal_walsh_verifier
./q4_metal_walsh_verifier q4-input.json q4-metal-report.json "$input_sha" "$source_sha"
shasum -a 256 MetalWalshVerifier_r015.swift q4-input.json q4_metal_walsh_verifier q4-metal-report.json > manifest-sha256.txt
