#!/bin/zsh
set -euo pipefail
root="$HOME/CodexRemoteMath/a28-selector/r014-gpu"
cache="$HOME/CodexRemoteMath/cache"
mkdir -p "$root" "$cache/swift-module-cache" "$cache/tmp"
export TMPDIR="$cache/tmp"
export PYTORCH_ENABLE_MPS_FALLBACK=0
cd "$root"
{ /usr/bin/swiftc --version; /usr/bin/sw_vers; echo "PYTORCH_ENABLE_MPS_FALLBACK=$PYTORCH_ENABLE_MPS_FALLBACK"; } > environment.txt
input_sha=$(shasum -a 256 selector-input-canonical.json | awk '{print $1}')
source_sha=$(shasum -a 256 MetalWalshVerifier.swift | awk '{print $1}')
/usr/bin/swiftc -O -module-cache-path "$cache/swift-module-cache" -framework Foundation -framework Metal MetalWalshVerifier.swift -o metal_walsh_verifier
./metal_walsh_verifier selector-input-canonical.json metal-walsh-report.json "$input_sha" "$source_sha"
shasum -a 256 MetalWalshVerifier.swift selector-input-canonical.json metal_walsh_verifier metal-walsh-report.json > manifest-sha256.txt
