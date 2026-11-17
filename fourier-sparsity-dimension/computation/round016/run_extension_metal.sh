#!/bin/zsh
set -euo pipefail
root="$HOME/CodexRemoteMath/a28-selector/r016/extension-metal"
cache="$HOME/CodexRemoteMath/cache"
mkdir -p "$root" "$cache/swift-module-cache" "$cache/tmp"
export TMPDIR="$cache/tmp"
export PYTORCH_ENABLE_MPS_FALLBACK=0
cd "$root"
{ /usr/bin/swiftc --version; /usr/bin/sw_vers; echo "PYTORCH_ENABLE_MPS_FALLBACK=$PYTORCH_ENABLE_MPS_FALLBACK"; } > environment.txt
input_sha=$(shasum -a 256 extension-seed.json | awk '{print $1}')
source_sha=$(shasum -a 256 ExtensionMetal.swift | awk '{print $1}')
/usr/bin/swiftc -O -module-cache-path "$cache/swift-module-cache" -framework Foundation -framework Metal ExtensionMetal.swift -o extension_metal_verifier
./extension_metal_verifier extension-seed.json extension-metal-report.json "$input_sha" "$source_sha"
shasum -a 256 ExtensionMetal.swift extension-seed.json extension_metal_verifier extension-metal-report.json > manifest-sha256.txt
