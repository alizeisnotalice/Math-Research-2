#!/bin/zsh
set -eu
seconds=${1:-105}
order=${2:-36}
suffix=${3:-main}
mount | rg -q ' on /Volumes/data '
workdir=${0:A:h}
cd "$workdir"
export TMPDIR="$workdir" TMP="$workdir" TEMP="$workdir"
export XDG_CACHE_HOME=/Volumes/data/CodexTemp/cache
export CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-clang
xcrun swiftc -O -module-cache-path "$CLANG_MODULE_CACHE_PATH" scan.swift -o scan
usage_line() {
    printf '%s %s ' "$1" "$(date -u +%FT%TZ)"
    ioreg -r -d 1 -c IOAccelerator |
        rg -o '"(Device Utilization %|Renderer Utilization %|Tiler Utilization %)"=[0-9]+' |
        paste -sd ' ' -
}
usage_line before > "ioreg-${suffix}.txt"
./scan "$seconds" "$order" "$suffix" "${SCAN_START_EPOCH:-0}" "${SCAN_START_INDEX:-0}" > "scan-${suffix}.log" 2>&1 &
job=$!
while kill -0 "$job" 2>/dev/null
do
    phase=during
    if [[ "$(date +%s)" -lt "${SCAN_START_EPOCH:-0}" ]]
    then
        phase=waiting
    fi
    usage_line "$phase" >> "ioreg-${suffix}.txt"
    sleep 0.5
done
result_code=0
wait "$job" || result_code=$?
usage_line after >> "ioreg-${suffix}.txt"
printf 'exit=%s pid=%s\n' "$result_code" "$job" >> "scan-${suffix}.log"
tail -5 "scan-${suffix}.log"
exit "$result_code"
