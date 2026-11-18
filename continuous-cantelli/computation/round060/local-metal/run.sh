#!/bin/zsh
set -eu
seconds=${1:-40}
order=${2:-96}
suffix=${3:-main}
mount | rg -q ' on /Volumes/data '
workdir=${0:A:h}
cd "$workdir"
export TMPDIR="$workdir" TMP="$workdir" TEMP="$workdir"
export XDG_CACHE_HOME=/Volumes/data/CodexTemp/cache
export CLANG_MODULE_CACHE_PATH=/Volumes/data/CodexTemp/cache/swift-clang
xcrun swiftc -O -module-cache-path "$CLANG_MODULE_CACHE_PATH" quartic.swift -o quartic
usage_line() {
    printf '%s %s ' "$1" "$(date -u +%FT%TZ)"
    ioreg -r -d 1 -c IOAccelerator |
        rg -o '"(Device Utilization %|Renderer Utilization %|Tiler Utilization %)"=[0-9]+' |
        paste -sd ' ' -
}
usage_line before > "ioreg-${suffix}.txt"
./quartic "$seconds" "$order" "$suffix" "${QUARTIC_START_EPOCH:-0}" "${QUARTIC_START_INDEX:-2147483648}" > "quartic-${suffix}.log" 2>&1 &
job=$!
while kill -0 "$job" 2>/dev/null
do
    phase=during
    if [[ "$(date +%s)" -lt "${QUARTIC_START_EPOCH:-0}" ]]
    then
        phase=waiting
    fi
    usage_line "$phase" >> "ioreg-${suffix}.txt"
    sleep 0.5
done
result_code=0
wait "$job" || result_code=$?
usage_line after >> "ioreg-${suffix}.txt"
printf 'exit=%s pid=%s\n' "$result_code" "$job" >> "quartic-${suffix}.log"
tail -5 "quartic-${suffix}.log"
exit "$result_code"
