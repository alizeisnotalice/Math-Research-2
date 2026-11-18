#!/bin/sh
set -eu
prefix=$1
shift
export TMPDIR="$HOME/CodexRemoteMath/cache/cc192-temp"
printf 'start %s\n' "$(date -u +%FT%TZ)" > "$prefix.lifecycle.log"
./marchaud_scan "$@" "$prefix.json" > "$prefix.stdout.log" 2>&1 &
worker=$!
printf 'pid %s\n' "$worker" >> "$prefix.lifecycle.log"
while kill -0 "$worker" 2>/dev/null; do
    usage=$(ioreg -r -c IOAccelerator -d 1 | perl -nle 'print $1 if /"Device Utilization %"=([0-9]+)/')
    printf '%s %s\n' "$(date -u +%FT%TZ)" "$usage" >> "$prefix.gpu.log"
    {
        printf 'sample %s\n' "$(date -u +%FT%TZ)"
        ps -axo pid,pcpu,comm | sort -k2,2nr | head -18
        ps -p "$worker" -o pid,pcpu,comm || true
    } >> "$prefix.ps.log"
    sleep 0.5
done
if wait "$worker"; then code=0; else code=$?; fi
printf 'end %s exit %s\n' "$(date -u +%FT%TZ)" "$code" >> "$prefix.lifecycle.log"
exit "$code"
