#!/bin/bash
set -u
out=both
mkdir "$out" || exit 1
printf '%s\n' "$$" > "$out/supervisor.pid"
date -u +%FT%TZ > "$out/start_utc"
values=(20 28 36)
pids=()
for L in "${values[@]}"; do
    "$HOME/.venvs/math/bin/python" certify_both.py "$L" "$out/$L.json" > "$out/$L.log" 2>&1 &
    pids+=("$!")
done
printf '%s\n' "${pids[@]}" > "$out/pids.txt"
pid_csv=$(IFS=,; printf '%s' "${pids[*]}")
sleep .3
ps -p "$pid_csv" -o pid,ppid,%cpu,stat,time,command > "$out/ps.txt"
failed=0
for i in 0 1 2; do
    wait "${pids[i]}"
    rc=$?
    printf '%s %s %s\n' "${values[i]}" "${pids[i]}" "$rc" >> "$out/exits.tsv"
    if [ "$rc" -ne 0 ]; then failed=1; fi
done
ps -p "$pid_csv" -o pid,ppid,%cpu,stat,time,command > "$out/post-ps.txt"
date -u +%FT%TZ > "$out/end_utc"
printf '%s\n' "$failed" > "$out/status.txt"
exit "$failed"
