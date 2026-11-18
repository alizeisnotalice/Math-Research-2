#!/bin/bash
set -u
out=run
mkdir "$out" || exit 1
printf '%s\n' "$$" > "$out/supervisor.pid"
date -u +%FT%TZ > "$out/start_utc"
lefts=(-1 0 .5)
rights=(0 1 1)
pids=()
tags=()
for i in 0 1 2; do
    for w in 1 -1; do
        tag="$i-$w"
        "$HOME/.venvs/math/bin/python" tail_cost.py "${lefts[i]}" "${rights[i]}" "$w" "$out/$tag.json" > "$out/$tag.log" 2>&1 &
        pids+=("$!")
        tags+=("$tag")
    done
done
printf '%s\n' "${pids[@]}" > "$out/pids.txt"
pid_csv=$(IFS=,; printf '%s' "${pids[*]}")
sleep .3
ps -p "$pid_csv" -o pid,ppid,%cpu,stat,time,command > "$out/ps.txt"
failed=0
for ((i=0;i<6;i++)); do
    wait "${pids[i]}"
    rc=$?
    printf '%s %s %s\n' "${tags[i]}" "${pids[i]}" "$rc" >> "$out/exits.tsv"
    if [ "$rc" -ne 0 ]; then failed=1; fi
done
ps -p "$pid_csv" -o pid,ppid,stat,command > "$out/post-ps.txt"
"$HOME/.venvs/math/bin/python" symbolic.py > "$out/symbolic.log" 2>&1
code=$?
printf '%s\n' "$code" > "$out/symbolic-exit.txt"
if [ "$code" -ne 0 ]; then failed=1; fi
date -u +%FT%TZ > "$out/end_utc"
printf '%s\n' "$failed" > "$out/status.txt"
exit "$failed"
