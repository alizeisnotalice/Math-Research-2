#!/bin/bash
set -u
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
out=broad
mkdir "$out" || exit 1
printf '%s\n' "$$" > "$out/supervisor.pid"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$out/queued_utc"
pids=()
for ((i=0; i<15; i++)); do
    "$HOME/.venvs/math/bin/python" broad.py "$i" "$out/shard-$i.json" > "$out/shard-$i.log" 2>&1 &
    pids+=("$!")
done
printf '%s\n' "${pids[@]}" > "$out/worker-pids.txt"
pid_csv=$(IFS=,; printf '%s' "${pids[*]}")
for ((sample=1; sample<=9; sample++)); do
    date -u '+%Y-%m-%dT%H:%M:%SZ' >> "$out/timeline.txt"
    ps -p "$pid_csv" -o pid,ppid,%cpu,stat,time,rss,command >> "$out/timeline.txt"
    sleep 7
    if (( sample == 5 )); then top -l 3 -n 0 -s 1 > "$out/top.txt"; fi
done
failed=0
for ((i=0; i<15; i++)); do
    wait "${pids[i]}"
    code=$?
    printf '%s %s %s\n' "$i" "${pids[i]}" "$code" >> "$out/exits.tsv"
    if (( code != 0 )); then failed=1; fi
done
ps -p "$pid_csv" -o pid,ppid,%cpu,stat,time,rss,command > "$out/post-exit-ps.txt"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$out/end_utc"
printf '%s\n' "$failed" > "$out/status.txt"
exit "$failed"
