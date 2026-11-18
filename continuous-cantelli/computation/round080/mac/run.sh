#!/bin/bash
set -u
cd "$HOME/CodexRemoteMath/cc79-cpu-nov18"
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
run="hermite"
mkdir "$run" || exit 1
printf '%s\n' "$$" > "$run/supervisor.pid"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$run/start_utc"
beg=$(date +%s)
pids=()
for ((i=0; i<15; i++)); do
    "$HOME/.venvs/math/bin/python" hermite_tail.py "$i" --exponent 4 --seconds 100 --out "$run/shard-$i.json" > "$run/shard-$i.log" 2>&1 &
    pids+=("$!")
done
printf '%s\n' "${pids[@]}" > "$run/worker-pids.txt"
pid_csv=$(IFS=,; printf '%s' "${pids[*]}")
for ((snapshot=1; snapshot<=5; snapshot++)); do
    sleep 3
    ps -p "$pid_csv" -o pid,ppid,%cpu,stat,time,rss,command > "$run/ps-$snapshot.txt"
done
top -l 3 -n 0 -s 2 > "$run/top.txt"
failed=0
for ((i=0; i<15; i++)); do
    wait "${pids[i]}"
    code=$?
    printf '%s %s %s\n' "$i" "${pids[i]}" "$code" >> "$run/exits.tsv"
    if (( code != 0 )); then failed=1; fi
done
ps -p "$pid_csv" -o pid,ppid,%cpu,stat,time,rss,command > "$run/post-exit-ps.txt"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$run/end_utc"
printf 'wall_seconds %s\nexit_code %s\n' "$(( $(date +%s) - beg ))" "$failed" > "$run/status.txt"
exit "$failed"
