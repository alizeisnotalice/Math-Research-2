#!/bin/bash
set -euo pipefail
host=$1
first=$2
workers=$3
start=$4
output="extended_$host"
mkdir "$output"

run_cohorts() {
    planned=$start
    for input in extended_input/cohort*.txt; do
        cohort=$(/usr/bin/basename "$input" .txt)
        ./r047_weight_census "$input" 27 "$first" "$workers" "$workers" "$planned" \
            "$output/$cohort-weight.json" "$output/$cohort-weight-usage.csv"
        ./r047_flat_census "$input" 27 "$first" "$workers" "$workers" 0 \
            "$output/$cohort-flat.json" "$output/$cohort-flat-usage.csv"
        planned=0
        printf '%s completed %s\n' "$(/bin/date -u +'%H:%M:%S')" "$cohort"
    done
}

run_cohorts >"$output/completion.txt" 2>"$output/errors.txt" &
runner=$!
while kill -0 "$runner" 2>/dev/null; do
    {
        /usr/bin/perl -MTime::HiRes=time -e 'printf "%.6f\n",time'
        ps -Ao pid=,pcpu=,comm= -r | /usr/bin/awk 'NR<=12'
    } >>"$output/other_processes.txt"
    sleep .30
done
wait "$runner"
