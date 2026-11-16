#!/bin/zsh
# Usage: monitor_agx_rss.sh PID CSV_PATH.  AGX values are whole-system,
# non-exclusive IOKit counters; powermetrics was unavailable without root.
set -eu
pid="$1"
out="$2"
printf 'epoch_s,pid,rss_kb,agx_device_util_pct,agx_renderer_util_pct,agx_tiler_util_pct\n' > "$out"
while kill -0 "$pid" 2>/dev/null; do
  now=$(date +%s)
  rss=$(ps -o rss= -p "$pid" 2>/dev/null | tr -d ' ')
  stats=$(ioreg -l -w0 -r -c AGXAcceleratorG18G 2>/dev/null | /usr/bin/sed -n '/"PerformanceStatistics"/p' | /usr/bin/head -n 1)
  dev=$(print -r -- "$stats" | /usr/bin/sed -n 's/.*"Device Utilization %"=\([0-9][0-9]*\).*/\1/p')
  renderer=$(print -r -- "$stats" | /usr/bin/sed -n 's/.*"Renderer Utilization %"=\([0-9][0-9]*\).*/\1/p')
  tiler=$(print -r -- "$stats" | /usr/bin/sed -n 's/.*"Tiler Utilization %"=\([0-9][0-9]*\).*/\1/p')
  printf '%s,%s,%s,%s,%s,%s\n' "$now" "$pid" "${rss:-}" "${dev:-}" "${renderer:-}" "${tiler:-}" >> "$out"
  sleep 1
done
