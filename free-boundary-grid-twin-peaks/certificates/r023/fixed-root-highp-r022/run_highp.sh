#!/bin/zsh
set -euo pipefail
n=11; P=262144; R=256; L=16; T=12
# Vertex IDs are generated from these requested coordinates, never hard-coded.
coords=('0,0' '1,0' '1,1')
prevx=-1; prevy=-1
for xy in $coords; do
  x=${xy%,*}; y=${xy#*,}; root=$((y*n+x))
  [[ "$root" -eq $((y*n+x)) ]]
  if (( prevx >= 0 )); then
    dx=$((x-prevx)); (( dx < 0 )) && dx=$((-dx))
    dy=$((y-prevy)); (( dy < 0 )) && dy=$((-dy))
    (( dx + dy == 1 ))
  fi
  prevx=$x; prevy=$y
  out="n11-highp-root${root}-p${P}-r${R}-l${L}.csv"
  [[ ! -e "$out" && ! -e "${out}.metadata.json" ]]
  # Blocks are separated by 2^36, vastly exceeding all within-run offsets.
  seed=$((0xA0228000 + root*0x1000000000))
  ./fixed_root_smc "$n" "$root" "$P" "$R" "$L" "$T" "$out" "$seed"
  printf 'complete root=%s xy=(%s,%s) seed=%s\n' "$root" "$x" "$y" "$seed"
done
