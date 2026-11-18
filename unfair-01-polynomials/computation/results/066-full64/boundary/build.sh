#!/bin/sh
set -eu
cd "$(dirname "$0")"
/usr/bin/clang++ -std=c++20 -O3 -Wall -Wextra -Werror \
  -fmodules-cache-path=/Volumes/data/CodexTemp/cache/shared-clang-modules \
  ordinary_remainder.cpp -o ordinary_remainder
/usr/bin/clang++ -std=c++20 -O3 \
  -fmodules-cache-path=/Volumes/data/CodexTemp/cache/shared-clang-modules \
  reference_source.cpp -o reference_small
