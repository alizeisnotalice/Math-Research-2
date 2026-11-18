#!/bin/sh
set -eu
cd "$(dirname "$0")"
/usr/bin/clang++ -std=c++20 -O3 -Wall -Wextra -Werror \
  -fmodules-cache-path=/Volumes/data/CodexTemp/cache/shared-clang-modules sieve.cpp -o sieve
