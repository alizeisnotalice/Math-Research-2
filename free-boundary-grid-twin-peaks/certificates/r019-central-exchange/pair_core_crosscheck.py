#!/usr/bin/env python3
"""Cross-check n=4 same/split totals using the project's prescribed-pair DP."""
import json
import subprocess
import sys

exe = sys.argv[1]
n = 4
m = n // 2
same = split = total = 0
for p in range(n * n):
    for q in range(p + 1, n * n):
        if abs(p % n - q % n) + abs(p // n - q // n) == 1:
            continue
        k = int(subprocess.check_output([exe, str(n), str(p), str(q)], text=True))
        total += k
        if ((p % n) < m) == ((q % n) < m):
            same += k
        else:
            split += k
assert (total, same, split) == (779299174080, 191585447680, 587713726400)
print(json.dumps({'status': 'passed', 'A2': str(total), 'same': str(same),
                  'split': str(split)}, indent=2))
