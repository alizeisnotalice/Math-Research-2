"""Row-transfer census for n=8 independent sets; it allocates no subset DP."""
import json

n = 8
valid = [s for s in range(1 << n) if not s & (s << 1)]
state = {(0, 0): 1}
for _ in range(n):
    nxt = {}
    for (previous, cardinality), count in state.items():
        for row in valid:
            if not row & previous:
                key = (row, cardinality + row.bit_count())
                nxt[key] = nxt.get(key, 0) + count
    state = nxt
histogram = {}
for (_, cardinality), count in state.items():
    histogram[cardinality] = histogram.get(cardinality, 0) + count
peak_k, peak_count = max(histogram.items(), key=lambda item: item[1])
print(json.dumps({
    'n': n,
    'valid_row_masks': len(valid),
    'independent_sets': sum(histogram.values()),
    'histogram': histogram,
    'largest_layer': {'k': peak_k, 'raw_count': peak_count,
                      'optimistic_D4_representatives_lower_bound': (peak_count + 7) // 8,
                      'U256_bytes_at_optimistic_D4_bound': 32*((peak_count + 7) // 8)},
}, indent=2))
