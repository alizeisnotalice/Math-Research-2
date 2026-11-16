#!/usr/bin/env python3
"""Independent direct-coordinate audit of the compact H_d carry ledger."""
from __future__ import annotations
import hashlib
import json
from math import gcd
from pathlib import Path


HERE = Path(__file__).resolve().parent


def direct_events(d: int, a: int, b: int) -> tuple[bool, list[list[int]]]:
    s = (0, 1, 2, d-1, d)
    ss = sorted({i+j for i in s for j in s})
    target: dict[int, tuple[int, int]] = {}
    injected = True
    for i in s:
        for j in s:
            exponent = a*i + b*j
            if exponent in target:
                injected = False
            target[exponent] = (i, j)
    if not injected:
        return False, []
    events = set()
    for r in ss:
        for v in ss:
            exponent = a*r + b*v
            if exponent not in target:
                continue
            i, j = target[exponent]
            delta = r - i
            assert delta % b == 0
            t = delta // b
            assert v == j - a*t
            if t:
                events.add((i, j, t, r, v))
    return True, [list(row) for row in sorted(events)]


def main() -> None:
    ledger_path = HERE/'Hd-frontier-rev2-small-ledger.json'
    carries_path = HERE/'Hd-frontier-rev2-all-carries.json'
    ledger = json.loads(ledger_path.read_text())
    carries = json.loads(carries_path.read_text())
    carry_map = {(d, a, b): events for d, a, b, events in carries['records']}
    assert len(carry_map) == len(carries['records']), 'duplicate carry-map key'
    columns = ledger['small_record_columns']
    assert columns[:7] == ['d', 'a', 'b', 'total_degree', 'target_injective', 'status', 'noncanonical_count']
    expected_keys = {(d, a, b) for d in range(6, 61) for a in range(2, d + 1)
                     for b in range(a + 1, d + 1) if gcd(a, b) == 1}
    seen_keys = set()
    injected_keys = set()
    checked = 0
    for row in ledger['records']:
        d, a, b, _, injected, _, count, neg, over, low, neg_count, over_count, low_count = row
        assert 6 <= d <= 60 and 2 <= a < b <= d and gcd(a, b) == 1
        key = (d, a, b)
        assert key not in seen_keys, f'duplicate ledger key: {key}'
        seen_keys.add(key)
        actual_injected, events = direct_events(d, a, b)
        assert actual_injected == injected
        if injected:
            injected_keys.add(key)
            assert events == carry_map[key]
            assert count == len(events)
            assert neg == any(event[2] < 0 for event in events)
            assert over == any(event[4] > d for event in events)
            assert low == any(event[4] in (0, 1, 2) for event in events)
            assert neg_count == sum(event[2] < 0 for event in events)
            assert over_count == sum(event[4] > d for event in events)
            assert low_count == sum(event[4] in (0, 1, 2) for event in events)
        else:
            assert key not in carry_map
        checked += 1
    assert seen_keys == expected_keys, (len(seen_keys), len(expected_keys), expected_keys - seen_keys, seen_keys - expected_keys)
    assert set(carry_map) == injected_keys, (len(carry_map), len(injected_keys), injected_keys - set(carry_map), set(carry_map) - injected_keys)
    receipt = {
        'status': 'all_pairs_direct_coordinate_audit_passed', 'pairs_checked': checked,
        'small_ledger_sha256': hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
        'all_carries_sha256': hashlib.sha256(carries_path.read_bytes()).hexdigest(),
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'float_arithmetic_used': False,
    }
    (HERE/'Hd-frontier-rev2-independent-verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__': main()
