#!/usr/bin/env python3
"""Exact carry-frontier ledger for H_d(x^a)H_d(x^b), b<=d.

For every injected target support S^2, enumerate every nonzero integer t with
(i+b*t, j-a*t) in (S+S)^2.  No floating arithmetic or algebraic root claim is
used.  The resulting flags only rank proof obstacles; they do not exclude a
candidate.
"""
from __future__ import annotations

import hashlib
import json
from math import gcd
from pathlib import Path


HERE = Path(__file__).resolve().parent


def support(d: int) -> tuple[int, ...]:
    return (0, 1, 2, d - 1, d)


def sumset(s: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted({i + j for i in s for j in s}))


def exact_events(d: int, a: int, b: int) -> list[dict[str, int]]:
    s = support(d)
    ss = set(sumset(s))
    # These bounds cover all t compatible with the finite r/s coordinate box.
    bound = 2 * d
    rows: list[dict[str, int]] = []
    for i in s:
        for j in s:
            for t in range(-bound, bound + 1):
                if not t:
                    continue
                r, v = i + b * t, j - a * t
                if r in ss and v in ss:
                    rows.append({"i": i, "j": j, "t": t, "u_degree": r, "v_degree": v})
    return rows


def target_is_injective(d: int, a: int, b: int) -> bool:
    s = support(d)
    values = [a * i + b * j for i in s for j in s]
    return len(set(values)) == len(values) == 25


def known_status(d: int, a: int, b: int) -> str:
    if (d, a, b) == (10, 3, 4):
        return "proved_by_parent_report_H10_a3_b4"
    return "not_covered_by_the_explicit_Hd_stretch_or_tensor_theorems"


def record(d: int, a: int, b: int) -> dict[str, object]:
    s = support(d)
    injected = target_is_injective(d, a, b)
    events = exact_events(d, a, b) if injected else []
    negative = [row for row in events if row["t"] < 0]
    over = [row for row in events if row["v_degree"] > d]
    low = [row for row in events if row["v_degree"] in (0, 1, 2)]
    return {
        "d": d, "a": a, "b": b, "total_degree": d * (a + b),
        "S": list(s), "S_plus_S": list(sumset(s)), "target_injective": injected,
        "status": known_status(d, a, b),
        "noncanonical_representations": events,
        "noncanonical_count": len(events),
        "has_negative_t": bool(negative), "has_vdegree_gt_d": bool(over),
        "has_low_v012": bool(low),
        "negative_t_count": len(negative), "vdegree_gt_d_count": len(over), "low_v012_count": len(low),
    }


def main() -> None:
    if not Path('/Volumes/data').is_mount():
        raise SystemExit('/Volumes/data is not mounted; refusing to write elsewhere')
    rows = []
    for d in range(6, 61):
        for a in range(2, d + 1):
            for b in range(a + 1, d + 1):
                if gcd(a, b) == 1:
                    rows.append(record(d, a, b))
    injected = [row for row in rows if row['target_injective']]
    difficult = [row for row in injected if row['status'].startswith('not_covered') and
                 (row['has_negative_t'] or row['has_vdegree_gt_d'] or row['has_low_v012'])]
    difficult.sort(key=lambda row: (row['total_degree'], row['d'], row['a'], row['b']))
    small_records = [[row['d'], row['a'], row['b'], row['total_degree'], row['target_injective'], row['status'],
                      row['noncanonical_count'], row['has_negative_t'], row['has_vdegree_gt_d'], row['has_low_v012'],
                      row['negative_t_count'], row['vdegree_gt_d_count'], row['low_v012_count']] for row in rows]
    all_carries = [[row['d'], row['a'], row['b'],
                    [[event['i'], event['j'], event['t'], event['u_degree'], event['v_degree']]
                     for event in row['noncanonical_representations']]] for row in injected]
    payload = {
        "scope": "d=6..60, 2<=a<b<=d, gcd(a,b)=1; exact integer support/carry audit",
        "coordinate_identity": "a*i+b*j = a*(i+b*t)+b*(j-a*t)",
        "target_support": "S^2 with S={0,1,2,d-1,d}",
        "source_coordinate_region": "(S+S)^2",
        "all_pairs_count": len(rows), "target_injective_count": len(injected),
        "target_noninjective_count": len(rows) - len(injected),
        "small_record_columns": ["d", "a", "b", "total_degree", "target_injective", "status", "noncanonical_count", "has_negative_t", "has_vdegree_gt_d", "has_low_v012", "negative_t_count", "vdegree_gt_d_count", "low_v012_count"],
        "records": small_records,
        "frontier_first_not_marked_closed": difficult[0] if difficult else None,
        "frontier_top_25": difficult[:25],
        "claim_boundary": "The flags rank carry/proof obstacles only.  They are not finite exclusions or proofs of nonuniform factorization.",
    }
    output = HERE / 'Hd-frontier-rev2-small-ledger.json'
    output.write_text(json.dumps(payload, separators=(',', ':')) + '\n')
    carry_output = HERE / 'Hd-frontier-rev2-all-carries.json'
    carry_output.write_text(json.dumps({"columns": ["d", "a", "b", "events[i,j,t,u_degree,v_degree]"], "records": all_carries}, separators=(',', ':')) + '\n')
    receipt = {
        "status": "completed_exact_integer_enumeration",
        "ledger": str(output), "ledger_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "all_carries": str(carry_output), "all_carries_sha256": hashlib.sha256(carry_output.read_bytes()).hexdigest(),
        "all_pairs_count": len(rows), "target_injective_count": len(injected),
        "first_not_marked_closed": difficult[0] if difficult else None,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "float_arithmetic_used": False,
    }
    (HERE / 'Hd-frontier-rev2-receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
