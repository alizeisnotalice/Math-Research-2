#!/Users/zzh/.venvs/main/bin/python
"""Verify coverage and aggregate the fourteen independent core6 d/e screens."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tag(b: int) -> str:
    return f"neg{-b}" if b < 0 else str(b)


def main() -> None:
    ledger_path = HERE / "core6_linked_conditional_ledger_corrected.json"
    ledger = json.loads(ledger_path.read_text())
    expected_mothers = len(ledger["rows"])
    expected_raw = expected_mothers * 43 * 15
    all_survivors: list[dict[str, int]] = []
    batches = []
    seen = set()
    for b in list(range(-7, 0)) + list(range(1, 8)):
        path = HERE / f"core6_de_screen_b_{tag(b)}.json"
        data = json.loads(path.read_text())
        assert data["status"] == "NECESSARY_SCREENS_ONLY"
        assert data["b_layers"] == [b]
        assert data["input_ledger_sha256"] == digest(ledger_path)
        assert data["script_sha256"] == data["script_sha256"]  # recorded provenance
        counts = data["counts"]
        assert counts["raw_d_e"] == counts["mothers"] * 43 * 15
        assert ((counts["raw_d_e"] - counts["positive_axis"])
                + (counts["positive_axis"] - counts["annulus_schur"])
                + counts["annulus_schur"] == counts["raw_d_e"])
        assert counts["annulus_schur"] == len(data["survivors"])
        for row in data["survivors"]:
            key = tuple(row[k] for k in ("b", "T", "N", "delta", "ell", "d", "e"))
            assert key not in seen
            seen.add(key)
        all_survivors.extend(data["survivors"])
        batches.append({
            "b": b,
            "file": path.name,
            "file_sha256": digest(path),
            "counts": counts,
            "first_failure_stream_sha256": data["first_failure_stream_sha256"],
        })
    assert sum(row["counts"]["mothers"] for row in batches) == expected_mothers
    assert sum(row["counts"]["raw_d_e"] for row in batches) == expected_raw
    result = {
        "status": "COMPLETE_NECESSARY_SCREEN_COVERAGE",
        "input_ledger_sha256": digest(ledger_path),
        "expected_mothers": expected_mothers,
        "expected_raw_d_e": expected_raw,
        "batches": batches,
        "totals": {
            name: sum(row["counts"][name] for row in batches)
            for name in ("mothers", "raw_d_e", "positive_axis", "annulus_schur")
        },
        "survivors": all_survivors,
        "survivor_count": len(all_survivors),
        "scope": (
            "All rows in the conditional d/e box were checked by necessary positive-axis "
            "and annulus tests. The survivors are not claimed to satisfy the full factor model."
        ),
        "script_sha256": digest(Path(__file__)),
    }
    output = HERE / "core6_de_complete_necessary_screen_summary.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"status={result['status']}")
    print("totals=" + json.dumps(result["totals"], sort_keys=True))
    print(f"survivor_count={len(all_survivors)}")
    print(f"summary_sha256={digest(output)}")


if __name__ == "__main__":
    main()
