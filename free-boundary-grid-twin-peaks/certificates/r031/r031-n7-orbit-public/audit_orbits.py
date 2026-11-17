#!/usr/bin/env python3
"""Portable finite CRT audit for the nine non-center n=7 D4 representatives.

The program verifies archived outputs; it does not run the historical C++
enumerator. Its command line has no machine-specific default paths.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from math import factorial, gcd, prod
from pathlib import Path
import json
import re

SOURCE_SHA = "9d757fcfd65729242d2742cc94b9c34e3e88ba23d1052b985c3b70f0b8a02d13"
HISTORICAL_BINARY_SHA = "e05c577bae9fb899a9c5789a6b8e2e6f4d0008ae0ae79fdf02cd64e17917ecc2"
HISTORICAL_CENTER_SHA = "02b46f419cbe8497855b336e959828e27908cfc54b0d4f9e1b941656af4b5293"
REPS = (0, 7, 8, 14, 15, 16, 21, 22, 23)
MULT = {0: 4, 7: 8, 8: 4, 14: 8, 15: 8, 16: 4, 21: 4, 22: 4, 23: 4, 24: 1}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_digest(value: object) -> str:
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def fraction(x: int, y: int) -> dict:
    z = Fraction(x, y)
    return {"numerator": str(z.numerator), "denominator": str(z.denominator)}


def pratt_roots(pratt_path: Path) -> tuple[list[int], int]:
    data = json.loads(pratt_path.read_text())
    nodes, seen = data["nodes"], set()

    def visit(p: int) -> None:
        if p in seen:
            return
        node = nodes.get(str(p))
        assert node is not None and int(node["p"]) == p
        if p == 2:
            assert node == {"p": "2", "base": "two"}
            seen.add(2)
            return
        fac = {int(q): int(e) for q, e in node["p_minus_1_factorization"].items()}
        value = 1
        for q, e in fac.items():
            assert q >= 2 and e >= 1
            visit(q)
            value *= q**e
        assert value == p - 1
        a = int(node["witness"])
        assert 1 < a < p and pow(a, p - 1, p) == 1
        for q in fac:
            assert gcd(pow(a, (p - 1) // q, p) - 1, p) == 1
        seen.add(p)

    roots = [int(x) for x in data["roots"]]
    for root in roots:
        visit(root)
    return roots, len(seen)


def direct_crt(res: list[int], primes: list[int]) -> int:
    q = prod(primes)
    return sum(a * (q // p) * pow(q // p, -1, p) for a, p in zip(res, primes)) % q


def garner_crt(res: list[int], primes: list[int]) -> int:
    x, m = 0, 1
    for a, p in zip(res, primes):
        x += m * (((a - x) % p) * pow(m, -1, p) % p)
        m *= p
    assert 0 <= x < m
    return x


def qe8(B: int, T: int, t: int) -> dict:
    lhs, rhs = t * (t + 1) * T, 8 * B
    return {"t": t, "T": str(T), "T_over_B": fraction(T, B),
            "target_8_over_t_tplus1": fraction(8, t * (t + 1)),
            "lhs_t_tplus1_T": str(lhs), "rhs_8B": str(rhs),
            "slack_rhs_minus_lhs": str(rhs - lhs), "passes": lhs <= rhs}


def centre_semantics(center: dict) -> dict:
    """Remove only documented I/O-location fields before comparing r030 forms."""
    result = json.loads(json.dumps(center))
    result["pratt"]["certificate_path"] = Path(result["pratt"]["certificate_path"]).name
    result["reference_selected_normalizer"]["path"] = Path(result["reference_selected_normalizer"]["path"]).name
    for row in result["raw_files"]:
        for key in ("json", "exit", "resource"):
            row[key] = Path(row[key]).name
    return result


def check_historical(current: dict, historical_path: Path) -> dict:
    """Compare mathematical fields while intentionally normalizing path strings."""
    old = json.loads(historical_path.read_text())
    assert old["status"] == "passed"
    assert current["crt"] == old["crt"]
    assert current["global"] == old["global"]
    assert current["local_qe8_t2_to_4"] == old["local_qe8_t2_to_4"]
    assert current["local_failures"] == old["local_failures"]
    assert current["pratt"]["roots"] == old["pratt"]["roots"]
    assert current["pratt"]["verified_nodes"] == old["pratt"]["verified_nodes"]
    assert current["pratt"]["certificate_sha256"] == old["pratt"]["certificate_sha256"]
    def orbit_payload(orbits: dict) -> dict:
        result = json.loads(json.dumps(orbits))
        for row in result.values():
            # The verified-file content hash is mathematical provenance; only
            # its historical absolute location is intentionally discarded.
            row.pop("reference_path")
        return result
    assert orbit_payload(current["orbits"]) == orbit_payload(old["orbits"])
    old_rows = {(x["v"], x["prime_index"]): {k: v for k, v in x.items() if k not in {"json", "exit", "resource"}}
                for x in old["raw_checkpoints"]}
    new_rows = {(x["v"], x["prime_index"]): {k: v for k, v in x.items() if k not in {"json", "exit", "resource"}}
                for x in current["raw_checkpoints"]}
    assert new_rows == old_rows
    return {"historical_audit_path": "historical/n7-short-tail-orbits-crt-audit.json",
            "historical_audit_sha256": digest(historical_path),
            "path_metadata_normalized": ["raw checkpoint json/exit/resource paths", "orbit verified-file paths"],
            "verified": ["CRT modulus, bounds, and both reconstruction methods",
                         "all ten orbit A2v/Bv/h/T/multiplicity values and QE8 rows",
                         "all 36 raw receipt SHA-256 values",
                         "Pratt roots, node count, and certificate SHA-256"]}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw", type=Path, required=True, help="directory containing exit.json, started.json, source/, and results/")
    p.add_argument("--r29-verified", type=Path, required=True, help="directory containing v*-verified.json and n7-orbits-verified.json")
    p.add_argument("--center", type=Path, required=True, help="r030 public center CRT JSON")
    p.add_argument("--pratt", type=Path, required=True, help="Pratt certificate JSON")
    p.add_argument("--output", type=Path, required=True, help="write the audited public JSON here")
    p.add_argument("--historical-audit", type=Path, help="historical r031 audit used for replay comparison")
    p.add_argument("--historical-center", type=Path, help="historical r030 center JSON for I/O-only publicization comparison")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    raw, r29, center_path, pratt_path, output = (args.raw.resolve(), args.r29_verified.resolve(), args.center.resolve(), args.pratt.resolve(), args.output.resolve())
    primes, pratt_nodes = pratt_roots(pratt_path)
    Q, bound = prod(primes), factorial(48)
    assert Q > bound and json.loads((raw / "exit.json").read_text())["rc"] == 0
    source_path = raw / "source" / "short_tail_mod.cpp"
    assert digest(source_path) == SOURCE_SHA
    raw_started = json.loads((raw / "started.json").read_text())
    assert raw_started["source"] == SOURCE_SHA and tuple(raw_started["v"]) == REPS

    outputs, local_qe8, provenance, layers, r29_sources = {}, {}, [], {}, {}
    for v in REPS:
        ref_path = r29 / f"v{v}-verified.json"
        ref = json.loads(ref_path.read_text())
        assert ref["verified"] and ref["v"] == v
        hist = list(map(int, ref["independent_set_layer_histogram"]))
        assert sum(hist) > 0
        expected_sets, expected_max = sum(hist), max(hist)
        assert ref["source_sha256"]
        layers[str(v)], r29_sources[str(v)] = hist, str(ref["source_sha256"])
        rows = []
        for pi, p in enumerate(primes):
            base = raw / "results" / f"v{v}" / f"p{pi}"
            jp, ep, rp = base.with_suffix(".json"), base.with_suffix(".exit.json"), base.with_suffix(".resource.txt")
            row, exit_row, resource = json.loads(jp.read_text()), json.loads(ep.read_text()), rp.read_text()
            assert exit_row["rc"] == 0 and exit_row["out"] == digest(jp) and exit_row["resource"] == digest(rp)
            assert re.search(r"\n\s*0\s+swaps\n", resource)
            assert row["n"] == 7 and row["v"] == v and row["prime_index"] == pi and int(row["prime"]) == p
            assert row["algorithm"] == "proposed-short-tail-modular-layer-rolling"
            assert row["formula_status"] == "pending independent proof"
            assert row["independent_sets"] == expected_sets and row["max_layer_count"] == expected_max
            rows.append(row)
            provenance.append({"v": v, "prime_index": pi, "json": f"raw/results/v{v}/p{pi}.json", "json_sha256": digest(jp),
                               "exit": f"raw/results/v{v}/p{pi}.exit.json", "exit_sha256": digest(ep),
                               "resource": f"raw/results/v{v}/p{pi}.resource.txt", "resource_sha256": digest(rp)})
        h = {}
        for key in ("h1", "h2", "h3"):
            residues = [int(row[f"{key}_residue"]) for row in rows]
            a, b = direct_crt(residues, primes), garner_crt(residues, primes)
            assert a == b and 0 <= a <= bound and all(a % p == residue for p, residue in zip(primes, residues))
            h[key] = a
        B, A = int(ref["Bv"]), int(ref["A2v"])
        T = {2: B - h["h1"], 3: B - h["h1"] - h["h2"], 4: B - h["h1"] - h["h2"] - h["h3"]}
        assert all(x >= 0 for x in T.values())
        outputs[v] = {"multiplicity": MULT[v], "A2v": A, "Bv": B, "h": h, "T": T,
                      "reference_path": f"r29-verified/v{v}-verified.json", "reference_sha256": digest(ref_path),
                      "independent_sets": expected_sets, "max_layer_count": expected_max}
        local_qe8[str(v)] = [qe8(B, T[t], t) for t in (2, 3, 4)]

    center = json.loads(center_path.read_text())
    assert center["status"] == "passed" and center["n"] == 7 and center["v"] == 24
    assert center["source_sha256_verified"] == SOURCE_SHA and center["pratt"]["certificate_sha256"] == digest(pratt_path)
    historical_center_path = (args.historical_center or output.parent / "historical" / "short-tail-n7-center-crt-audit-original.json").resolve()
    assert digest(historical_center_path) == HISTORICAL_CENTER_SHA
    historical_center = json.loads(historical_center_path.read_text())
    # This equality is deliberately stronger than checking selected numerical
    # fields: all non-path JSON material, including the CRT receipts, agrees.
    assert centre_semantics(center) == centre_semantics(historical_center)
    ref24 = json.loads((r29 / "v24-verified.json").read_text())
    h24 = {key: int(center["h"][key]) for key in ("h1", "h2", "h3")}
    B24, A24 = int(center["reference_selected_normalizer"]["Bv"]), int(center["reference_selected_normalizer"]["A2v"])
    assert B24 == int(ref24["Bv"]) and A24 == int(ref24["A2v"])
    T24 = {2: B24 - h24["h1"], 3: B24 - h24["h1"] - h24["h2"], 4: B24 - h24["h1"] - h24["h2"] - h24["h3"]}
    outputs[24] = {"multiplicity": 1, "A2v": A24, "Bv": B24, "h": h24, "T": T24,
                   "reference_path": "center/short-tail-n7-center-crt-audit-public.json", "reference_sha256": HISTORICAL_CENTER_SHA,
                   "center_reused_from_r30": True}
    local_qe8["24"] = [qe8(B24, T24[t], t) for t in (2, 3, 4)]
    assert sum(x["multiplicity"] for x in outputs.values()) == 49
    global_h = {key: sum(x["multiplicity"] * x["h"][key] for x in outputs.values()) for key in ("h1", "h2", "h3")}
    global_A = sum(x["multiplicity"] * x["A2v"] for x in outputs.values())
    global_B = sum(x["multiplicity"] * x["Bv"] for x in outputs.values())
    table = json.loads((r29 / "n7-orbits-verified.json").read_text())
    assert global_A == int(table["sum_A2v"]) and global_B == int(table["sum_Bv"])
    global_T = {2: global_B - global_h["h1"], 3: global_B - global_h["h1"] - global_h["h2"], 4: global_B - global_h["h1"] - global_h["h2"] - global_h["h3"]}
    assert all(x >= 0 for x in global_T.values())

    out = {"status": "passed", "scope": "exact finite n=7 orbit aggregation only; QE8 checked only t=2,3,4 and no all-t/asymptotic claim",
           "source_provenance": {"raw_source_path": "raw/source/short_tail_mod.cpp", "raw_source_sha256": SOURCE_SHA,
             "raw_started_receipt": "raw/started.json", "raw_started_declared_source_sha256": raw_started["source"],
             "historical_binary_hash_receipt": {"sha256": HISTORICAL_BINARY_SHA, "included_in_public_package": False, "verified_during_this_replay": False,
               "note": "The historical executable is intentionally not distributed; this replay verifies its archived source and receipts, not a missing binary."},
             "raw_formula_status_preserved": "pending independent proof", "raw_formula_status_note": "Historical raw label retained; the r29 identity is used without rewriting archived Mac outputs.",
             "r29_verified_module_source_sha256": r29_sources},
           "pratt": {"certificate_path": "pratt/pratt-certificates.json", "certificate_sha256": digest(pratt_path), "roots": [str(p) for p in primes], "verified_nodes": pratt_nodes},
           "crt": {"methods": ["direct cofactor CRT", "incremental Garner CRT"], "modulus": str(Q), "modulus_bits": Q.bit_length(), "count_bound_48_factorial": str(bound), "count_bound_bits": bound.bit_length(), "unique_under_bound": True},
           "orbits": {str(v): {"multiplicity": x["multiplicity"], "A2v": str(x["A2v"]), "Bv": str(x["Bv"]), "h": {k: str(z) for k, z in x["h"].items()}, "T": {str(k): str(z) for k, z in x["T"].items()}, "independent_sets": x.get("independent_sets"), "max_layer_count": x.get("max_layer_count"), "reference_path": x["reference_path"], "reference_sha256": x["reference_sha256"], "center_reused_from_r30": x.get("center_reused_from_r30", False)} for v, x in outputs.items()},
           "global": {"A2": str(global_A), "B": str(global_B), "h": {k: str(z) for k, z in global_h.items()}, "T": {str(k): str(z) for k, z in global_T.items()}, "qe8_t2_to_4": [qe8(global_B, global_T[t], t) for t in (2, 3, 4)]},
           "local_qe8_t2_to_4": local_qe8,
           "local_failures": [{"v": v, "t": row["t"], "slack_rhs_minus_lhs": row["slack_rhs_minus_lhs"]} for v, rows in local_qe8.items() for row in rows if not row["passes"]],
           "verified_r29_layer_histograms": layers, "raw_checkpoints": provenance,
           "center_publicization": {"historical_original_center_path": "historical/short-tail-n7-center-crt-audit-original.json", "historical_original_center_sha256": HISTORICAL_CENTER_SHA, "public_center_path": "center/short-tail-n7-center-crt-audit-public.json", "public_center_sha256": digest(center_path), "io_only_changed_fields": ["pratt.certificate_path", "reference_selected_normalizer.path", "raw_files[*].json/exit/resource"], "normalized_semantic_sha256": canonical_digest(centre_semantics(center)), "note": "The public r030 form changes only documented path strings relative to the historical center JSON; this replay asserts equality after exactly that normalization."}}
    if args.historical_audit:
        out["historical_comparison"] = check_historical(out, args.historical_audit.resolve())
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"status": out["status"], "output": output.name, "sha256": digest(output)}, indent=2))


if __name__ == "__main__":
    main()
