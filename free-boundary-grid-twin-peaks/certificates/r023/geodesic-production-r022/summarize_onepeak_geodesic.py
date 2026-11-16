#!/usr/bin/env python3
"""Strict Decimal80 aggregation for complete r022 one-peak geodesic JSONL.

This tool deliberately accepts no partial run: the sidecar fixes the source
hash, all metadata parameters, and the deterministic seed formula.  A zero
Double ``Z`` is permitted for a non-dead replicate when its finite ``logZ``
shows that this is only underflow; weights are reconstructed from logZ.
"""
import argparse
import hashlib
import json
from decimal import Decimal, InvalidOperation, getcontext
from pathlib import Path

getcontext().prec = 80
ZERO, ONE = Decimal(0), Decimal(1)
DOUBLE_MIN_SUBNORMAL = Decimal("5e-324")
DOUBLE_ROUNDING_ABS = DOUBLE_MIN_SUBNORMAL / Decimal(2)
DOUBLE_RELATIVE_TOLERANCE = Decimal("2e-12")
REQUIRED_META = {
    "implementation", "base_source_sha256", "n", "particles", "replicates",
    "tau_lo", "tau_hi", "domain", "method", "move_steps",
}


class ValidationError(ValueError):
    pass


def fail(message):
    raise ValidationError(message)


def decimal(value, label):
    if value is None:
        fail(f"{label}: expected a finite decimal, got null")
    try:
        out = value if isinstance(value, Decimal) else Decimal(str(value))
    except (InvalidOperation, ValueError) as exc:
        raise ValidationError(f"{label}: malformed decimal {value!r}") from exc
    if not out.is_finite():
        fail(f"{label}: non-finite decimal")
    return out


def close(a, b, label, relative=DOUBLE_RELATIVE_TOLERANCE):
    """Double-compatible comparison, including subnormal rounding allowance."""
    scale = max(abs(a), abs(b))
    allowance = DOUBLE_ROUNDING_ABS + relative * scale
    if abs(a - b) > allowance:
        fail(f"{label}: {a} != {b}; tolerance is {allowance}")


def strict_integer(value, label):
    """Counts must be JSON integers, never bools, strings, or Decimal 1.5."""
    if type(value) is not int:
        fail(f"{label}: expected a JSON integer")
    return value


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_jsonl(path):
    records = []
    for number, line in enumerate(Path(path).read_text().splitlines(), start=1):
        if not line.strip():
            fail(f"line {number}: blank line is not permitted")
        try:
            item = json.loads(line, parse_float=Decimal, parse_int=int)
        except json.JSONDecodeError as exc:
            raise ValidationError(f"line {number}: invalid JSON") from exc
        if not isinstance(item, dict) or "kind" not in item:
            fail(f"line {number}: object with kind is required")
        records.append(item)
    if not records:
        fail("JSONL is empty")
    return records


def load_sidecar(path, source_path):
    sidecar = json.loads(Path(path).read_text(), parse_float=Decimal, parse_int=int)
    if set(sidecar) != {"source_sha256", "metadata", "seed_formula"}:
        fail("sidecar must contain exactly source_sha256, metadata, seed_formula")
    source_hash = sidecar["source_sha256"]
    if not isinstance(source_hash, str) or len(source_hash) != 64:
        fail("sidecar source_sha256 must be 64 hexadecimal characters")
    if sha256(source_path) != source_hash:
        fail("source SHA-256 does not match sidecar")
    expected = sidecar["metadata"]
    if not isinstance(expected, dict) or not REQUIRED_META.issubset(expected):
        fail("sidecar metadata omits required parameters")
    formula = sidecar["seed_formula"]
    if set(formula) != {"base_hex", "n_multiplier", "tau_multiplier"}:
        fail("sidecar seed_formula must specify base_hex, n_multiplier, tau_multiplier")
    try:
        int(formula["base_hex"], 16)
    except (TypeError, ValueError) as exc:
        raise ValidationError("seed_formula base_hex is invalid") from exc
    if not all(isinstance(formula[key], int) and formula[key] > 0
               for key in ("n_multiplier", "tau_multiplier")):
        fail("seed formula multipliers must be positive integers")
    return sidecar


def expected_seed(meta, formula, replicate):
    base = int(formula["base_hex"], 16)
    domain = int(str(meta["domain"]), 16)
    raw = (base + domain + meta["n"] * formula["n_multiplier"]
           + meta["tau_lo"] * formula["tau_multiplier"] + replicate)
    return raw & ((1 << 64) - 1)


def validate_inverse(record, particles, groups, where):
    for name in ("inverse_m2_groups", "inverse_m3_groups"):
        values = record.get(name)
        if not isinstance(values, list) or len(values) != groups:
            fail(f"{where}: {name} must have {groups} entries")
        integers = [strict_integer(x, f"{where}.{name}") for x in values]
        if any(x < 0 for x in integers):
            fail(f"{where}: {name} has a negative count")
        mean_name = "mean_m2" if name.endswith("m2_groups") else "mean_m3"
        scaled_name = "scaled_mean_m2" if name.endswith("m2_groups") else "scaled_mean_m3"
        zmean_name = "Zmean_m2" if name.endswith("m2_groups") else "Zmean_m3"
        mean = decimal(record.get(mean_name), f"{where}.{mean_name}")
        close(mean, Decimal(sum(integers)) / Decimal(particles), f"{where}.{mean_name}")
        close(decimal(record.get(scaled_name), f"{where}.{scaled_name}"), mean,
              f"{where}.{scaled_name}")
        z = decimal(record.get("Z"), f"{where}.Z")
        close(decimal(record.get(zmean_name), f"{where}.{zmean_name}"), z * mean,
              f"{where}.{zmean_name}")


def validate_replicate(record, meta, formula, seen):
    n, particles, reps = meta["n"], meta["particles"], meta["replicates"]
    groups, m = (particles + 255) // 256, 2 * n - 1
    where = f"replicate {record.get('replicate')!r}"
    if record.get("tau") != n * n + 1:
        fail(f"{where}: tau must equal N+1")
    r = record.get("replicate")
    if not isinstance(r, int) or not 0 <= r < reps or r in seen:
        fail(f"{where}: duplicate or out-of-domain replicate index")
    seen.add(r)
    actual_seed = record.get("seed_hex")
    try:
        actual_seed_value = int(actual_seed, 16) if isinstance(actual_seed, str) else None
    except ValueError:
        actual_seed_value = None
    if actual_seed_value != expected_seed(meta, formula, r):
        fail(f"{where}: seed is outside the declared deterministic domain")
    dead = record.get("dead")
    if not isinstance(dead, bool):
        fail(f"{where}: dead must be boolean")
    geo_groups, scaled = record.get("geodesic_prefix_group_counts"), record.get("scaled_geodesic_prefix_p")
    if dead:
        # On a connected grid, a one-birth terminal word has a boundary at
        # every unfinished prefix.  Therefore a onepeak SMC run cannot die.
        fail(f"{where}: dead=true violates the connected-grid onepeak boundary invariant")
    logz = decimal(record.get("logZ"), f"{where}.logZ")
    # A non-dead Double underflow has Z=0 but finite logZ.  Accept exactly the
    # range where exp(logZ) rounds below half the smallest Double subnormal;
    # otherwise verify Z against logZ with relative plus subnormal tolerance.
    raw_z = decimal(record.get("Z"), f"{where}.Z")
    if raw_z < ZERO:
        fail(f"{where}.Z: negative raw Z")
    reconstructed_z = logz.exp()
    if raw_z == ZERO:
        if reconstructed_z >= DOUBLE_ROUNDING_ABS:
            fail(f"{where}.Z=0 is not justified by finite logZ underflow")
    else:
        close(raw_z, reconstructed_z, f"{where}.Z versus exp(logZ)")
    if not isinstance(geo_groups, list) or len(geo_groups) != groups:
        fail(f"{where}: geodesic group count has wrong group dimension")
    raw = [0] * m
    for group, values in enumerate(geo_groups):
        cap = min(256, particles - 256 * group)
        if not isinstance(values, list) or len(values) != m:
            fail(f"{where}: group {group} has wrong d dimension")
        values = [strict_integer(x, f"{where}.geodesic group {group}") for x in values]
        if any(x < 0 or x > cap for x in values):
            fail(f"{where}: group {group} is outside [0,{cap}]")
        if any(values[k + 1] > values[k] for k in range(m - 1)):
            fail(f"{where}: group {group} is not monotone in d")
        for k, value in enumerate(values):
            raw[k] += value
    if raw[0] != particles or raw[1] != particles:
        fail(f"{where}: non-dead onepeak run requires d=1,d=2 equal P")
    if any(raw[k + 1] > raw[k] for k in range(m - 1)):
        fail(f"{where}: aggregate geodesic counts are not monotone")
    if not isinstance(scaled, list) or len(scaled) != m:
        fail(f"{where}: scaled_geodesic_prefix_p has wrong dimension")
    reported = [decimal(x, f"{where}.scaled[{k}]") for k, x in enumerate(scaled)]
    # Use the exact integer numerator/P in aggregation, after checking that
    # the GPU's printed Double agrees.  This avoids importing print-rounding
    # error into the Decimal80 weighted ratios.
    points = [Decimal(value) / Decimal(particles) for value in raw]
    for k, value in enumerate(reported):
        close(value, points[k], f"{where}.scaled[{k}]")
    validate_inverse(record, particles, groups, where)
    return {"replicate": r, "logZ": logz, "p": points, "raw": raw}


def diagnostic_weighted_summary(valid, m):
    if not valid:
        return [{"d": d + 1, "observed": False, "reason": "all_replicates_dead"}
                for d in range(m)]
    max_logz = max(row["logZ"] for row in valid)
    weights = [(row["logZ"] - max_logz).exp() for row in valid]
    total_weight = sum(weights, ZERO)
    alpha = [weight / total_weight for weight in weights]
    R = len(valid)
    out = []
    means = []
    for k in range(m):
        values = [row["p"][k] for row in valid]
        mean = sum((a * value for a, value in zip(alpha, values)), ZERO)
        nonzero = sum(value != ZERO for value in values)
        event_weights = [weight * value for weight, value in zip(weights, values)]
        event_total = sum(event_weights, ZERO)
        item = {"d": k + 1, "nonzero_replicates": nonzero}
        if event_total == ZERO:
            # A literal weighted zero after no observed events is not an
            # estimate of a true zero probability, so do not emit it as one.
            item.update({"observed": False, "weighted_ratio_decimal": None, "estimate": None,
                         "weighted_event_ess": None, "weighted_event_maxshare": None,
                         "diagnostic_se_decimal": None})
            means.append(None)
        else:
            shares = [value / event_total for value in event_weights]
            ess = ONE / sum((share * share for share in shares), ZERO)
            maxshare = max(shares)
            if R > 1:
                # Ratio-of-means delta diagnostic.  This is deliberately not
                # a Kish-ESS weighted-sample-variance approximation.
                variance = (Decimal(R) / Decimal(R - 1)
                            * sum(((a * (value - mean)) ** 2
                                   for a, value in zip(alpha, values)), ZERO))
                se = variance.sqrt()
            else:
                variance, se = None, None
            item.update({"observed": True, "weighted_ratio_decimal": format(mean, "f"),
                         "estimate": format(mean, "f"),
                         "weighted_event_ess": format(ess, "f"),
                         "weighted_event_maxshare": format(maxshare, "f"),
                         "diagnostic_se_decimal": None if se is None else format(se, "f")})
            means.append(mean)
        out.append(item)
    for k in range(m - 1):
        item = out[k]
        if means[k] is None or means[k + 1] is None or means[k] == ZERO:
            item.update({"qbar_next_over_current": None, "qbar_delta_se_decimal": None})
            continue
        qbar = means[k + 1] / means[k]
        delta = None
        if R > 1:
            # Direct delta residual for qbar=p_{d+1}/p_d:
            # R/(R-1) sum [alpha (p_next-qbar*p_current)]^2 / pbar_d^2.
            variance = (Decimal(R) / Decimal(R - 1)
                        * sum(((a * (row["p"][k + 1] - qbar * row["p"][k])) ** 2
                               for a, row in zip(alpha, valid)), ZERO)
                        / (means[k] ** 2))
            delta = variance.sqrt()
        item.update({"qbar_next_over_current": format(qbar, "f"),
                     "qbar_delta_se_decimal": None if delta is None else format(delta, "f")})
    out[-1].update({"qbar_next_over_current": None, "qbar_delta_se_decimal": None})
    return out


def summarize(input_path, sidecar_path, source_path):
    sidecar = load_sidecar(sidecar_path, source_path)
    expected, formula = sidecar["metadata"], sidecar["seed_formula"]
    records = load_jsonl(input_path)
    metadata = [row for row in records if row["kind"] == "metadata"]
    endings = [row for row in records if row["kind"] == "run_end"]
    replicates = [row for row in records if row["kind"] == "replicate"]
    if len(metadata) != 1 or len(endings) != 1 or records[-1]["kind"] != "run_end":
        fail("exactly one metadata first and one final run_end are required")
    if records[0]["kind"] != "metadata" or len(records) != len(replicates) + 2:
        fail("JSONL contains unknown records or records in the wrong position")
    meta = metadata[0]
    if any(meta.get(key) != value for key, value in expected.items()):
        fail("actual metadata does not match explicit sidecar metadata")
    if meta["tau_lo"] != meta["tau_hi"] or meta["tau_lo"] != meta["n"] ** 2 + 1:
        fail("only complete onepeak tau=N+1 inputs are accepted")
    if len(replicates) != meta["replicates"]:
        fail("replicate count does not equal metadata R")
    seen, valid = set(), []
    for row in replicates:
        checked = validate_replicate(row, meta, formula, seen)
        if checked is not None:
            valid.append(checked)
    if seen != set(range(meta["replicates"])):
        fail("missing replicate index in complete R domain")
    m = 2 * meta["n"] - 1
    return {
        "kind": "strict-onepeak-geodesic-summary",
        "input": str(Path(input_path)),
        "source_sha256": sidecar["source_sha256"],
        "metadata": meta,
        "replicates_total": meta["replicates"],
        "replicates_nondead": len(valid),
        "replicates_dead": meta["replicates"] - len(valid),
        "decimal_precision": getcontext().prec,
        "independence_note": "SE and delta-SE are independent-replicate diagnostics only; no mixing claim is made.",
        "d_summary": diagnostic_weighted_summary(valid, m),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="complete JSONL run")
    parser.add_argument("--sidecar", required=True, help="explicit source/metadata/seed JSON")
    parser.add_argument("--source", required=True, help="actual source file to SHA-256 verify")
    parser.add_argument("--out", required=True, help="summary JSON path")
    args = parser.parse_args()
    try:
        result = summarize(args.input, args.sidecar, args.source)
    except ValidationError as exc:
        raise SystemExit(f"VALIDATION FAILED: {exc}")
    Path(args.out).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"passed": True, "out": args.out, "nondead": result["replicates_nondead"]}))


if __name__ == "__main__":
    main()
