"""Join the frozen whole-leaf witnesses for every selection of the four point sets."""
from collections import defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import sys
from zipfile import ZipFile


HERE = Path(__file__).resolve().parent
R2 = "public-replay-r2-083750/"


def ledger(archive, filename):
    return json.loads(archive.read(filename))


def add_certificate(proofs, leaf, inputs):
    inputs = frozenset(inputs)
    proofs[leaf].add(inputs)
    # The new P30 dominates every old closed support row used by these kernels.
    if "P30" in inputs:
        proofs[leaf].add((inputs - {"P30"}) | {"newP30"})


def summarize():
    with ZipFile(HERE / "wave2.zip") as archive:
        source = archive.read(R2 + "joined-r2-plus-five-index.json")
        assert sha256(source).hexdigest() == "65b4ab38358276aed105897686075377468b4443bd6d74556b272d28df973966"
        old = json.loads(source)
        replacement = ledger(archive, R2 + "new-p30/joined-ab.json")
        scope = {f"{part}/{i}" for part, indices in ledger(archive, R2 + "stage4-retained.json").items() for i in indices}
    proofs = defaultdict(set)
    remaining = set()
    for part, lists in old["parts"].items():
        for i in lists["five"]:
            add_certificate(proofs, f"{part}/{i}", ())
        for i in lists["P30"]:
            add_certificate(proofs, f"{part}/{i}", ("P30",))
        q_alone = set(lists["Q11_without_P30"])
        assert q_alone <= set(lists["Q11"])
        for i in lists["Q11"]:
            add_certificate(proofs, f"{part}/{i}", ("oldQ11",) if i in q_alone else ("oldQ11", "P30"))
        remaining.update(f"{part}/{i}" for i in lists["remaining"])
    assert set(proofs).isdisjoint(remaining) and set(proofs) | remaining == scope
    assert len(remaining) == 873 and len(scope) == 5442
    replacement_only = {leaf for leaf, label in replacement["proof_by_leaf"].items() if label == "continued_new_P30"}
    assert len(replacement_only) == 10 and replacement_only <= remaining
    for leaf in replacement_only:
        add_certificate(proofs, leaf, ("newP30",))
    with ZipFile(HERE / "additional-eleven.zip") as archive:
        additions = ledger(archive, "joint-manifest.json")["leaves"]
    assert len(additions) == 59 and additions.keys() <= remaining
    assert replacement_only.isdisjoint(additions)
    added_nodes = sum(record["suffix"]["visited"] for record in additions.values())
    roots_and_children = sum(record["whole"]["visited"] for record in additions.values())
    for leaf, record in additions.items():
        assert record["whole"]["visited"] == 2*record["whole"]["terminals"] - 1
        inputs = record["whole"]["extra_inputs"]
        assert set(inputs) <= {"P30", "oldQ11", "newQ11"}
        add_certificate(proofs, leaf, inputs)
    models = []
    point_sets = ("P30", "newP30", "oldQ11", "newQ11")
    for size in range(5):
        for family in combinations(point_sets, size):
            available = frozenset(family)
            closed = {leaf for leaf, choices in proofs.items() if any(row <= available for row in choices)}
            pending = sorted(scope - closed, key=lambda s: tuple(map(int, s.split("/"))))
            digest = sha256(("\n".join(pending) + "\n").encode()).hexdigest()
            models.append({"point_sets": list(family), "closed": len(closed), "pending": len(pending),
                           "pending_id_list_sha256": digest})
    return {"baseline_curves": ["disk", "R3", "R5", "R7", "R9"], "threshold": "67/80",
            "scope": len(scope), "joined_roots": len(additions), "joined_gap_visits": added_nodes,
            "joined_total_visits": roots_and_children, "models": models}


if __name__ == "__main__":
    Path(sys.argv[1]).write_text(json.dumps(summarize(), indent=2) + "\n")
