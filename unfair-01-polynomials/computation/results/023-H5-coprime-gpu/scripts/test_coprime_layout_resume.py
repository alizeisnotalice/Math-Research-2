#!/usr/bin/env python3
"""CPU-only schema, layout, explicit-mask, and resume-ledger checks."""
import copy
import importlib.util
from pathlib import Path


SOURCE = Path(__file__).with_name("stream_coprime_ds_gpu_v1.py")
spec = importlib.util.spec_from_file_location("coprime_engine", SOURCE)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def option(mode, indices):
    return {"mode": mode, "selected_group_indices": indices}


def copies(count, block, mode="interior"):
    return [option(mode, [100 * block + i]) for i in range(count)]


# a=2,b=3 is a compact exact instance.  Feasible canonical E keys are
# 2,3,4,5 and the terminal E=6; no literal E=1 or E=4 smoke assumption is
# permitted by the engine.
blocks = []
for block in range(3):
    blocks.append({"block": block, "options": {
        "2": copies(3, block),       # C(2,0) C(3,1)
        "3": copies(2, block),       # C(2,1) C(3,0)
        "4": copies(3, block),       # C(2,0) C(3,2)
        "5": copies(6, block),       # C(2,1) C(3,1)
        "6": [option("all_a", [100 * block]), option("all_b", [100 * block + 1])],
    }})
data = {"a": 2, "b": 3, "blocks": blocks}
layouts = module.option_layouts(data, 2, 3)
assert list(sorted(layouts)) == [2, 3, 4, 5, 6]
assert layouts[2] == ([3, 3, 3], 27)
assert layouts[5] == ([6, 6, 6], 216)
assert layouts[6] == ([1, 2, 2], 4)
assert module.verify_coprime_layout_count(layouts, 2, 3) == 282  # F_2 F_3 / 2 + 2
assert module.choose_options(data, 6, 2, 3)[0][0]["mode"] == "all_a"
mask, records = module.selected_groups(blocks, 6, [0, 1, 0])
assert [r["mode"] for r in records] == ["all_a", "all_b", "all_a"]
assert mask == (1 << 0) + (1 << 101) + (1 << 200)


def entry(E, hi):
    radices, total = layouts[E]
    return {"E": E, "radices_block0_first": radices, "total_codes": total,
            "coverage_completed_half_open": [0, hi],
            "counts": {"codes_seen": hi, "fractional_masks": 0,
                       "closed_unit_interval_fractional_masks": 0}, "top64": []}


old = [entry(2, 8), entry(3, 8), entry(4, 8)]
prior = module.validate_resume_entries(old, layouts)
assert module.resume_start(prior[2], 12) == 8
report = {"E_progress": copy.deepcopy(old)}
report["E_progress"] = [row for row in report["E_progress"] if row["E"] != 2] + [entry(2, 12)]
assert next(row for row in report["E_progress"] if row["E"] == 3) == old[1]
try:
    module.resume_start(prior[2], 7)
except ValueError:
    pass
else:
    raise AssertionError("coprime resume must refuse coverage shrink")

h = [1, 1, 0, 1]
assert module.exact_coprime_target(h, 2, 3) == [1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 0, 0, 0, 1]
print("coprime CPU layout/resume test: PASS (actual E keys, terminal all_a, Franel total, explicit mask)")
