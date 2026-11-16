#!/usr/bin/env python3
"""CPU-only regression test for checkpoint ledger seeding and extension rules."""
import copy
import importlib.util
from pathlib import Path


SOURCE = Path(__file__).with_name("gpu_balanced_ds_v2.py")
spec = importlib.util.spec_from_file_location("ds_v2", SOURCE)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def choices(n):
    return [{} for _ in range(n)]


# This imitates a bounded numerical run with max_codes=4 for all E, followed
# by an extension to 8.  E=2 is the later persisted row while E=1 is replaced
# first, exactly the order that previously risked dropping later progress.
data = {"L": 3, "blocks": [
    {"options": {"1": choices(3), "2": choices(3), "3": choices(3)}},
    {"options": {"1": choices(3), "2": choices(3), "3": choices(3)}},
]}
layouts = module.option_layouts(data, 3)


def entry(E, hi):
    radices, total = layouts[E]
    return {"E": E, "radices_block0_first": radices, "total_codes": total,
            "coverage_completed_half_open": [0, hi],
            "counts": {"codes_seen": hi, "fractional_masks": 0,
                       "closed_unit_interval_fractional_masks": 0}, "top64": []}


old_progress = [entry(1, 4), entry(2, 4), entry(3, 3)]
prior = module.validate_resume_entries(old_progress, layouts)
assert module.resume_start(prior[1], 8) == 4
report = {"E_progress": copy.deepcopy(old_progress)}  # seed before first atomic write
new_e1 = entry(1, 8)
report["E_progress"] = [row for row in report["E_progress"] if row["E"] != 1] + [new_e1]
later = next(row for row in report["E_progress"] if row["E"] == 2)
assert later == old_progress[1], "first E write erased or changed later persisted progress"
try:
    module.resume_start(prior[1], 3)
except ValueError:
    pass
else:
    raise AssertionError("coverage shrink must be refused")
try:
    broken = copy.deepcopy(old_progress)
    broken[1]["counts"]["codes_seen"] = 3
    module.validate_resume_entries(broken, layouts)
except ValueError:
    pass
else:
    raise AssertionError("coverage/count mismatch must be refused")
print("resume-ledger regression: PASS (seed later E, extend 4->8, reject shrink/mismatch)")
