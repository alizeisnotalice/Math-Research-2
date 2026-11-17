#!/usr/bin/env python3
"""Exact cover-DP checks for the MM affine-completion capacity jump."""
import argparse,hashlib,json,sys,time
from pathlib import Path
from functools import lru_cache


def capacity(code, full):
    def dim_CA(A):
        return sum((c & ~A) == 0 for c in code).bit_length() - 1
    @lru_cache(None)
    def F(R):
        if R == 0:
            return 0
        anchor = (R & -R).bit_length() - 1
        return max(dim_CA(A) + 1 + F(R ^ A)
                   for A in code if A and (A & ~R) == 0 and ((A >> anchor) & 1))
    return F(full)


def main():
    started=time.monotonic()
    parser=argparse.ArgumentParser();parser.add_argument("--out",type=Path,required=True);args=parser.parse_args()
    rows = []
    for q in (4, 8):
        full = (1 << q) - 1
        T = full ^ 1 ^ (1 << 1)
        C = frozenset(mask for mask in range(1 << q) if (bin(mask & T).count("1") & 1) == 0)
        ell = sum(((y >> 0) & 1) << y for y in range(q))
        Cp = frozenset(set(C) | {c ^ ell for c in C})
        rows.append({
            "q": q, "dim_C": len(C).bit_length() - 1,
            "ell_in_C": ell in C,
            "dim_completed": len(Cp).bit_length() - 1,
            "kappa_C": capacity(C, full),
            "kappa_completed": capacity(Cp, full),
        })
    assert rows == [
        {"q": 4, "dim_C": 3, "ell_in_C": False, "dim_completed": 4,
         "kappa_C": 6, "kappa_completed": 8},
        {"q": 8, "dim_C": 7, "ell_in_C": False, "dim_completed": 8,
         "kappa_C": 10, "kappa_completed": 16},
    ]
    result={"status":"passed","rows":rows,"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"python":sys.version,"elapsed_seconds":time.monotonic()-started}
    args.out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
