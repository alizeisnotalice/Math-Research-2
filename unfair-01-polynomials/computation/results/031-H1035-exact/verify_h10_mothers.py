"""Mother-only Fraction/Rouche verifier for the H10(3,5) coefficient phase.

The source JSON was made for the older (2,3) task; this checker deliberately
reads no lift records and re-verifies only its ten H10 mother-root disks from
the H10 coefficients themselves.
"""
from fractions import Fraction as F
from exact_complex import cj, pf, mother, dis

H10 = [1,1,1,0,0,0,0,0,0,1,1]


def verify(data: dict) -> dict:
    assert data["format"] == "H10-mother-balls-from-026-v1"
    assert data["H"] == H10 and len(data["source_root_balls_sha256"]) == 64
    assert data["coverage"] == {"upper_half_mother_disks":5,"conjugate_root_disks":10,
                                "old_lifts_intentionally_omitted":True}
    scale = int(data["scale"])
    assert scale > 0 and isinstance(data["mothers"],list) and len(data["mothers"]) == 5
    mothers = sorted(data["mothers"], key=lambda x:x["block"])
    assert [m["block"] for m in mothers] == list(range(5))
    disks = []
    for block,m in enumerate(mothers):
        assert m["id"] == f"mother-{block}" and m["kind"] == "mother"
        z,r = cj(m["disk"]["center"]),pf(m["disk"]["radius"])
        assert r > 0 and z[1] > r
        lhs,rhs = mother(H10,z,r,scale)
        assert lhs < rhs and pf(m["lhs"]) == lhs and pf(m["rhs"]) == rhs
        disks.extend(((z,r),((z[0],-z[1]),r)))
    assert len(disks) == 10
    for i,(z,r) in enumerate(disks):
        for w,s in disks[:i]: assert dis(z,r,w,s)
    return {"status":"H10_10_mother_roots_exact_fraction_rouche_passed",
            "mother_pairs":5,"covered_roots":10,
            "source_scope":"only old mother disks were used; old 2/3 lifts ignored"}
