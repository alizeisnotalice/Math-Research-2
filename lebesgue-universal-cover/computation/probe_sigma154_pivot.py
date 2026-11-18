"""Verify physical extent and exact area after enlarging the AH reference grid."""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
import io
import json
from pathlib import Path
import subprocess
import sys
from contextlib import redirect_stdout

from flint import ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from certify_three_sigma_area_budget import ah_clearance, enclosure, fixed_ah, geometries, pack, setup, subtract


def run(comp, config):
    sys.path.insert(0, str(comp))
    import audit_optimized_eh_chord as g
    from certify_ah_margin_envelope import replace_line
    from certify_sigma_extended_graph import geometry

    template, initial, _ = geometry(comp, [(F(".429"), F(".010"))], sigma="1.54",
                                     far=F(".292"), negative=F(".329"), eps=F(".0114"))
    frozen = json.loads((HERE / "sigma154-ah429-eh-chain997-composition.json").read_bytes())
    batches = {}
    rows = []
    for sigma, ps, ds, fs, ns, es, end in config:
        p, turn, far, negative, eps, domain = map(F, (ps, ds, fs, ns, es, end))
        assert F(".429") <= p < domain <= F(".48") and int(domain*1000) == domain*1000
        assert negative*5000 == int(negative*5000) and p*5000 == int(p*5000)
        assert far*1000 == int(far*1000)
        code = replace_line(initial, 'sigma_degrees = "1.54"', f'sigma_degrees = "{sigma}"')
        code = replace_line(code, 'F("429/1000"), F("1/100")', f'F("{p}"), F("{turn}")')
        code = replace_line(code, 'F("429/1000"), F("73/250"),F("12/25"), F("57/5000")',
                            f'F("{p}"), F("{far}"),F("12/25"), F("{eps}")')
        code = replace_line(code, 'pbound = F("329/1000")', f'pbound = F("{negative}")')
        code = replace_line(code, 'range(401,431)', f'range(401,{int(domain*1000)+1})')
        chart = {"__file__": template["__file__"], "__name__": "__main__"}
        with redirect_stdout(io.StringIO()):
            exec(compile(code, "pivot_release_geometry", "exec"), chart)
        A = chart["A"]
        assert chart["knots"][-1] == domain and chart["cuts"][-1] == (p,turn)
        assert all(v<p and d<turn for v,d in chart["cuts"][:-1])
        assert A(negative)>chart["physical_rights"][-1]

        def reference_gap(w):
            left = chart["chart"](chart["add"](chart["z"], chart["mul"](A(w),chart["D"])))[0]
            right = chart["coord_hull"](chart["psi"], left)
            return min(min((l-left).lower(), (right-l).lower())
                       for pair in chart["parameter_boxes"] for l in pair[:2])

        enlarged_gap = reference_gap(domain)
        assert enlarged_gap>0
        old_gap = reference_gap(F(".43"))
        max_deficit = max(chart["deficit_by_bin"]).upper()
        ctx.prec=160
        if sigma not in batches:
            data = json.loads(subprocess.check_output([sys.executable, "-E", "-B",
                        str(comp/"certify_gibbs_area_refined.py"), sigma, "16", "64"]))
            assert F(data["sigma_degrees"])==F(sigma) and data["passed"]
            base=tuple(map(F,data["area_interval"]))
            e = geometries(sigma, frozen["old_parameters"], frozen["EH_chain_params"], g)
            assert sorted(e["facets"])==list(range(11))
            eh = enclosure(e["new_eh"]-e["old_eh"])
            batches[sigma]=base,eh
        base,eh=batches[sigma]
        cap, cross, derivative = fixed_ah(comp, sigma, p, turn, far)
        assert cross[1]<far and derivative>0
        geo = setup(sigma, g)
        U,C,T=(geo[k] for k in ("U","C3","T"))
        radial=g.minus(U,C)
        c,s=g.unit(g.arb(frozen["EH_deepest"]))
        V=g.plus(C,g.minus(g.scale(c,radial),g.scale(s,g.rot(radial))))
        clear=ah_clearance(geo,V,T,g,p,turn)
        rows.append(dict(sigma=sigma,pivot=ps,turn=ds,reference_end=end,
            geometry_sha256=sha256(code.encode()).hexdigest(),reference_gap=pack(enlarged_gap),
            old_043_reference_gap=pack(old_gap), physical_rows=len(chart["physical_deficits"]),
            upper_deficit=pack(max_deficit), AH_crossing=pack(cross),
            base_area=pack(base,35),AH_area=pack(cap,35),EH_area=pack(eh,35),
            AH_EH_clearance=pack(clear),candidate_area=pack(subtract(subtract(base,cap),eh),35)))
    return dict(scope="Physical gate and true area. Separate EH and AH placement proofs are required at each angle.",
                parameters=config,rows=rows, source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == "__main__":
    assert __debug__
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("comp",type=Path)
    ap.add_argument("config",type=Path)
    a=ap.parse_args()
    print(json.dumps(run(a.comp.resolve(),json.loads(a.config.read_bytes())),separators=(",",":")))
