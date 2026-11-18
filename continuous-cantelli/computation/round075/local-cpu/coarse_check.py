"""Independent direct integration of the exact coarse-gap decomposition."""
import json
from pathlib import Path
import mpmath as mp


def expected(case,row):
    b,A,alpha,c=(case[k] for k in ("base","digits","alpha","c"))
    m=row["exponent"]
    n=int(row["n"])
    h=(2*mp.pi**2*n*n*mp.mpf(c))**(-1/mp.mpf(alpha))
    power=mp.mpf(1)/alpha
    wave=2*mp.pi*n*h
    rs=[]
    for ell in range(m):
        for d,e in zip(A,A[1:]):
            if e<=d+1:
                continue
            width=mp.mpf(e-d-1)/b**(ell+1)
            cap=min(mp.mpf(130),(width/(2*h))**alpha)
            cuts=[0]+[j for j in (1,4,12,24,48,75) if j<cap]+[cap]
            f=lambda x:x**(power-1)*mp.exp(-x)*mp.cos(wave*x**power) if x else mp.mpf(0)
            rs.append(2*len(A)**ell*h*power*mp.quad(f,cuts))
    return mp.fsum(rs)


def main():
    mp.mp.dps=72
    data=json.loads(Path("suppressed-resonance-results.json").read_text())
    all_cases=[c for shard in data["shards"] for c in shard["tests"]]
    cases={}
    for case in sorted(all_cases,key=lambda q:q["profile"]):
        if case["alpha"]==.55 or (case["base"]==9 and case["alpha"] in (1.95,2.)):
            cases.setdefault((case["base"],case["alpha"]),case)
    records=[]
    for (base,alpha),case in sorted(cases.items()):
        r=next(r for r in reversed(case["rows"]) if r["base_harmonic"]==1)
        target=expected(case,r)
        got=mp.mpc(*[mp.mpf(x) for x in r["coefficient"]])
        cross=abs(target-got)
        rel=cross/max(mp.mpf("1e-50"),abs(target))
        assert cross<mp.mpf("2e-11")
        records.append({"base":base,"alpha":alpha,"n":r["n"],
                        "independent":mp.nstr(target,35),
                        "relative_difference":mp.nstr(rel,15),"abs_difference":mp.nstr(cross,15)})
    Path("coarse-check.json").write_text(json.dumps({"dps":mp.mp.dps,"rows":records},indent=2)+"\n")


if __name__=="__main__":
    main()
