#!/usr/bin/env python3
"""Read-only finite plug-in consistency audit for completed tau=2 distance runs."""
from decimal import Decimal, localcontext
from fractions import Fraction
import gzip, hashlib, json, math, os, pathlib, sys

OUT = pathlib.Path(__file__).resolve().parent
CERT_ROOT = pathlib.Path(os.environ.get('A022_CERTIFICATES_ROOT', OUT.parent.parent))
SOURCE = CERT_ROOT / "r024" / "r024-local-distance" / "source" / "tau2_distance_terminal_r023.swift"
CONFIGS = [
  ("n8_macbook_P16384_R32", CERT_ROOT / "r024" / "r024-macbook" / "n8" / "n8-distance-p16384-r32-l16.jsonl", CERT_ROOT / "r024" / "r024-macbook" / "n8" / "n8-distance-p16384-r32-l16.sidecar.json"),
  ("n12_macbook_P131072_R128", CERT_ROOT / "r024" / "r024-macbook" / "n12" / "n12-distance-p131072-r128-l16.jsonl", CERT_ROOT / "r024" / "r024-macbook" / "n12" / "n12-distance-p131072-r128-l16.sidecar.json"),
  ("n16_local_P131072_R128", CERT_ROOT / "r024" / "r024-local-distance" / "raw" / "n16-p131072-r128-l16.jsonl", CERT_ROOT / "r024" / "r024-local-distance" / "receipts" / "n16-distance.sidecar.json"),
  ("n16_macbook_P262144_R128", CERT_ROOT / "r025" / "r025-macbook-gpu" / "n16-distance-p262144-r128-l16.jsonl.gz", CERT_ROOT / "r025" / "r025-macbook-gpu" / "n16-distance-p262144-r128-l16.sidecar.json"),
]

def sha(p):
    h=hashlib.sha256()
    opener=gzip.open if str(p).endswith('.gz') else open
    with opener(p,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()
def D(x): return Decimal(str(x))
def dstr(x): return format(x, ".70g")
def fstr(x): return dstr(Decimal(x.numerator)/Decimal(x.denominator))
def load(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        return [json.loads(line, parse_float=Decimal) for line in f if line.strip()]

def audit(name, raw_path, side_path):
    lines=load(raw_path); side=json.loads(pathlib.Path(side_path).read_text())
    meta=[x for x in lines if x.get("kind")=="metadata"]
    end=[x for x in lines if x.get("kind")=="run_end"]
    reps=[x for x in lines if x.get("kind")=="replicate" and x.get("tau")==2]
    if len(meta)!=1 or len(end)!=1: raise ValueError(f"{name}: metadata/run_end count")
    m=meta[0]; expected=side["expected_metadata"]
    if sha(raw_path)!=side["raw_sha256"]: raise ValueError(f"{name}: raw SHA mismatch")
    for key,value in expected.items():
        if m.get(key)!=value: raise ValueError(f"{name}: metadata mismatch {key}")
    if side["source_sha256"]!=sha(SOURCE): raise ValueError(f"{name}: source SHA mismatch")
    R=int(m["replicates"]); n=int(m["n"])
    if sorted(int(x["replicate"]) for x in reps)!=list(range(R)): raise ValueError(f"{name}: replicate IDs")
    seeds=[x["seed_hex"] for x in reps]
    if len(set(seeds))!=R: raise ValueError(f"{name}: duplicate seeds")
    if any(x.get("dead") is not False for x in reps): raise ValueError(f"{name}: dead replicate")
    bins=2*n-1
    with localcontext() as ctx:
        ctx.prec=100
        maxlog=max(x["logZ"] for x in reps)
        weighted_Z=[]; weighted_L=[]; weighted_hist=[Fraction(0) for _ in range(bins)]
        group_checks=True
        for x in reps:
            z=Fraction((x["logZ"]-maxlog).exp())
            Lgroups=x["tau2_lift_L_groups"]
            Hgroups=x["tau2_distance_L_hist_groups"]
            if len(Lgroups)!=len(Hgroups) or any(len(h)!=bins for h in Hgroups): raise ValueError(f"{name}: group shape")
            if any(sum(int(v) for v in h)!=int(l) for h,l in zip(Hgroups,Lgroups)): group_checks=False
            # L_groups stores the raw sum over `particles` terminal samples.
            # The plug-in B/A uses the per-particle mean L; the distance
            # numerator and denominator retain raw units and cancel together.
            L=sum(int(v) for v in Lgroups)
            weighted_Z.append(z); weighted_L.append(z*Fraction(L, int(m["particles"])))
            for d in range(bins): weighted_hist[d]+=z*sum(int(h[d]) for h in Hgroups)
        if not group_checks: raise ValueError(f"{name}: hist/L group mismatch")
        denom_mean_L=sum(weighted_L); ratio=sum(weighted_Z)/denom_mean_L
        if weighted_hist[0] or weighted_hist[1]: raise ValueError(f"{name}: nonzero d=0/1 mass")
        denom_raw_L=sum(weighted_hist)
        if not denom_raw_L: raise ValueError(f"{name}: zero distance denominator")
        p=[x/denom_raw_L for x in weighted_hist]
        g=sum(Fraction(d,2*(d-1))*p[d] for d in range(2,bins))
        slack=g-ratio
        defect=max(Fraction(0),-slack)
        bound=defect/(Fraction(2)-g) if defect else Fraction(0)
        return {"name":name,"raw_certificate_path":str(pathlib.Path(raw_path).relative_to(CERT_ROOT)),"raw_sha256":sha(raw_path),"sidecar_certificate_path":str(pathlib.Path(side_path).relative_to(CERT_ROOT)),
          "source_file":SOURCE.name,"source_sha256":sha(SOURCE),"metadata":expected,"run_end":str(end[0]["wall_seconds"]),
          "strict": {"metadata_one":True,"run_end_one":True,"raw_sha_matches_sidecar":True,"source_sha_matches_sidecar":True,"replicate_ids_complete":True,"seeds_unique":True,"all_dead_false":True,"distance_group_L_exact":True,"distance_global_L_exact":True,"distance_0_1_zero":True},
          "fraction_backend":"Decimal(100 digits) exp(logZ-maxlog), then Decimal-to-Fraction; all comparisons on Fraction",
          "B_over_A_plugin":fstr(ratio),"distance_rhs_gD_plugin":fstr(g),"slack_gD_minus_ratio":fstr(slack),
          "negative_defect":fstr(defect),"if_negative_max_abs_Bratio_or_TV_lower_bound":fstr(bound),
          "distance_support":[d for d,v in enumerate(p) if v],"P2_mass_exact_one":sum(p)==1}

def main():
    seen=set(); out=[]
    for cfg in CONFIGS:
        r=audit(*cfg)
        if r["raw_sha256"] in seen: raise ValueError("duplicate raw SHA")
        seen.add(r["raw_sha256"]); out.append(r)
    result={"scope":"separate finite plug-in checks; batches are never pooled and no significance or bug claim follows", "formula":"B/A2 <= g_D=sum_{d>=2} d/(2(d-1)) P2(d); B/A2 plugin=sum Z/sum Z*mean(L), mean(L)=raw_Lsum/particles; P2 plugin=sum Z*raw_Lhist/sum Z*raw_Lsum", "reproduce":"python3 audit_scale.py > scale-audit.json", "batches":out}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=="__main__": main()
