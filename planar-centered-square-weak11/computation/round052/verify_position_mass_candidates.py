#!/usr/bin/env python3
"""Independently sum each retained R52 bbox from sparse source cells, no prefix sum."""
import argparse,hashlib,json,os,platform,sys,time
from fractions import Fraction as F
from pathlib import Path

def main():
    a=argparse.ArgumentParser();a.add_argument("slab");a.add_argument("source_slab");a.add_argument("mass_velocity");a.add_argument("--output",required=True);v=a.parse_args()
    s=time.perf_counter();paths=[Path(v.slab),Path(v.source_slab),Path(v.mass_velocity)];m,source,mu=[json.loads(p.read_text()) for p in paths];mu=list(map(F,mu));xs={(F(x),int(v)):i for i,(x,v) in enumerate(source["x_classes"])};ys={F(y):i for i,y in enumerate(source["y_classes"])}
    # A dictionary of occupied cells, unlike the producer's two prefix sums.
    cells={}
    for q in source["source_catalog"]:
        key=(xs[F(q["x"]),q["v"]],ys[F(q["y"])])
        w,d=cells.get(key,(F(0),F(0)));cells[key]=(w+F(q["w"]),d+mu[q["label"]])
    rows={}
    for (i,j),wd in cells.items():rows.setdefault(i,[]).append((j,*wd))
    count=nonzero=0;lam=F(m["lambda"]);rad=F(m["R"]);T=F(m["dt_radius"])
    for box,r in zip(m["boxes"],m["box_records"],strict=True):
        i,k,j,l=r["bbox"];w=d=F(0)
        for a in range(i,k+1):
            for b,ww,dd in rows.get(a,[]):
                if j<=b<=l:w+=ww;d+=dd
        assert w==F(r["m0"]) and d==F(r["mu"])
        q0,q1=F(r["q0"]),F(r["q1"])
        assert q0==w/(4*lam) and q1==d/(4*lam)
        assert F(r["q_full_domain_lower"])==q0-abs(q1)*T>0
        assert F(r["q_full_domain_upper"])==q0+abs(q1)*T
        assert all(F(e[3])==q0 and F(e[4])==q1 for e in box if e[2])
        assert r["radius_branch"]=="root" and q0+abs(q1)*T<rad*rad
        count+=1;nonzero+=d!=0
    assert count==m["counts"]["candidate_boxes"]==len(m["boxes"])
    assert sum(value for key,value in m["counts"].items() if key.startswith("reject_"))+count==m["counts"]["boxes_visited"]
    out={"status":"passed","method":"independent sparse occupied-cell direct sums, no prefix sum","candidate_mass_and_derivative_checks":count,"nonzero_derivative_candidates":nonzero,"empty_ledger_count":m["counts"]["boxes_visited"]-count,"scope":"every retained actual midpoint bbox; omission inequalities separately audited","python":sys.version,"platform":platform.platform(),"elapsed_seconds":time.perf_counter()-s,"input_sha256":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},"code_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    assert os.path.ismount("/Volumes/data");Path(v.output).write_text(json.dumps(out,indent=2)+"\n");print(json.dumps(out,indent=2))
if __name__=="__main__":main()
