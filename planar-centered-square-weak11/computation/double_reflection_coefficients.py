#!/usr/bin/env python3
"""R54 exact true-axis left-seam and full double-closure coefficient ledger."""
import argparse, hashlib, json, os, platform, sys
from fractions import Fraction as F
from pathlib import Path


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def interval(z):
    return {"lower": str(z[0]), "upper": str(z[1]),
            "lower_float_display": float(z[0]), "upper_float_display": float(z[1])}


def main():
    ap=argparse.ArgumentParser()
    for k in ["original", "reflected", "prior", "axis", "old", "new", "axis-error", "output"]:
        ap.add_argument("--"+k, required=True)
    a=ap.parse_args()
    assert os.path.ismount("/Volumes/data")
    paths={k:Path(getattr(a,k.replace("-","_"))) for k in ["original","reflected","prior","axis","old","new","axis-error"]}
    d={k:json.loads(p.read_text()) for k,p in paths.items()}
    original,reflected=d["original"],d["reflected"]
    oh,rh=sha(paths["original"]),sha(paths["reflected"])
    X=max(F(p[0]) for p in original["pts"])
    assert len(original["pts"])==len(reflected["pts"])
    assert all((X-F(x),F(y),F(m))==tuple(map(F,q)) for (x,y,m),q in zip(original["pts"],reflected["pts"]))
    assert all(0<=F(x)<F(original["Px"]) for x,y,m in reflected["pts"])
    assert all(original[k]==reflected[k] for k in ["Px","Py","lam"])
    assert F(original["lam"])==1 and d["prior"]["input_sha256"]==oh
    for k in ["axis","old","new","axis-error"]:
        assert d[k]["status"]=="PASS" and d[k]["input_sha256"]==rh
    assert d["old"]["state"]=="seam-old" and d["new"]["state"]=="seam-new"
    axh=sha(paths["axis"])
    assert all(d[k]["seam_sha256"]==axh for k in ["old","new","axis-error"])
    axis=d["axis"]
    lo,hi=F(axis["c_edge_lower"]),F(axis["c_edge_upper"])
    assert hi-lo==F(1,1<<64) and lo>X
    for k in ["old","new"]:
        assert F(d[k]["rational_axis"])==hi
        assert list(map(F,d[k]["window"]))==[hi-15,hi,F(0),F(original["Py"])]
        assert F(d[k]["tail_bound"])<1
    states={row["state"]:row["area"] for row in d["prior"]["sources"]}
    A=F(states["bulk"]["lower"]),F(states["bulk"]["upper"])
    bx=(F(states["left"]["lower"])+F(states["right"]["lower"])-4*A[1],
        F(states["left"]["upper"])+F(states["right"]["upper"])-4*A[0])
    glrat=(F(d["new"]["inner_area"])-F(d["old"]["outer_area"]),
           F(d["new"]["outer_area"])-F(d["old"]["inner_area"]))
    el=F(d["axis-error"]["gamma_change_upper"])
    gl=glrat[0]-el,glrat[1]+el
    grrat=d["prior"]["coefficient_intervals"]["gamma_rational_axis"]
    er=F(d["prior"]["axis_error_source"]["gamma_error"])
    gr=F(grrat["lower"])-er,F(grrat["upper"])+er
    K=bx[0]+gl[0]+gr[0],bx[1]+gl[1]+gr[1]
    assert K[0]<=K[1]<-F(12,5)
    px,py=F(original["Px"]),F(original["Py"])
    m=sum(F(p[2]) for p in original["pts"])
    rho=m/(px*py)
    short=rho*(1+px/15)*(1+py/30)
    long=rho*F(3,2)*(1+py/(6*px))
    assert 6*px>30 and short<1 and long<1
    # For N>=7, left-seam window reaches no omitted column i>=N.
    assert 7*px-(X-lo)-30>0
    result={"status":"PASS","template":"R41", "python":sys.version,"platform":platform.platform(),
        "arithmetic":"Fraction; exact interval arithmetic and strict rational comparisons",
        "random_seed":None,"floating_tolerance":None,"script_sha256":sha(__file__),
        "input_sha256":oh,"reflected_input_sha256":rh,"cell_mass":str(m),"X":str(X),
        "original_true_left_axis_interval":interval((X-hi,X-lo)),
        "short_R15_tail_bound":str(short),"long_radius_tail_bound":str(long),
        "six_Px_minus_30":str(6*px-30),
        "coefficients":{name:interval(z) for name,z in [("A",A),("b_x",bx),("gamma_left_rational_axis",glrat),("gamma_left_true_axis",gl),("gamma_right_true_axis",gr),("K_double_true_axis",K)]},
        "left_axis_error":str(el),"right_axis_error":str(er),
        "strict_K_upper":"-12/5","formula":"J_double(N,lambda=1)=A/m+K/(m*N), every integer N>=7",
        "strict_threshold_one_deficit":"12/(5*m*N)",
        "scope":"Fixed canonical R41 source phase, true left/right halo axes, double-reflected period2W; fixed threshold1; no assertion about all thresholds, other phases, N<=6, R35, or arbitrary smoothing.",
        "sources":[{"role":k,"path":str(p),"sha256":sha(p),"size":p.stat().st_size} for k,p in paths.items()],
        "new_run_counts":{k:d[k]["counts"] for k in ["axis","old","new"]}}
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":"PASS","gamma_left":interval(gl),"K":interval(K),"strict_K_upper":"-12/5"}))


if __name__=="__main__":
    main()
