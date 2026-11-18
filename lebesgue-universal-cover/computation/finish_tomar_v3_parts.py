"""Reverify remaining original leaf slots at the higher target, with exact boxes and fresh Arb."""
import argparse
from fractions import Fraction as F
import hashlib,json,os
from pathlib import Path
import subprocess,sys,time
import numpy as np
from flint import arb

p=argparse.ArgumentParser(description=__doc__)
for key in ("checkout","certdir","stage","outdir"):p.add_argument(key,type=Path)
p.add_argument("--parts",default=None);p.add_argument("--portable",action="store_true");p.add_argument("--depth",default=12,type=int)
a=p.parse_args()
sys.path.insert(0,str(Path(__file__).resolve().parent))
from tomar_leaf_box import leaf_box
from verify_tomar_fan_table import PINNED_MANIFEST, verified_table
up=a.checkout.resolve();sys.path.insert(0,str(up/"checker"))
import lbcheck
def sha(path):
    with path.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
pinned="10af98606c98ca9ffd0d9f402cae1e2efe74b792"
assert subprocess.check_output(["git","-C",str(up),"rev-parse","HEAD"],text=True).strip()==pinned
assert subprocess.check_output(["git","-C",str(up),"show",pinned+":checker/lbcheck.py"])==(up/"checker/lbcheck.py").read_bytes()
man=a.certdir/"SHA256SUMS";assert sha(man)==PINNED_MANIFEST
expected={name.removeprefix("./"):value for value,name in (l.split() for l in man.read_text().splitlines())}
a.outdir.mkdir(parents=True,exist_ok=True)
parts=sorted(int(x.stem[4:]) for x in a.stage.glob("part*.json")) if a.parts is None else sorted(set(map(int,a.parts.split(","))))
rows=[];start=time.monotonic()
for part in parts:
    stem=f"part{part:04}"
    input_receipt=a.stage/(stem+".json");prior=json.loads(input_receipt.read_bytes());prev=a.stage/(stem+".npz")
    assert prior["part"]==part and prior["remaining_sha256"]==sha(prev) and prior["prior_unresolved_sha256"] is None
    assert prior["tested"]==prior["all_original"]==prior["closed"]+prior["remaining"]
    assert prior["initial_mask"]=="all original nodes with code 253, 254 or 255"
    here=Path(__file__).parent
    driver="run_tomar_portable_parts.py" if a.portable else "run_tomar_integer_parts.py"
    assert prior["library_sha256"]=="ae3814e8832165e2ead77af3def4f03948dddf78603f06f51ac8b23cbe7be525"
    assert prior["table_sha256"]=="0c5bcf5dc69d8e8398dfbdb209f1563d5623523ffd5bca21ae6da58b6d92d88c"
    assert prior["target"]=="0.83661"
    for name,field in [(driver,"driver_sha256"),("tomar_integer_fan.cpp","cpp_source_sha256"),
                       ("verify_tomar_fan_table.py","verifier_sha256"),("tomar_fan_table.py","generator_sha256")]:
        assert sha(here/name)==prior[field]
    assert prior["input_manifest_sha256"]==PINNED_MANIFEST
    target=prior["target"];pos=np.load(prev,allow_pickle=False)["pos"]
    assert len(pos)==prior["remaining"] and np.array_equal(pos,np.unique(pos))
    root=a.certdir/(stem+".json");codefile=a.certdir/(stem+".npy")
    assert sha(root)==expected[root.name]==prior["header_sha256"] and sha(codefile)==expected[codefile.name]==prior["codes_sha256"]
    meta=json.loads(root.read_bytes());codes=np.load(codefile,allow_pickle=False)
    prefix=np.r_[1,1+np.cumsum(np.where(codes<11,1,-1))]
    assert (prefix[:-1]>0).all() and prefix[-1]==0 and np.count_nonzero(codes>=253)==prior["tested"]
    out=a.outdir/(stem+".json")
    basis={"target":target,"part":part,"source_sha256":sha(Path(__file__)),
           "decoder_sha256":sha(Path(__file__).with_name("tomar_leaf_box.py")),
           "input_receipt_sha256":sha(input_receipt),"prior_unresolved_sha256":sha(prev),"depth":a.depth,
           "upstream_checker_sha256":sha(up/"checker/lbcheck.py")}
    if out.exists():
        record=json.loads(out.read_bytes())
        assert all(record[k]==v for k,v in basis.items());rows.append(record);continue
    side=a.certdir/(stem+"_multi.npz")
    if side.is_file():
        assert sha(side)==expected[side.name]
        with np.load(side,allow_pickle=False) as z:
            ids=z["pos"];offset=np.r_[0,np.cumsum(z["nstruct"])];sh=z["shifts"];lam=z["lam"]
            assert np.array_equal(ids,np.flatnonzero(codes==253)) and offset[-1]==len(sh)==len(lam)
    else:
        assert not (codes==253).any()
    final=[]
    for idx0 in pos:
        idx=int(idx0);tag=int(codes[idx]);c,h=leaf_box(prefix,codes,part,idx,meta)
        checker=lbcheck.Checker(meta["ns"],float(target),meta["N"],meta["eps"])
        assert repr(float(target))==target
        calls=0;deepest=0
        if tag==253:
            match=np.searchsorted(ids,idx)
            assert ids[match]==idx
            sl=slice(int(offset[match]),int(offset[match+1]))
        def prove(cc,hh,dep):
            global calls,deepest
            calls+=1;deepest=max(deepest,dep)
            if tag==253:
                ok,msg=checker.check_multi_leaf(cc,hh,sh[sl],lam[sl])
                if not ok and dep:
                    adjusted=(c-cc+sh[sl]*h)/hh
                    ok,msg=checker.check_multi_leaf(cc,hh,adjusted,lam[sl])
                if not ok:ok,msg=checker.check_bound_leaf(cc,hh)
            elif tag==254:ok,msg=checker.check_bound_leaf(cc,hh)
            else:
                assert tag==255
                # The exact dyadic endpoints below match their binary64 sums.
                l=cc-hh;u=cc+hh
                assert all(F(x)-F(y)==F(z) for x,y,z in zip(cc,hh,l))
                assert all(F(x)+F(y)==F(z) for x,y,z in zip(cc,hh,u))
                ok=checker.check_far_leaf(l,u,meta["T"])
            if ok:return True
            if dep==a.depth:return False
            k=int(np.argmax(hh*np.array([1,1,.5,1,1,.5,1,1,.5,1,1])));half=hh.copy();half[k]/=2
            assert F(half[k])*2==F(hh[k])
            for sign in (-1,1):
                child=cc.copy();child[k]+=sign*half[k]
                assert F(child[k])==F(cc[k])+sign*F(half[k])
                if not prove(child,half,dep+1):return False
            return True
        ok=prove(c,h,0)
        final.append([idx,tag,calls,deepest,ok])
    record={**basis,"all_original":prior["tested"],"closed_by_stage":prior["closed"],
            "arb_attempts":len(final),"arb_closed":sum(row[4] for row in final),
            "unresolved":[row[0] for row in final if not row[4]],"leaves":final}
    tmp=out.with_suffix(".partial");tmp.write_text(json.dumps(record,indent=2)+"\n");tmp.replace(out)
    rows.append(record)
print(json.dumps({"parts":parts,"total_original":sum(r["all_original"] for r in rows),
    "closed_by_integer":sum(r["closed_by_stage"] for r in rows),"arb_checked":sum(r["arb_attempts"] for r in rows),
    "arb_passed":sum(r["arb_closed"] for r in rows),"remaining":sum(len(r["unresolved"]) for r in rows),
    "new_wall_seconds":round(time.monotonic()-start,3),"rows":[[r["part"],r["all_original"],r["closed_by_stage"],r["arb_closed"],len(r["unresolved"])] for r in rows]},indent=2))
