"""Independent run of the two-angle stationary definition and the old and new local minimum basins."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):os.environ[key]="1"
import flint,numpy,scipy
from flint import arb,ctx

def sha(p):return sha256(p.read_bytes()).hexdigest()
def main():
    if not __debug__:raise SystemExit("Assertions required")
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("scratch",type=Path)
    p.add_argument("--contact",type=Path,default=Path(__file__).resolve().parent)
    p.add_argument("--wide",type=Path,default=Path(__file__).resolve().parent)
    a=p.parse_args();Q=a.contact.resolve();P=a.wide.resolve();ctx.prec=192
    jobs=[
        (Q/"certify_five_body_contact_stationary.py",Q/"five-body-contact-stationary.json",[]),
        (Q/"certify_five_body_contact_transverse.py",Q/"five-body-contact-transverse.json",[]),
        (Q/"certify_five_body_contact_basin.py",Q/"five-body-contact-basin.json",[]),
        (P/"certify_five_wide_contact_tube.py",P/"five-wide-contact-tube.json",[Q,"--radius",".000012","--transverse",".00007","--cube",".00012","--cells","8","--adaptive-max","128"])]
    output=[]
    with tempfile.TemporaryDirectory(dir=a.scratch)as t:
        t=Path(t)
        for script,ref,opts in jobs:
            raw=subprocess.check_output([sys.executable,"-B",script,*opts],cwd=script.parent)
            act=json.loads(raw);expected=json.loads(ref.read_bytes())
            aa={k:v for k,v in act.items()if k!="python"};bb={k:v for k,v in expected.items()if k!="python"}
            identical=aa==bb
            if script.name!="certify_five_wide_contact_tube.py":assert identical,script.name
            else:
                for k in ("center","scale","halfwidth","transverse","cube","degree","complete","initial","code_sha256"):assert aa[k]==bb[k]
                assert aa["complete"] and sum(Fraction(c,int(n)**2)for n,c in aa["accepted_cells"].items())==1
                assert arb(aa["minima"][0])>0 and Fraction(aa["minima"][1])>Fraction(".00012")
                assert arb(aa["minima"][2])>0 and arb(aa["minima"][3])>arb(".000000034")
                local_tube=t/"tube.json";local_tube.write_bytes(raw)
            output.append(dict(name=script.name,source_sha256=sha(script),reference_sha256=sha(ref),
                numerical_fields_identical=identical,independent_output_sha256=sha256(raw).hexdigest()))
            if script.name=="certify_five_body_contact_stationary.py":
                assert arb(act["area_ball"])<arb(".83729500622102525") and arb(act["area_ball"])>arb(".83729500622102524")
                assert arb(act["Hessian_determinant"])>0 and arb(act["contraction_norm_upper"])<1
                stationary={k:act[k]for k in ("stationary_angle_box","contraction_norm_upper","verified_topology","area_ball")}
        script=P/"certify_five_wide_basin_bridge.py"
        raw=subprocess.check_output([sys.executable,"-B",script,Q,"--tube",local_tube],cwd=script.parent)
        act=json.loads(raw);expected=json.loads((P/"five-wide-basin-bridge.json").read_bytes())
        assert {k:v for k,v in act.items()if k!="tube_sha256"}=={k:v for k,v in expected.items()if k!="tube_sha256"}
        output.append(dict(name=script.name,source_sha256=sha(script),reference_sha256=sha(P/"five-wide-basin-bridge.json"),
            numerical_fields_identical=True,independent_output_sha256=sha256(raw).hexdigest()))
    print(json.dumps(dict(complete=True,architecture=platform.machine(),python=platform.python_version(),flint=flint.__version__,
        numpy=numpy.__version__,scipy=scipy.__version__,stationary=stationary,independent_local_proofs=output),indent=2))
if __name__=="__main__":main()
