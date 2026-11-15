"""Identify the two raw limiter owners and certify disjoint fixed y supports."""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime,timezone
import argparse,hashlib,json,sys

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--slab',type=Path,required=True);ap.add_argument('--slope',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=json.loads(a.slab.read_text());s=json.loads(a.slope.read_text());pair=s['sweep']['limiting_pair'];owners=[]
    assert s['slab_sha256']==sha(a.slab)=='c0d262c39a01db8d50cd2fee1a00c7867aa8e08417c4fc85522103b70ce13b25'
    for endpoint in pair:
        rows=[(j,k)for j,box in enumerate(d['boxes'])for k in(0,1)if box[k]==endpoint];assert len(rows)==1;owners.append(rows[0])
    assert owners==[(25244,0),(25243,1)]
    r1,r2=[d['box_records'][j]for j,k in owners];assert r1['bbox']==[137,138,25,25]and r2['bbox']==[137,137,36,36]
    y1,y2=F(d['y_classes'][25]),F(d['y_classes'][36]);assert y1==F('35136509207263900091257231773/64000000000000000000000000000')<F(3,5)and y2==F('316077158213425564705302455727/64000000000000000000000000000')>F(24,5)
    assert 0<F(r1['radius_squared'])<F(1,100)and 0<F(r2['radius_squared'])<F(1,100)
    out={'status':'passed','scope':'limiter identity and disjoint y support only; no nearest true area wall claim','slab_sha256':sha(a.slab),'slope_sha256':sha(a.slope),'owners_zero_based_box_and_horizontal_endpoint':owners,'bbox':[r1['bbox'],r2['bbox']],'y_centers':[str(y1),str(y2)],'radius_squared':[r1['radius_squared'],r2['radius_squared']],'first_y_support_strictly_below':'7/10','second_y_support_strictly_above':'47/10','strict_support_gap_lower_bound':'4','script_sha256':sha(Path(__file__)),'python':sys.version,'command':[sys.executable,*sys.argv],'verified_utc':datetime.now(timezone.utc).isoformat()}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
