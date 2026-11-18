"""Verify the original P30 at the five-body placement fixed by the hard-leaf control."""
import argparse,json,sys
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_angles import verify

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('comp','pose','certificate','output'):p.add_argument(n,type=Path)
    a=p.parse_args();c=json.loads(a.certificate.read_bytes())
    center=json.loads(a.pose.read_bytes())['hard']['center']
    assert c['box']==[[v,v] for v in center] and c['five_only'] and F(c['target'])==F('.8378')
    result=verify(a.comp,a.comp/'finite-point-leaf-input.json',c,source_p30=True)
    audit=dict(fixed_curved_host_center=center,point_cloud='original P30',strict_minimum_area_lower=str(F(c['target'])),
               angle_cells=result['cases'],minimum_margin=result['minimum_normalized_strict_margin'],
               source_sha256=result['points_sha256'],certificate_sha256=sha256(a.certificate.read_bytes()).hexdigest())
    a.output.write_text(json.dumps(audit,indent=2)+'\n');print(audit)
if __name__=='__main__':main()
