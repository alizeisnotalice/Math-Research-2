"""Construct or replay an all-isometry Q11 obstruction at the six-body escape pose."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from certify_multi_host_points import SupportHull, points, q, verify_rows

TARGET = Fraction(".0022")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('pose','cloud','receipt'):
        p.add_argument(name,type=Path)
    p.add_argument('--generate',action='store_true')
    p.add_argument('--padding',type=Fraction,default=TARGET)
    a=p.parse_args()
    pose_raw=a.pose.read_bytes()
    points_raw=a.cloud.read_bytes()
    pose=json.loads(pose_raw)[:11]
    data=json.loads(points_raw)
    cloud,max_d2=points(data)
    host=SupportHull(list(map(q,pose)))
    inputs={
        'host_sha256': sha256(pose_raw).hexdigest(),
        'Q11_sha256': sha256(points_raw).hexdigest(),
        'threshold':str(TARGET),
    }
    if a.generate:
        import numpy as np
        from generate_multi_host_points import Sampler, suggest
        sampled=Sampler([np.asarray(pose,dtype=float)],720,720).host_table[0]
        raw=np.array(data['points'],dtype=float)/data['point_denominator']
        cert={str(sign):suggest(host,sampled,raw,cloud,sign,TARGET) for sign in (1,-1)}
        assert all(cert.values())
        a.receipt.write_text(json.dumps({**inputs,'hands':cert},separators=(',',':'))+'\n')
    cert=json.loads(a.receipt.read_bytes())
    assert all(cert[k]==v for k,v in inputs.items())
    assert a.padding >= TARGET
    results={str(s):verify_rows(host,cloud,s,cert['hands'][str(s)],a.padding) for s in (1,-1)}
    print(json.dumps({**inputs,'Q11_max_pair_square_numerator':max_d2,'hands':results,
                      'certificate_sha256':sha256(a.receipt.read_bytes()).hexdigest()},indent=2))

if __name__=='__main__':
    main()
