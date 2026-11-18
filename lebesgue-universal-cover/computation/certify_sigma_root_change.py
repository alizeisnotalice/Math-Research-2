"""An angle-dependent change from two compatible roots to none at the same AH parameters."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
from certify_sigma_paired_rigidity import rig


def verify():
    w0,w1=F('.42849'),F('.42851')
    left,right=F('1.5399'),F('1.5401')
    high0,high1=F('1.5599'),F('1.5601')
    angles=list(map(F,['53.4','113.6','114.1','54.3']))
    probe=rig(left,right,w0,w1,probe_angles=angles)
    assert min(probe)>0
    num,bounds=rig(high0,high1,w0,w1,endcap=F('0.1'),min_cell=F(1,1000000))
    return {'paired_parameter':list(map(str,(w0,w1))),
            'two_root_sigma_degrees':list(map(str,(left,right))),
            'probed_degrees_with_sixty_degree_sign_flip':list(map(str,angles)),
            'signed_margins_lower':list(map(str,probe)),
            'certified_open_root_brackets_degrees':[['53.4','53.6'],['54.1','54.3']],
            'no_extra_root_sigma_degrees':list(map(str,(high0,high1))),
            'interior_angle_cells':num,'root_slope_first_last_tenth_degree':str(bounds[-2]),
            'interior_positive_lower':str(bounds[-1])}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    args.output.write_text(json.dumps(verify(),indent=2)+'\n')
