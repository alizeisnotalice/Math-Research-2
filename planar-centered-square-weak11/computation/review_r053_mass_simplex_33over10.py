#!/usr/bin/env python3
"""Gzip wrapper for the frozen independent R53 whole-tree checker."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,gzip,hashlib,json,platform,sys,time
from review_r053_mass_simplex import run

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--certificate',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();t=time.perf_counter()
    with gzip.open(a.certificate,'rt') as f:cert=json.load(f)
    assert Q(cert['target'])==Q(33,10)
    result=run(cert)
    b=Q(result['uniform_area_upper_exact']);assert b<Q(33,10)
    assert Q(33,10)-b==Q(cert['strict_margin_below_target'])
    # Explicit local-to-expanded-bbox class: D_bbox^2 <= (35/33)M gives <3.5M.
    assert b*Q(35,33)<Q(7,2)
    # A simple margin-paying nine-tight-block class, delta <= sqrt(M)/20000.
    delta=Q(1,20000);error=4*511*(delta+delta*delta)
    assert b+error<Q(7,2)
    result.update(target='33/10',strict_margin_below_target=str(Q(33,10)-b),expanded_bbox_area_ratio='35/33',tight_nine_block_delta_ratio='1/20000',tight_nine_block_contraction_error_ratio=str(error),python=sys.version,platform=platform.platform(),elapsed_seconds=time.perf_counter()-t,source_sha256={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [a.certificate,Path(__file__).resolve(),Path(__file__).with_name('review_r053_mass_simplex.py')]},scope='Whole closed finite mass simplex, exact stronger target33/10; arbitrary-bbox general C2 remains open. Additional rational checks support stated conditional geometric extensions only.')
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
