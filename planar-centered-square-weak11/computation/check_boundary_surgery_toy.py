#!/usr/bin/env python3
"""Compare full surgery-collar output with the isolated-lattice closed form."""
from pathlib import Path
from fractions import Fraction as F
import argparse,math,json,os,hashlib,sys

def main():
    ap=argparse.ArgumentParser();ap.add_argument('area');ap.add_argument('--output',required=True);a=ap.parse_args()
    assert os.path.ismount('/Volumes/data');p=Path(a.area);d=json.loads(p.read_text())
    assert d['side']=='left' and F(d['t'])==F(1,2)
    den=1<<100;k=math.isqrt(2*den*den);lo,hi=F(k,den),F(k+1,den)
    expected=F(15,4)-hi/2,F(15,4)-lo/2
    actual=F(d['true_axis_area_lower']),F(d['true_axis_area_upper'])
    assert actual[0]<=expected[0]<=expected[1]<=actual[1]
    K=2*(actual[0]-3),2*(actual[1]-3);assert K[0]>F(2,25)
    result={'status':'PASS','python':sys.version,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':d['input_sha256'],'area_output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'formula':'half-collar area=3+(3-2sqrt(2))/4; balanced K=(3-2sqrt(2))/2>2/25','expected_exact_interval':list(map(str,expected)),'computed_true_area_interval':list(map(str,actual)),'K_interval':list(map(str,K)),'numerical_display':list(map(float,K)),'scope':'Analytically isolated Px=Py=10 unit-atom toy at t=1/2. TailR2 excludes all but singleton and seam adjacent pair; no real-template conclusion.','random_seed':None,'floating_tolerance':None}
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'PASS','K_display':result['numerical_display']}))
if __name__=='__main__':main()
