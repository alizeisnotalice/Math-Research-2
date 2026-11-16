import json,hashlib,argparse,sys
from fractions import Fraction
from pathlib import Path
sys.set_int_max_str_digits(100000)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 a=argparse.ArgumentParser();a.add_argument('--input',type=Path,required=True);a.add_argument('--output',type=Path,required=True);a.add_argument('--roots-input',type=Path,required=True);a.add_argument('--balls',type=Path,required=True);x=a.parse_args();d=json.loads(x.input.read_text());D=1<<40;o=[]
 for r in d['negative_certificates']:
  lo,hi=Fraction(r['lower_bound']),Fraction(r['upper_bound']);L=lo.numerator*D//lo.denominator;U=-((-hi.numerator*D)//hi.denominator)
  if U>=0:raise SystemExit('grid loses sign')
  o.append({k:r[k] for k in ('E','code','side','coefficient_index','digits')}|{'lo_n':str(L),'hi_n':str(U)})
 z={'format':'H8-a2b3-negative-Q-compact-v2','source_sha256':sha(x.input),'denominator_power':40,'negative_certificates':o,'unknown':d['unknown'],'input_sha256':sha(x.roots_input),'balls_sha256':sha(x.balls),'H':[1,1,1,0,0,0,0,1,1],'a':2,'b':3,'codes_total':1482};x.output.write_text(json.dumps(z,indent=2)+'\n');print(len(o),sha(x.output))
if __name__=='__main__':main()
