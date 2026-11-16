#!/usr/bin/env python3
"""Adversarial tests for parent-owned verify_H8_exact.py; parent files stay read-only."""
import copy,hashlib,importlib.util,json,shutil,tempfile,argparse,sys
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
def run(name,mut,P,m,OUT):
 with tempfile.TemporaryDirectory(dir=OUT.parent) as d:
  d=Path(d)
  for n in ('H8-a2b3-inputs.json','H8-a2b3-root-balls.json','H8-a2b3-certificate.json','H8-mother-certificate.json'):
   shutil.copy2(P/n,d/n)
  x={n:json.loads((d/n).read_text()) for n in ('H8-a2b3-inputs.json','H8-a2b3-root-balls.json','H8-a2b3-certificate.json','H8-mother-certificate.json')};mut(x)
  for n,v in x.items():write(d/n,v)
  # Permit deep semantic paths after input/balls mutations.
  c=x['H8-a2b3-certificate.json'];c['input_sha256']=sha(d/'H8-a2b3-inputs.json');c['balls_sha256']=sha(d/'H8-a2b3-root-balls.json');write(d/'H8-a2b3-certificate.json',c)
  try:m.verify(d/'H8-a2b3-inputs.json',d/'H8-a2b3-root-balls.json',d/'H8-a2b3-certificate.json',d/'H8-mother-certificate.json')
  except (AssertionError,KeyError,ValueError,IndexError):return True
  return False
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--package',type=Path,required=True);ap.add_argument('--checker',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();P=a.package;OUT=a.output
 sys.path.insert(0,str(P));spec=importlib.util.spec_from_file_location('parent_verify',a.checker);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 tests={
 'fixed_H':lambda x:x['H8-a2b3-inputs.json']['candidates'][0].__setitem__('H_coefficients_ascending',[1]*9),
 'mother_duplicate':lambda x:x['H8-a2b3-root-balls.json']['mothers'][1].__setitem__('block',0),
 'mother_nonpositive_imag':lambda x:x['H8-a2b3-root-balls.json']['mothers'][0]['disk']['center'].__setitem__('im','0/1'),
 'lift_radius_huge':lambda x:x['H8-a2b3-root-balls.json']['lifts'][0]['disk'].__setitem__('radius','1/1'),
 'lift_center_shift':lambda x:x['H8-a2b3-root-balls.json']['lifts'][0]['disk']['center'].__setitem__('re','0/1'),
 'lift_wrong_family':lambda x:x['H8-a2b3-root-balls.json']['lifts'][0].__setitem__('family','b_lift'),
 'lift_wrong_slot':lambda x:x['H8-a2b3-root-balls.json']['lifts'][0].__setitem__('lift_index',1),
 'negative_duplicate':lambda x:x['H8-a2b3-certificate.json']['negative_certificates'].append(copy.deepcopy(x['H8-a2b3-certificate.json']['negative_certificates'][0])),
 'negative_missing':lambda x:x['H8-a2b3-certificate.json']['negative_certificates'].pop(),
 'wrong_digits':lambda x:x['H8-a2b3-certificate.json']['negative_certificates'][0].__setitem__('digits',[2,2,2,2]),
 'fake_positive':lambda x:x['H8-a2b3-certificate.json']['negative_certificates'][0].__setitem__('hi_n','0'),
 'pseudo_standard':lambda x:x['H8-a2b3-certificate.json']['unknown'][0].__setitem__('digits',[1,0,0,0]),
 'extra_unknown':lambda x:x['H8-a2b3-certificate.json']['unknown'].append({'E':2,'code':0,'digits':[0]*4,'terminal_standard':False}),
 'galois_wrong_factor':lambda x:x['H8-mother-certificate.json']['certificates'][0].__setitem__('coefficients_ascending',[1,1]),
 }
 ans={k:run(k,v,P,m,OUT) for k,v in tests.items()};OUT.write_text(json.dumps({'parent_source_sha256':sha(a.checker),'rejected':ans,'all_rejected':all(ans.values())},indent=2)+'\n');print(OUT.read_text())
if __name__=='__main__':main()
