"""Semantic coefficient-certificate mutations; source hashes are not a gate."""
from pathlib import Path
import argparse,copy,json,traceback
from verify_h7_certificate import check

def main():
 p=argparse.ArgumentParser()
 for key in ('balls','ledger','certificate','galois','scratch'):p.add_argument('--'+key,type=Path,required=True)
 a=p.parse_args();a.scratch.mkdir(exist_ok=True,parents=True);C=json.loads(a.certificate.read_text());rows=[]
 def run(name,mut,expected):
  c=copy.deepcopy(C);mut(c);path=a.scratch/(name+'.json');path.write_text(json.dumps(c));out=a.scratch/(name+'-receipt.json')
  try:check(a.balls,a.ledger,path,a.galois,out)
  except (ValueError,AssertionError) as e:
   assert str(e)==expected,(name,str(e));rows.append({'name':name,'guard':str(e),'line':traceback.extract_tb(e.__traceback__)[-1].lineno})
  else:raise AssertionError('accepted '+name)
 run('duplicate_mask',lambda c:c['negative_certificates'].__setitem__(1,copy.deepcopy(c['negative_certificates'][0])),'unique complete mask records')
 run('missing_negative',lambda c:c['negative_certificates'].pop(),'record counts')
 run('invalid_index',lambda c:c['negative_certificates'][0].__setitem__('index',999),'coefficient index')
 run('invalid_side',lambda c:c['negative_certificates'][0].__setitem__('side','C'),'side')
 run('false_negative_at_constant',lambda c:c['negative_certificates'][0].__setitem__('index',0),'independent exact negative coefficient')
 run('nonnegative_proposed_upper',lambda c:c['negative_certificates'][0].__setitem__('hi_n','0'),'proposed negative interval')
 run('partial_marked_certificate',lambda c:c.__setitem__('status','running'),'full scan status')
 def nonstandard(c):
  r=next(r for r in c['negative_certificates'] if r['E']==12);r['mask'],c['unknown'][0]['mask']=c['unknown'][0]['mask'],r['mask']
 run('fake_standard_complete_cover',nonstandard,'exact standard root set')
 (a.scratch/'negative-tests.json').write_text(json.dumps({'status':'passed','tests':rows,'count':len(rows)},indent=2)+'\n');print(json.dumps({'status':'passed','count':len(rows)}))
if __name__=='__main__':main()
