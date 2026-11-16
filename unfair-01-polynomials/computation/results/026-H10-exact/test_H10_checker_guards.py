"""Semantic tamper tests for H10 exact checker; writes only under --scratch."""
from pathlib import Path
import argparse,copy,hashlib,json,traceback
from verify_H10_exact import verify

def write(p,x):Path(p).write_text(json.dumps(x)+"\n")
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def reject(fn):
 try:fn()
 except (AssertionError,KeyError,TypeError,ValueError) as e:return {'exception_type':type(e).__name__,'reason':str(e) or 'assertion guard','actual_last_frame':{'file':Path(traceback.extract_tb(e.__traceback__)[-1].filename).name,'line':traceback.extract_tb(e.__traceback__)[-1].lineno}}
 raise AssertionError('tampered artifact accepted')
def neg(s):return '-' + s if not s.startswith('-') else s[1:]
def main():
 p=argparse.ArgumentParser()
 for n in ('input','balls','certificate','galois','scratch'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();a.scratch.mkdir(parents=True,exist_ok=True)
 I=json.loads(a.input.read_text());B=json.loads(a.balls.read_text());C=json.loads(a.certificate.read_text());G=json.loads(a.galois.read_text())
 def run(tag,i=I,b=B,c=C,g=G):
  ip=a.scratch/(tag+'-input.json');bp=a.scratch/(tag+'-balls.json');cp=a.scratch/(tag+'-cert.json');gp=a.scratch/(tag+'-galois.json')
  write(ip,i);b=copy.deepcopy(b);c=copy.deepcopy(c);b['input_sha256']=digest(ip);write(bp,b);c['input_sha256']=digest(ip);c['balls_sha256']=digest(bp);write(cp,c);write(gp,g)
  return verify(ip,bp,cp,gp)
 R={}
 b=copy.deepcopy(B);b['mothers'][0]['disk']['radius']='-1/1';R['negative_radius']=reject(lambda:run('negative-radius',b=b))
 b=copy.deepcopy(B);b['lifts'][0]['lift_index']=99;R['lift_slot_mismatch']=reject(lambda:run('slot',b=b))
 c=copy.deepcopy(C);c['unknown'].append({'E':2,'code':0,'digits':[0]*5,'terminal_standard':False});R['extra_unknown']=reject(lambda:run('extra-unknown',c=c))
 c=copy.deepcopy(C);c['negative_certificates'].append(copy.deepcopy(c['negative_certificates'][0]));R['duplicate_negative']=reject(lambda:run('duplicate-negative',c=c))
 c=copy.deepcopy(C);c['unknown'].append(copy.deepcopy(c['unknown'][0]));R['duplicate_unknown']=reject(lambda:run('duplicate-unknown',c=c))
 b=copy.deepcopy(B);b['mothers'][0]['disk']['radius']='0/1';R['nonpositive_radius']=reject(lambda:run('zero-radius',b=b))
 b=copy.deepcopy(B);b['lifts'][0]['family']='b_lift';R['lift_family_mismatch']=reject(lambda:run('family',b=b))
 b=copy.deepcopy(B);b['mothers'][0]['disk']['center']['im']=neg(b['mothers'][0]['disk']['center']['im']);R['mother_nonpositive_imag']=reject(lambda:run('mother-imag',b=b))
 i=copy.deepcopy(I);b=copy.deepcopy(B);root=i['candidates'][0]['blocks'][0]['parent_root_positive_imag'];i['candidates'][0]['blocks'][1]['parent_root_positive_imag']=root
 for g in i['candidates'][0]['groups']:
  if g['block']==1:g['parent_root_positive_imag']=root
 b['mothers'][1]['disk']=copy.deepcopy(b['mothers'][0]['disk']);b['mothers'][1]['lhs']=b['mothers'][0]['lhs'];b['mothers'][1]['rhs']=b['mothers'][0]['rhs'];R['duplicate_mother_disk']=reject(lambda:run('mother-duplicate',i=i,b=b))
 c=copy.deepcopy(C);c['negative_certificates'].pop();R['negative_omission']=reject(lambda:run('negative-omission',c=c))
 c=copy.deepcopy(C);c['negative_certificates'][0]['digits']=[99]*5;R['wrong_digits']=reject(lambda:run('wrong-digits',c=c))
 c=copy.deepcopy(C);c['negative_certificates'][0]['hi_n']='0';R['fake_positive_bound']=reject(lambda:run('fake-positive',c=c))
 c=copy.deepcopy(C);c['unknown'][0]['terminal_standard']=False;R['pseudo_standard']=reject(lambda:run('pseudo-standard',c=c))
 g=copy.deepcopy(G);g['certificates'][0]['reductions'][0]['factors'][0]['coefficients_ascending'][0]=0;R['galois_wrong_factor']=reject(lambda:run('galois-factor',g=g))
 labels={'negative_radius':'disk radius must be positive','lift_slot_mismatch':'lift slot matches source','extra_unknown':'only canonical standard unknown permitted','duplicate_negative':'duplicate (E,code) rejected','duplicate_unknown':'duplicate unknown (E,code) rejected','nonpositive_radius':'disk radius must be strictly positive','lift_family_mismatch':'ball family/slot must equal source group','mother_nonpositive_imag':'declared parent disk must have positive imaginary part','duplicate_mother_disk':'mother-disk separation after input/binding synchronization','negative_omission':'complete 8309 negative-record cover required','wrong_digits':'certificate digits must equal independently regenerated layout','fake_positive_bound':'stored upper endpoint must be strictly negative','pseudo_standard':'sole unknown must be flagged exact terminal standard','galois_wrong_factor':'finite-field factor product/irreducibility gate'}
 for key,value in R.items():value['expected_semantic_gate']=labels[key]
 out={'status':'passed_semantic_tamper_regressions','guards':list(R),'rejection_reasons':R};write(a.scratch/'tamper-regressions.json',out);print(json.dumps(out))
if __name__=='__main__':main()
