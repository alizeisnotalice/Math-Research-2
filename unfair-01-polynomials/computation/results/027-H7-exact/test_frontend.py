"""Deep semantic tests without source hash gates; all mutations are in memory."""
import argparse,copy,json,traceback
from pathlib import Path
from verify_h7_frontend import verify,layouts,ALL

def main():
 p=argparse.ArgumentParser();p.add_argument('--balls',type=Path,required=True);p.add_argument('--ledger',type=Path,required=True);p.add_argument('--record',type=Path,required=True);a=p.parse_args();B=json.loads(a.balls.read_text());L=json.loads(a.ledger.read_text());rows=[]
 def check(name,mut,gate):
  b=copy.deepcopy(B);l=copy.deepcopy(L);mut(b,l)
  try:verify(b,l)
  except (ValueError,AssertionError,KeyError) as e:
   assert str(e)==gate,(name,str(e),gate);rows.append({'name':name,'semantic_guard':str(e),'frame_line':traceback.extract_tb(e.__traceback__)[-1].lineno})
  else:raise AssertionError('accepted '+name)
 check('negative_radius',lambda b,l:b['mothers'][0]['disk'].__setitem__('radius','-1/1'),'positive radius')
 check('mother_kind_swap',lambda b,l:b['mothers'][3].__setitem__('kind','complex_mother'),'mother kind binding')
 check('duplicate_mother_block',lambda b,l:b['mothers'][1].__setitem__('block',0),'unique mother block')
 check('family_mismatch',lambda b,l:b['lifts'][0].__setitem__('family','b'),'group mapping binding')
 check('real_linear_as_quadratic',lambda b,l:b['lifts'][21].__setitem__('kind','quadratic'),'group mapping binding')
 check('duplicate_lift_slot',lambda b,l:b['lifts'][1].__setitem__('lift_index',b['lifts'][0]['lift_index']),'unique lift slot')
 def duplicate_disk(b,l):
  for key in ('disk','lhs','rhs'):b['lifts'][1][key]=copy.deepcopy(b['lifts'][0][key])
 check('duplicate_disk_synchronized_rouche',duplicate_disk,'disk separation')
 check('repeated_E_ledger',lambda b,l:l['rows'].__setitem__(1,copy.deepcopy(l['rows'][0])),'unique ledger rows')
 check('wrong_real_choice_same_count',lambda b,l:l['rows'][0].__setitem__('real_choices',[[0,1,0]]),'ledger local choices')
 check('terminal_complement_convention',lambda b,l:l['rows'][-1].__setitem__('terminal_fixed_first_complex_all_a',False),'terminal convention')
 check('wrong_lift_rouche_bound',lambda b,l:b['lifts'][0].__setitem__('lhs','0/1'),'lift Rouche')
 d,m,_=verify(B,L);assert len(m)==12158 and len(set(mask for E,mask in m))==12158
 out={'status':'passed','tests':rows,'negative_tests':len(rows),'valid_frontend':d,'coverage_note':'Independent 25-bit masks include real slots; complement quotient exact.'};a.record.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'passed','negative_tests':len(rows),'layouts':len(m)}))
if __name__=='__main__':main()
