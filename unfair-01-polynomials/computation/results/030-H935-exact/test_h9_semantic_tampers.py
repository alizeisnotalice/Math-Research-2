"""Direct production verify(data) regressions for H9 root-cover semantics."""
import copy,json
from pathlib import Path
from verify_h9_roots import verify

SOURCE=Path(__file__).with_name('H9-35-root-balls.json')
def rejected(label,data,triggered):
 try:verify(data)
 except (ValueError,KeyError,TypeError) as error:
  triggered[label]=str(error);return
 raise AssertionError(f'accepted: {label}')
def main():
 base=json.loads(SOURCE.read_text());got={}
 x=copy.deepcopy(base)
 for key in ('disk','lhs','rhs'):x['mothers'][3][key],x['mothers'][4][key]=x['mothers'][4][key],x['mothers'][3][key]
 rejected('mother_conjugate_binding',x,got)
 x=copy.deepcopy(base);x['groups'][1]['index']=0;rejected('duplicate_group_index',x,got)
 x=copy.deepcopy(base);x['groups'][0]['source']='unknown';rejected('unknown_group_source',x,got)
 x=copy.deepcopy(base);next(g for g in x['groups'] if g['source']=='G5_lift' and g['exponent']==3)['exponent']=5;rejected('wrong_lift_exponent',x,got)
 x=copy.deepcopy(base);x['groups'][0]['disk']['radius']='1/10';rejected('expanded_cyclotomic_radius',x,got)
 out={'status':'H9_production_verify_semantic_tampers_passed','triggers':got};Path(__file__).with_name('H9-semantic-tampers.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
