"""Semantic mutations through the actual modular verification entry point."""
import copy,json
from pathlib import Path
from verify_h9_modular import verify
HERE=Path(__file__).resolve().parent
def main():
    data=json.loads((HERE/'G5-cycle-results.json').read_text());checks=[]
    for label in ['wrong_product','reducible_factor','wrong_multiplicity','wrong_cycle_degree']:
        x=copy.deepcopy(data)
        if label=='wrong_product':x['G5_substitution_factors']['307']['G_x3'][0]['ascending'][0]+=1
        if label=='reducible_factor':x['irreducible5']['factors'][0]['ascending']=[0,0,0,0,0,1]
        if label=='wrong_multiplicity':x['G5_substitution_factors']['1201']['Phi25'][0]['multiplicity']=2
        if label=='wrong_cycle_degree':x['G5_substitution_factors']['307']['G_x3'][0]['degree']=3
        try:verify(x)
        except ValueError as err:checks.append({'test':label,'rejected':str(err)})
        else:raise RuntimeError(label+' accepted')
    (HERE/'H9-modular-negative-tests.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(checks)
if __name__=='__main__':main()
