import copy,json
from pathlib import Path
from verify_h10_modular import verify
HERE=Path(__file__).resolve().parent
def main():
    m=json.loads((HERE/'H10-mother-certificate.json').read_text());c=json.loads((HERE/'H10-double-cycles.json').read_text());tests=[]
    for name in ['bad_M','wrong_factor_product','reducible','wrong_multiplicity']:
        mm,cc=copy.deepcopy(m),copy.deepcopy(c)
        if name=='bad_M':cc['certificates'][0]['residual_degree_lcm_M']=8
        if name=='wrong_factor_product':cc['certificates'][0]['families']['3'][0]['coefficients_ascending_mod_ell'][0]+=1
        if name=='reducible':mm['certificates'][0]['reductions'][0]['factors'][0]['coefficients_ascending']=[0]*10+[1]
        if name=='wrong_multiplicity':cc['certificates'][1]['families']['5'][0]['multiplicity']=2
        try:verify(mm,cc)
        except ValueError as ex:tests.append({'test':name,'rejected':str(ex)})
        else:raise RuntimeError('mutation accepted')
    (HERE/'H10-modular-negative-tests.json').write_text(json.dumps(tests,indent=2)+'\n');print(tests)
if __name__=='__main__':main()
