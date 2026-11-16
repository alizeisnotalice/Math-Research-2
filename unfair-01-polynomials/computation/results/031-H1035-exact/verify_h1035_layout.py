"""Independent exhaustive 2^10 fiber-choice coverage and encoding cross-check."""
import json,itertools,copy
from pathlib import Path
HERE=Path(__file__).resolve().parent
def expected():
    raw=[];canonical=set()
    for v in itertools.product((0,1),repeat=10):
        eps=list(v[:5]);eta=list(v[5:]);sums=[x+y for x,y in zip(eps,eta)]
        if len(set(sums))!=1:continue
        e=sums[0];raw.append((tuple(eps),tuple(eta)))
        if e==2 or (e==1 and eps[0]==0):
            eps=[1-x for x in eps];eta=[1-x for x in eta];e=2-e
        canonical.add((e,sum(x<<i for i,x in enumerate(eta))))
    assert len(raw)==34 and len(canonical)==17
    return canonical
def verify(layouts,coefficients,search):
    rows=layouts['layouts'];keys=[(x['e'],x['p5_mask']) for x in rows]
    assert len(keys)==17 and len(set(keys))==17 and set(keys)==expected()
    transformed={("e0",0) if e==0 else ("e1",((~mask)&30)>>1) for e,mask in keys}
    actual={(x['e'],x['pair_3_bits']) for x in coefficients['records']}
    assert transformed==actual and len(coefficients['records'])==17
    negative_search={x['mask']^15 for x in search['rows'] if x['status']=='negative'}
    negative_verify={x['pair_3_bits'] for x in coefficients['records'] if x['status']=='rejected'}
    assert negative_search==negative_verify==set(range(15))
    return {'all_1024_raw_bit_choices_checked':True,'balanced_raw_choices':34,'exchange_classes':17,
            'layout_set_equals_independent_enumeration':True,'search_and_fraction_encodings_crosschecked':True,
            'standard_search_mask':0,'standard_fraction_pair3_bits':15}
def main():
    l=json.loads((HERE/'H10-17-layouts.json').read_text());c=json.loads((HERE/'parent-coefficient-verification.json').read_text());s=json.loads((HERE/'H1035-search-witnesses.json').read_text())
    out=verify(l,c,s);tests=[]
    for mode in ['missing','duplicate','wrong_encoding']:
        bad=copy.deepcopy(l)
        if mode=='missing':bad['layouts'].pop()
        if mode=='duplicate':bad['layouts'][-1]=copy.deepcopy(bad['layouts'][-2])
        if mode=='wrong_encoding':bad['layouts'][1]['p5_mask']=1
        try:verify(bad,c,s)
        except AssertionError:tests.append(mode)
        else:raise RuntimeError('mutation accepted')
    out['production_mutations_rejected']=tests
    (HERE/'parent-layout-verification.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
if __name__=='__main__':main()
