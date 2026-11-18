"""Construct the deepest enclosing EH chord as an extra row for all AH root tests."""
from hashlib import sha256
import json


def add_row(ns, archive):
    from certify_eh_upper_chain import front,b,g
    old_path=archive/'eh-two-line-cover-certificate.json'
    parameter_path=archive/'eh-upper-chain-parameters.json'
    old,parameters=(json.loads(p.read_bytes()) for p in (old_path,parameter_path))
    assert len(parameters)==9
    g.SIGMA=g.arb('1.5497')*g.PI/180
    shape=front(old,parameters)
    n=b.unit_outer(shape['v'],shape['t'])
    row=('EH-deep-chord',n,g.dot(n,shape['t']))
    assert all(name!=row[0] for name,_,_ in ns['rows'])
    ns['rows']=list(ns['rows'])+[row]
    return {'label':row[0],'parameters_sha256':sha256(parameter_path.read_bytes()).hexdigest(),
            'R126_sha256':sha256(old_path.read_bytes()).hexdigest(),
            'number_of_EH_facets':len(shape['active']),
            'geometry_construction':'unit_outer(deepest_enclosing_upper_endpoint,old_lower_endpoint)'}
