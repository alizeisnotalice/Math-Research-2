"""Binary ranks of restricted Boolean monomials, for the nine R32 classes."""
import json
from pathlib import Path
p=Path(__file__).resolve().parent
def rank(rows):
    basis={}
    for row in rows:
        while row:
            pivot=row.bit_length()
            if pivot not in basis:
                basis[pivot]=row
                break
            row^=basis[pivot]
    return len(basis)
data=json.loads((p.parent/'round032/rm48_weight32_representatives.json').read_text())
result={}
for rep in data['representatives']:
    m=rep['affine_span_dimension']
    terms=[sum(1<<i for i in mon) for mon in rep['monomials_zero_based']]
    S=[x for x in range(1<<m) if (rep['constant_term']+sum(x&t==t for t in terms))%2]
    assert len(S)==32
    ranks=[rank(sum(1<<i for i,x in enumerate(S) if x&t==t) for t in range(1<<m) if t.bit_count()<=d) for d in range(5)]
    result[rep['id']]=ranks
assert result['span7_code11']==[1,8,24,31,32]
(p/'evaluation_ranks.json').write_text(json.dumps(result,indent=2)+'\n')
print(result)
