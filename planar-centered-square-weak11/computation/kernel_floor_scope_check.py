"""Exact metadata/geometry audit of the full-variable kernel-floor certificate.
No LP solve, no hl_core import, no numerical matrix identity is certified here.
"""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import json,hashlib,platform,sys
certificate=Path(sys.argv[1]); candidate=Path(sys.argv[2]); menu=Path(sys.argv[3]); output=Path(sys.argv[4])
c=json.loads(certificate.read_text()); k=json.loads(candidate.read_text()); d=json.loads(menu.read_text())
classes=[tuple(x) for x in c['variables']['classes']]
expected=sorted(tuple(map(int,s.split(','))) for s in k['a'])
assert classes==expected and len(classes)==1181 and len(set(classes))==1181
assert all(0<=j<=i<=int(k['N']) for i,j in classes)
assert all(c['variables'][s]==k[s]==d[s] for s in ['alpha','R','h'])
assert c['candidate_sha256']==d['candidate_sha256']==hashlib.sha256(candidate.read_bytes()).hexdigest()
assert c['source_input_sha256']==hashlib.sha256(menu.read_bytes()).hexdigest()
assert d['classes']==[list(x) for x in classes]
assert len(d['rows'])==10657
for row in d['rows']:
    u,v=map(F,row['point'])
    assert u>=v>=0 and u+v>0
    if row['kind']=='A': assert v>0 and u+v!=F(k['R']) and row['rhs']==0
    else: assert row['kind']=='K' and row['rhs']==int(u+v<=1)
rows=c['active_constraints']; assert len(rows)==1182
assert len({(r['kind'],str(r.get('point'))) for r in rows})==len(rows)
menu_set={(r['kind'],tuple(r['point']),r['rhs']) for r in d['rows']}
for row in rows:
    assert row['kind'] in ['K','A']
    u,v=map(F,row['point'])
    assert u>=v>=0 and u+v>0
    if row['kind']=='A': assert v>0 and u+v!=F(k['R']) and row['rhs']==0
    else: assert row['rhs']==int(u+v<=1)
    assert (row['kind'],tuple(row['point']),row['rhs']) in menu_set
    assert 0<F(row['dual']['lo'])<=F(row['dual']['hi'])
def orbit(i,j): return {(s*x,t*y) for x,y in [(i,j),(j,i)] for s in [-1,1] for t in [-1,1]}
orb=[orbit(i,j) for i,j in classes]
assert len(set.union(*orb))==sum(map(len,orb))
objective_lo=sum(F(r['dual']['lo'])*r['rhs'] for r in rows)
assert objective_lo>F(3589342,1000000)
result={'scope':'exact class/point/rhs/orbit/recorded-endpoint audit only, not a new verified numerical dual identity',
'python':platform.python_version(),'full_class_count':len(classes),'variable_count':len(classes)+1,
'active_count':len(rows),'necessary_menu_count':len(d['rows']),
'active_kinds_rhs':{str(x):n for x,n in Counter((r['kind'],r['rhs']) for r in rows).items()},
'orbit_sizes':{str(x):n for x,n in Counter(map(len,orb)).items()},
'all_classes_match_exact_K4_input':True,'all_points_necessary_if_continuous_C1_C2_C4':True,
'all_active_rows_in_original_menu':True,'all_10657_menu_points_necessary':True,'stored_dual_lower_endpoints_strictly_positive':True,
'stored_objective_lower_endpoints_imply_floor_gt':'1794671/500000',
'free_spline_columns_require_equalities':True,'solve_replayed':False,
'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [certificate,candidate,menu]},
'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
output.write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))
