import pathlib,json,hashlib,importlib.util,time,sys,os
from fractions import Fraction as F
assert os.path.ismount('/Volumes/data')
b=pathlib.Path(__file__).resolve().parent
src=pathlib.Path('/Volumes/data/CodexProjects/T2/.local/compute/frontier063/exact_source_frontier.py')
inp=pathlib.Path('/Volumes/data/CodexProjects/T2/.local/workspaces/pcw-astra-covering/planar-centered-square-weak11/computation/inputs/r064-fixed-position-positive21.json')
assert hashlib.sha256(src.read_bytes()).hexdigest()=='5e90e5f06926b6716f56c2194b0cedeadf15d761b85173e3829d707ccc8cdf36'
assert hashlib.sha256(inp.read_bytes()).hexdigest()=='46acda723e8590054912a18f71c67efc95b5a0bcfd449a9c60ce615413ac41dd'
spec=importlib.util.spec_from_file_location('frozen_r63',src);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
d=json.loads(inp.read_text());scale=F(16);pts=[[F(x)*4*scale for x in p] for p in d['points']];m=[F(x)*16*scale**2 for x in d['initial_masses']];assert all(x.denominator==1 for p in pts for x in p)and all(x.denominator==1 for x in m)
pts=[[int(x) for x in p]for p in pts];m=list(map(int,m));t=time.monotonic();bits=40
i,o,c,g=mod.rectangles(pts,m,bits);z=mod.deletions(pts,m,bits,i,o,0,True);comp=mod.eligible_components(pts,m,o)
assert z['all_maximum_strictly_eligible_masks']==[d['S_mask']]
assert len(comp)==1 and not comp[0]['terminal_35_33']
r={'status':'PASS_INPUT_SPECIFIC_SEED_DELETION','python':sys.version,'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'geometry_source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'driver_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'spatial_scale':str(scale),'mass_scale':str(scale**2),'scaling':'x->16x and m->256m preserve threshold1 and all area/mass ratios','call_scope':'Reusable rectangle/deletion functions only; no allcut arrays or run_case N<=20 CLI; 21 source geometry uses arbitrary Python integer masks.','bits':bits,'closure_count':c,'eligible_rectangles':len(o),'Gamma_original_area':str(F(g,2**(2*bits))/scale**2),'components':comp,'deletions':z,'wall_seconds':time.monotonic()-t,'scope':'One fixed 21-positive-source seed; no optimization or universal deletion bound.'}
(b/'seed-deletion.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'status':r['status'],'closures':c,'ratios':[z['max_mass_min_lower'],z['max_mass_min_upper']],'ratio_display':float(F(z['max_mass_min_upper'])),'seconds':r['wall_seconds']},indent=2))
