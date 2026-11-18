"""Restrict transformed root tests to reference cells containing a physical first hit in at least one orientation."""
import argparse,inspect,json,sys
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from flint import arb
from certify_sigma_ah_full_chart import geometry
from certify_ah_same_fan_cancellation import support_bound
from certify_ah_local_physical_support import validate as ordinary
from certify_ah_deep_fan_chord import validate as complete
from certify_ah_affine_root_class import validate as conditional


def reach(ns,j):
    A=ns['A'];x,y=ns['knots'][j:j+2];step=(y-x)/75;used=[]
    for i in range(75):
        ends=(x+i*step,x+(i+1)*step);ell=[ns['chart'](ns['add'](ns['z'],ns['mul'](A(w),ns['D'])))[0] for w in ends]
        ell=ell[0].union(ell[1]);pair=ns['coord_hull'](ns['psi'],ell)
        if any(phys.overlaps(ell) or phys.overlaps(pair) for phys,_ in ns['physical_deficits']):used.append(i)
    starts=[i for i in used if i-1 not in used];ends=[i+1 for i in used if i+1 not in used]
    return list(zip(starts,ends,strict=True))


def run(archive,inp,out):
    cfg=json.loads(inp.read_bytes());p,d=map(F,cfg['line']);far,negative,eps=map(F,cfg['physical_parameters'])
    ns,src,_=geometry(archive,[(p,d)],sigma=cfg['sigma_degrees'],far=far,negative=negative,eps=eps)
    env={**ordinary.__globals__,'support_bound':support_bound};exec(compile(inspect.getsource(ordinary),'reachable_sigma_conditioning','exec'),env);direct=env['validate']
    rows=[]
    for j,l,r,t,s,sg in cfg['bins']:
        ranges=reach(ns,j);sub=[]
        if sg==2:
            assert ranges==[(0,75)];sub=[complete(ns,j,F(l),F(r),t,s,sg,angular=12,depth=4)]
        else:
            for a,b in ranges:sub.append(direct(ns,j,F(l),F(r),t,s,sg,begin=a,end=b,angular=12,depth=4))
        passed=all(arb(q['root_margin_lower'])>0 and arb(q['placement_margin_lower'])>0 for q in sub)
        row={'reference_bin':j,'reachable_first_hit_subcell_ranges':ranges,'tested':sub,'passed':passed};rows.append(row)
        print(j,ranges,[(q['root_margin_lower'],q['placement_margin_lower']) for q in sub],flush=True)
    obj={'parameters':cfg,'rows':rows,'max_physical_deficit':str(max(ns['deficit_by_bin'])),'all_displayed_bins_passed':all(q['passed'] for q in rows),'physical_geometry_sha256':sha256(src.encode()).hexdigest(),'input_sha256':sha256(inp.read_bytes()).hexdigest(),'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'flint':ns['flint_version'],'bits':128}
    out.write_text(json.dumps(obj,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();run(opt.archive.resolve(),opt.input,opt.output)
