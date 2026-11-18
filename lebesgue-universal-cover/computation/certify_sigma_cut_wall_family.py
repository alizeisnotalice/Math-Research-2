"""Project genuine support points of the obstructed pentagon onto the two physical AH coordinates."""
import argparse,inspect,json,sys
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from certify_sigma_perturbed_pentagon_obstruction import run as witness


def run(archive,inp,rec,out):
    doc=json.loads(rec.read_bytes());cfg=json.loads(inp.read_bytes());assert doc['parameters']==cfg and cfg['cut_only'] and doc['source_sha256']==sha256((archive/'certify_sigma_perturbed_pentagon_obstruction.py').read_bytes()).hexdigest();src=inspect.getsource(witness);src=src[:src.index('    bounds=[]')]+'''    data=[]
    for row in doc['dihedral_witnesses']:
        if row['edge']!='AH-cut':continue
        t,s=row['turn'],row['sign'];pre=ns['rotate']((cu[0],s*cu[1]),trig(A(-t)));active=[k for k,x in enumerate(edges) if cross(x,pre)>0 and cross(pre,edges[(k+1)%10])>0];assert len(active)==1
        k,arc=ns['sector'][active[0]];point=add(V[k],mul(arc,pre));image=ns['rotate'](point,trig(A(t)));image=(image[0],s*image[1]);r=sub(image,ns['z']);a=dot(r,ns['direct'])+A(p)*ns['length'];b=dot(r,ns['normal']);base=A(d).tan()*a-b
        data.append({'turn':t,'sign':s,'point':list(map(str,image)),'unit_w_coordinate':str(-dot(r,ns['direct'])/ns['length']),'height_relative_base':str(b),'pivot_minus_w_times_length':str(a),'selected_baseline_violation_over_cos':str(base),'slope_threshold_tan_delta':str(b/a) if not a.contains(0) else None,'w_less_than_pivot':bool(a>0)});assert base>0,(a,base)
    assert len(data)==2;return data
''';env={**witness.__globals__,'doc':doc};exec(compile(src,'true_AH_support_point_for_sector','exec'),env);rows=env['run'](archive,inp,None);r={'witness':cfg,'witness_sha256':sha256(rec.read_bytes()).hexdigest(),'fixed_points':rows,'all_fixed_points_before_pivot':all(q['w_less_than_pivot'] for q in rows),'consequence':'At least one displayed height must exceed or equal |D|(p-w_i) tan(delta) for its own row. If both fixed points are before the baseline pivot, every greater or equal pivot and greater or equal angle below pi/2 is excluded.','source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()};out.write_text(json.dumps(r,indent=2)+'\n');print(rows)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('input',type=Path);p.add_argument('witness',type=Path);p.add_argument('output',type=Path);o=p.parse_args();run(o.archive.resolve(),o.input,o.witness,o.output)
