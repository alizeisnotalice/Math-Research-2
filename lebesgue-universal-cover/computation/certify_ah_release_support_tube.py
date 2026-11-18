"""Quantify a Hausdorff and arbitrary curvature-density neighborhood for the tapered AH host."""
import argparse
from hashlib import sha256
import inspect,json,sys
from pathlib import Path
from flint import arb,ctx
import certify_ah_tapered_release as base


def run(archive):
    src=inspect.getsource(base.run)
    marker='h>((1-t)*val+t*tri)'
    assert src.count(marker)==1
    eps='0.00000075'
    new=src.replace(marker,f"h>((1-t)*val+t*tri+A('{eps}'))")
    ns=dict(vars(base))
    exec(compile(new,'uniform_true_support_radius','exec'),ns)
    j=ns['run'](archive)
    ctx.prec=160
    r3=arb(3).sqrt()
    simplex=arb(eps)*arb('.06')
    disk=arb('.0046845')
    assert arb(j['disk_nontrivial_slack'])>disk
    near_disk_slack=disk-arb('.06')*(1/r3-arb('.5'))-simplex
    assert near_disk_slack>0
    # Norms of the six Dirichlet Green kernels, each on a pi/3 angular interval.
    max_green=1/(2*r3)
    positive_mass=2/r3-1
    return dict(single_sign_Hausdorff_radius=eps,whole_simplex_Hausdorff_radius='0.000000045',
                disk_Hausdorff_radius='.0046845',near_disk_join_margin=str(near_disk_slack),
                w_interval=j['w_interval'],negative_triangle_leaves=j['negative_triangle_leaves'],
                negative_visited=j['negative_visited'],positive_triangle_leaves=j['positive_triangle_leaves'],
                positive_visited=j['positive_visited'],normalized_hexagon_contacts=True,
                single_sign_curvature_L1_or_primitive_oscillation_limit=str(arb(eps)/max_green),
                simplex_curvature_L1_or_primitive_oscillation_limit=str(simplex/max_green),
                single_sign_curvature_Linfty_limit=str(arb(eps)/positive_mass),
                simplex_curvature_Linfty_limit=str(simplex/positive_mass),
                host_area=j['true_area'],base_sha256=sha256(Path(base.__file__).read_bytes()).hexdigest(),
                specialization_sha256=sha256(new.encode()).hexdigest(),universal_placement_proved=False)

if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('archive',type=Path)
    cli.add_argument('output',type=Path)
    opt=cli.parse_args()
    path=opt.archive.resolve()
    sys.path.insert(0,str(path))
    obj=run(path)
    opt.output.write_text(json.dumps(obj,indent=2)+'\n')
    print(obj)
