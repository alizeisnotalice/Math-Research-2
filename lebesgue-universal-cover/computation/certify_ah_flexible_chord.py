"""Complete an AH chord using fine physical budgets and, when needed, actual contact-class brackets."""
from pathlib import Path
import argparse,inspect,sys,json
from hashlib import sha256
from flint import arb
from certify_ah_chain_chord import run as original
from certify_ah_fine_floor_chart import geometry
from certify_ah_local_physical_support import validate as direct
from certify_ah_flexible_two_roots import either_root


def validate(ns,j,left,right,turn,sign,orientation,angular=12,depth=4):
    if orientation!=2:return direct(ns,j,left,right,turn,sign,orientation,angular=angular,depth=depth)
    return either_root(ns,j)



def run(archive,inp,out):
    source=inspect.getsource(original);env={**original.__globals__,'geometry':geometry,'validate':validate,'__file__':__file__}
    exec(compile(source,'actual_contact_families','exec'),env);env['run'](archive,inp,out)
    data=json.loads(out.read_bytes())
    for name in ('certify_ah_two_hit_budgets.py','certify_ah_pinned_normals.py','certify_ah_root_pose_choices.py','certify_ah_physical_floor.py','certify_ah_fine_floor_chart.py','certify_ah_fine_contact_chart.py','certify_ah_local_physical_support.py','certify_ah_affine_root_class.py','certify_ah_flexible_two_roots.py'):
        data['source_hashes'][name]=sha256((archive/name).read_bytes()).hexdigest()
    out.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path)
    opt=ap.parse_args();sys.path.insert(0,str(opt.archive.resolve()));run(opt.archive.resolve(),opt.input,opt.output)
