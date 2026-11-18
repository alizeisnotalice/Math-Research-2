"""Reuse the eleven EH facets and directly recheck the narrow flexible AH root classes."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import inspect
import json
from pathlib import Path
import sys

from flint import arb,ctx


def read(path):
    return json.loads(path.read_bytes())


def verify_weak(ah_archive,ah,normal,offset,cost):
    import certify_ah_flexible_chord as flexible

    cfg=ah['parameters']
    source=ah_archive/'certify_ah_flexible_chord.py'
    assert sha256(source.read_bytes()).hexdigest()==ah['script_sha256']
    extra=[tuple(map(F,row)) for row in cfg['lines']]
    assert extra==[(F('.427'),F('.006'))]
    assert 'root_conditioned_classes' not in ah
    closed=ah['bins']+ah['common_classes']
    weak={r['reference_bin'] for r in ah['common_classes'] if not arb(r['placement_margin_lower'])>cost}
    ordinary=[arb(r['placement_margin_lower'])-cost for r in closed if arb(r['placement_margin_lower'])>cost]
    assert len(ordinary)+len(weak)==len(closed) and weak
    far,negative,eps=map(F,cfg['physical_parameters'])
    old_precision=ctx.prec
    try:
        ns,full,_=flexible.geometry(ah_archive,extra,far=far,negative=negative,eps=eps)
        assert sha256(full.encode()).hexdigest()==ah['geometry_sha256']
        ns['rows']=list(ns['rows'])+[('EH-deep-chord',normal,offset)]
        out=[]
        for low,high,left,right,turn,sign,orientation in cfg['classes']:
            for j in sorted(weak.intersection(range(low,high))):
                row=flexible.validate(ns,j,F(left),F(right),turn,sign,orientation,angular=12,depth=4)
                assert arb(row['root_margin_lower'])>0 and arb(row['placement_margin_lower'])>0,row
                out.append(row)
        assert sorted(r['reference_bin'] for r in out)==sorted(weak)
    finally:
        ctx.prec=old_precision
    limits=[v.lower() for v in ordinary]+[arb(r['placement_margin_lower']).lower() for r in out]
    return min(limits),out


def run(archive,ah_archive,ah_file,parameters,weights,out):
    from certify_eh_upper_chain import run as original
    from certify_ah_margin_envelope import replace_line

    source=inspect.getsource(original)
    source=replace_line(source,'        from certify_eh_pinned_weak_classes import verify',
                         '        verify=verify_flexible_weak')
    closure=lambda archive,ah,n,c,cost: verify_weak(ah_archive,ah,n,c,cost)
    env={**original.__globals__,'verify_flexible_weak':closure}
    exec(compile(source,'EH_chain_with_flexible_AH_transfer','exec'),env)
    data=env['run'](archive,ah_file,parameters,read(weights),False,None)
    data['weights_sha256']=sha256(weights.read_bytes()).hexdigest()
    data['script_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    data['specialized_EH_verifier_sha256']=sha256(source.encode()).hexdigest()
    out.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    for key in ('archive','ah_archive','ah','parameters','weights','out'):
        cli.add_argument(key,type=Path)
    opt=cli.parse_args()
    archive,ah_archive=opt.archive.resolve(),opt.ah_archive.resolve()
    sys.path.insert(0,str(archive))
    sys.path.insert(0,str(ah_archive))
    run(archive,ah_archive,opt.ah,opt.parameters,opt.weights,opt.out)
