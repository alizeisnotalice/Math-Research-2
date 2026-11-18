"""Complete the AH chord after exact subadditivity cancellation of support coefficients on common fans."""
from fractions import Fraction as F
from pathlib import Path
import argparse,inspect,json,sys
from hashlib import sha256
from certify_ah_flexible_chord import run as base_run
from certify_ah_flexible_two_roots import either_root as split
from certify_ah_same_fan_cancellation import validate as conditioned,support_bound
from certify_ah_local_physical_support import validate as ordinary
from certify_ah_taller_floor_chart import geometry


def validate(ns,j,left,right,turn,sign,orientation,angular=12,depth=4):
    if orientation==2:
        src=inspect.getsource(split)
        marker='    rows=[]'
        addition="    options=[(F('22.5'),F('36.5'),F('.00002')),(F('23'),F('36'),F('.00002')),(F('23.6'),F('35.3'),F('.00001')),(F('23.2'),F(36),F('.00003'))]+options\n"
        src=src.replace(marker,addition+marker,1);env={**split.__globals__,'conditional':conditioned}
        exec(compile(src,'same_fan_conditional_brackets','exec'),env);return env['either_root'](ns,j)
    env={**ordinary.__globals__,'support_bound':support_bound}
    exec(compile(inspect.getsource(ordinary),'same_fan_single_brackets','exec'),env)
    return env['validate'](ns,j,left,right,turn,sign,orientation,angular=angular,depth=depth)


def run(archive,inp,out):
    src=inspect.getsource(base_run);env={**base_run.__globals__,'geometry':geometry,'validate':validate,'__file__':__file__}
    exec(compile(src,'global_fan_chord','exec'),env);env['run'](archive,inp,out)
    doc=json.loads(out.read_bytes())
    for name in ('certify_ah_same_fan_cancellation.py','certify_ah_taller_floor_chart.py'):
        doc['source_hashes'][name]=sha256((archive/name).read_bytes()).hexdigest()
    out.write_text(json.dumps(doc,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path)
    opt=ap.parse_args();run(opt.archive.resolve(),opt.input,opt.output)
