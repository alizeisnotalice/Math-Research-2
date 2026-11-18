"""Use a lower left endpoint for the first branch of one common-fan reference bin."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import inspect,json
from pathlib import Path

from certify_ah_breakpoint_fan_chord import geometry,run as original,validate as inherited
from certify_ah_flexible_two_roots import either_root as split
from certify_ah_same_fan_cancellation import validate as conditioned
from certify_ah_margin_envelope import replace_line


def validate(ns,j,left,right,turn,sign,orientation,angular=12,depth=4):
    if orientation!=4:
        return inherited(ns,j,left,right,turn,sign,orientation,angular=angular,depth=depth)
    src=inspect.getsource(split)
    src=replace_line(src,'((F(14),anchor,300,-1),(other,F(46),180,1))',
                     f'((F("{left}"),anchor,300,-1),(other,F(46),180,1))')
    additions=(('22.5','36.5','.00002'),('23','36','.00002'),
               ('23.6','35.3','.00001'),('23.2','36','.00003'))
    assignment=f'    options=[tuple(map(F,t)) for t in {additions!r}]+options\n'
    src=replace_line(src,'    rows=[]',assignment+'    rows=[]')
    env={**split.__globals__,'conditional':conditioned}
    exec(compile(src,'lower_first_root_branch','exec'),env)
    return env['either_root'](ns,j)


def run(archive,inp,out):
    src=inspect.getsource(original)
    env={**original.__globals__,'validate':validate,'__file__':__file__}
    exec(compile(src,'lower_split_fan_driver','exec'),env)
    env['run'](archive,inp,out)
    result=json.loads(out.read_bytes())
    dep='certify_ah_breakpoint_fan_chord.py'
    result['source_hashes'][dep]=sha256((archive/dep).read_bytes()).hexdigest()
    out.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ('archive','input','output'):
        parser.add_argument(k,type=Path)
    a=parser.parse_args()
    run(a.archive.resolve(),a.input,a.output)
