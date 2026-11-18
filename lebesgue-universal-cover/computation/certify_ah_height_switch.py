"""Reassign the fold bin to the high root when the low-root left sign becomes inconclusive."""
import argparse,inspect,sys
from pathlib import Path
from certify_ah_height_adjusted_chord import run as parent_run
from certify_ah_root_conditioned_support import root_conditioned
from certify_ah_margin_envelope import replace_line


def run(archive,inp,out):
    raw=inspect.getsource(root_conditioned);raw=replace_line(raw,'if j<34 else','if j<35 else');sc=dict(root_conditioned.__globals__);exec(compile(raw,'restored_high_root','exec'),sc)
    code=inspect.getsource(parent_run);env={**parent_run.__globals__,'root_conditioned':sc['root_conditioned'],'__file__':__file__};exec(compile(code,'restored_high_root_area','exec'),env);env['run'](archive,inp,out)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();sys.path.insert(0,str(opt.archive.resolve()));run(opt.archive.resolve(),opt.input,opt.output)
