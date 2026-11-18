"""Verify the full single-chord host using chart-deficit derivatives in physical coordinates."""
import argparse,inspect,json,sys
from hashlib import sha256
from pathlib import Path
from certify_ah_chain_chord import run as previous_run
from certify_ah_centered_deficit import geometry


def run(archive,inp,out):
    src=inspect.getsource(previous_run);env={**previous_run.__globals__,'geometry':geometry,'__file__':__file__};exec(compile(src,'direct_chart_complete','exec'),env);env['run'](archive,inp,out)
    result=json.loads(out.read_bytes());result['source_hashes']['certify_ah_centered_deficit.py']=sha256(Path(inspect.getfile(geometry)).read_bytes()).hexdigest();out.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('archive',type=Path);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path);opt=ap.parse_args();sys.path.insert(0,str(opt.archive.resolve()));run(opt.archive.resolve(),opt.input,opt.output)
