"""Exact decimal sample and window bookkeeping audit (not search exclusion)."""
import argparse,hashlib,importlib.util,json,sys
from fractions import Fraction
from pathlib import Path

def main():
    if not __debug__:raise RuntimeError('assertion checks require no -O')
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--verifier',type=Path,required=True);p.add_argument('--worker',type=Path,required=True);p.add_argument('--prefix-plan',type=Path);p.add_argument('--base-source',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.run
    sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
    c=json.loads((r/'config.json').read_text());result=json.loads((r/'result.json').read_text());release=json.loads((r/'release.json').read_text());x=json.loads((r/'window-candidates.json').read_text())
    plan=a.prefix_plan or Path(c['prefix_plan_path']);base=a.base_source or Path(c['frozen_source_path']);pool_json=json.loads(plan.read_text());pool=set(pool_json['survivor_codes']);assert len(pool)==46011
    assert sha(plan)==c['prefix_plan_sha256'] and sha(base)==c['frozen_source_sha256'] and sha(a.worker)==c['script_sha256']
    assert result['complete'] and result['last_completed_step']==36000 and c['seed']==x['seed']==2026111704
    assert release['terminal_status']=='complete' and release['cleanup_called'] is True and release['timed_out'] is False and release['last_completed_step']==36000
    assert x['starts']==8192 and x['restart_every']==6000 and x['prefix_codes_sha256']==c['prefix_codes_sha256']
    assert [w['window_epoch'] for w in x['windows']]==list(range(6))
    spec=importlib.util.spec_from_file_location('fraction_verifier',a.verifier);v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
    output=[];seen=set()
    for epoch,w in enumerate(x['windows']):
        assert w['window_offset_mod_count']==epoch*8192%len(pool) and len(w['candidates'])==16
        assert (epoch+1)*6000-100<=w['closing_step']<=(epoch+1)*6000
        window={pool_json['survivor_codes'][(epoch*8192+i)%len(pool)] for i in range(8192)};seen|=window
        for case in w['candidates']:
            res,delta,comp=v.validate_case(case,64,104,15,Fraction(1,1000),window)
            assert res>0
            output.append({'epoch':epoch,'exact_decimal_residual':str(res),'float_residual':float(res),'complement_error':str(comp)})
    assert seen==pool and len(output)==96
    report={'status':'passed_96_saved_decimal_samples_and_six_window_coverage','scope':'No saved sample is an exact counterexample; no finite exclusion or trajectory-wide verification. Prefix pool completeness inherited from certified input.','samples':len(output),'prefixes':len(seen),'min_residual':min(z['float_residual'] for z in output),'residuals':output,'bindings':{n:sha(r/n) for n in ['config.json','result.json','release.json','window-candidates.json']},'verifier_sha256':sha(a.verifier),'script_sha256':sha(Path(__file__)),'python':sys.version}
    a.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['status','samples','prefixes','min_residual']}))

if __name__=='__main__':main()
