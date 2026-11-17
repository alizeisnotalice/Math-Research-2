#!/usr/bin/env python3
"""CPU-only adversarial tests for the d14 resume preflight."""
from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import shutil
import tempfile
from pathlib import Path

import torch

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('resume_d14', HERE / 'resume_gpu_factor_search_prefix.py')
assert SPEC and SPEC.loader
resume = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(resume)


def args(run: Path, state: Path) -> argparse.Namespace:
    return argparse.Namespace(m=8, n=10, starts=16, steps=12, fractional_index=2, margin=.1,
                              restart_every=4, seed=77, out=run.parent / 'out', resume_run=run,
                              resume_state=state, validate_resume_only=True)


def expect_reject(label: str, a: argparse.Namespace, old, state: Path) -> str:
    try:
        resume.load_and_validate_resume(a, old, state)
    except ValueError:
        return label
    raise AssertionError(f'accepted malformed resume fixture: {label}')


def main() -> None:
    old = resume.load_original()
    with tempfile.TemporaryDirectory(dir='/Volumes/data/CodexTemp/T2/unfair-01-polynomials/agents/resume_d14_39') as raw:
        root = Path(raw); run = root / 'run'; run.mkdir(); state_path = run / 'state.pt'; a = args(run, state_path)
        code_sha = resume.sha(Path(old.prefixes.__file__))
        free = torch.zeros((16, 12), dtype=torch.float32); theta = torch.zeros((16, 1), dtype=torch.float32)
        opt = torch.optim.Adam([torch.nn.Parameter(free.clone()), torch.nn.Parameter(theta.clone())], lr=.035)
        loss = (opt.param_groups[0]['params'][0]**2).sum() + (opt.param_groups[0]['params'][1]**2).sum(); loss.backward(); opt.step()
        codes = old.survivor_codes(2); digest = __import__('hashlib').sha256(json.dumps(codes, separators=(',', ':')).encode()).hexdigest()
        p = {'m':8,'n':10,'starts':16,'steps':12,'fractional_index':2,'margin':.1,'restart_every':4,'seed':77}
        config = {**p, 'script_sha256':resume.ORIGINAL_SHA, 'prefix_module_sha256':code_sha, 'prefix_codes_sha256':digest}
        (run/'config.json').write_text(json.dumps(config))
        (run/'checkpoint.json').write_text(json.dumps({'history':[{'step':0},{'step':4},{'step':8},{'step':9}], 'cases':[]}))
        state = {'step':8, 'free':free, 'theta':theta, 'optimizer':old.cpu_tree(opt.state_dict()),
                 'prefix_epoch':2, 'prefix_offset_mod_count':2 % len(codes), 'prefix_codes_sha256':digest, 'parameters':p}
        torch.save(state, state_path)
        for restart_step in (4, 8):
            old_epoch = restart_step // 4 - 1
            (run/f'restart_{restart_step:07d}.json').write_text(json.dumps({
                'old_prefix_epoch':old_epoch, 'old_prefix_offset_mod_count':old_epoch % len(codes),
                'next_prefix_epoch':old_epoch+1, 'next_prefix_offset_mod_count':(old_epoch+1) % len(codes),
                'history':[], 'cases':[]}))
        positive = resume.load_and_validate_resume(a, old, state_path)
        assert positive['step']==8 and [x['step'] for x in positive['discarded_history']]==[9]
        provenance_out = root / 'provenance-out'; provenance_out.mkdir()
        resume.write_resume_provenance(a, positive, provenance_out)
        provenance = json.loads((provenance_out/'resume-provenance.json').read_text())
        assert provenance['next_step']==9 and (provenance_out/'restart_0000004.json').read_bytes()==(run/'restart_0000004.json').read_bytes()
        bad_shape = copy.deepcopy(state); bad_shape['free']=torch.zeros((15,12),dtype=torch.float32); p1=root/'bad-shape.pt';torch.save(bad_shape,p1)
        bad_config=json.loads((run/'config.json').read_text()); bad_config['script_sha256']='0'*64; (run/'config.json').write_text(json.dumps(bad_config)); rejected_hash=expect_reject('config_source_hash',a,old,state_path)
        (run/'config.json').write_text(json.dumps(config))
        rejected_shape=expect_reject('free_shape',a,old,p1)
        bad_params=copy.deepcopy(state);bad_params['parameters']['m']=9;p2=root/'bad-parameters.pt';torch.save(bad_params,p2)
        rejected_params=expect_reject('state_parameters',a,old,p2)
        bad_range=copy.deepcopy(state);bad_range['free'][0,0]=2.;p3=root/'bad-range.pt';torch.save(bad_range,p3)
        rejected_range=expect_reject('free_out_of_range',a,old,p3)
        bad_adam=copy.deepcopy(state);bad_adam['optimizer']['param_groups'][0]['lr']=.1;p4=root/'bad-adam.pt';torch.save(bad_adam,p4)
        rejected_adam=expect_reject('adam_hyperparameter',a,old,p4)
        (run/'restart_0000008.json').unlink()
        rejected_restart=expect_reject('missing_restart_record',a,old,state_path)
        output={'status':'passed_cpu_resume_preflight_negative_tests','positive_step':positive['step'],
                'discarded_history_steps':[x['step'] for x in positive['discarded_history']],
                'negative_tests_rejected':[rejected_hash,rejected_shape,rejected_params,rejected_range,rejected_adam,rejected_restart],
                'provenance_copy_checked':True, 'torch':torch.__version__}
        (HERE/'test-resume-result.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,sort_keys=True))


if __name__ == '__main__':
    main()
