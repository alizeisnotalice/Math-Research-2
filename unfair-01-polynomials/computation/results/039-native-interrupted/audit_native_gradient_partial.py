#!/usr/bin/env python3
"""Audit only the recorded progress prefix of the interrupted native GPU run."""
from __future__ import annotations
import argparse,base64,copy,hashlib,json,os,re,sys,time
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(path):
    return [json.loads(x) for x in path.read_text().splitlines() if x.startswith('{')]
def payload_wl(payload):
    m=re.search(r"WL64=('[^']*')",payload);assert m
    return base64.b64decode(eval(m.group(1))).decode()
def audit(root,events,payload):
    wl=(root/'production.wl').read_text();assert payload_wl(payload)==wl
    assert events[0]['event']=='controller_started' and events[0]['kernel_pid']>0
    start=next(x for x in events if x['event']=='started')
    assert start=={'event':'started','pid':start['pid'],'mode':'production','m':60,'n':100,'dft_length':256,'starts':8192,'steps_requested':12000,'seed':2026111617,'scope':'finite nonconvex optimization, not exhaustive search/proof/certificate'}
    rows=[x for x in events if x['event']=='progress']; assert rows
    assert [x['step'] for x in rows]==list(range(100,100*(len(rows)+1),100))
    assert all(x['score_count']==8192 and x['minimum_loss']>=0 and x['seconds']>0 and len(x['best_starts'])==4 for x in rows)
    assert all(rows[i]['seconds']<rows[i+1]['seconds'] for i in range(len(rows)-1))
    assert all(x['event'] not in {'final','complete','controller_exit'} for x in events)
    return {'steps_recorded_through':rows[-1]['step'],'progress_rows':len(rows),'minimum_recorded_loss':min(x['minimum_loss'] for x in rows),'last_recorded_loss':rows[-1]['minimum_loss'],'last_recorded_seconds':rows[-1]['seconds'],'no_terminal_coefficient_rows':True,'no_worker_exit_receipt':True}
def reject(label,f):
    try:f()
    except (AssertionError,KeyError,TypeError,IndexError):return label
    raise AssertionError('accepted mutation '+label)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();assert os.path.ismount('/Volumes/data');beg=time.monotonic()
    events=parse(a.root/'production.stdout');payload=(a.root/'production-payload.py').read_text();main=audit(a.root,events,payload)
    bad_step=copy.deepcopy(events);next(x for x in bad_step if x['event']=='progress')['step']=99
    bad_count=copy.deepcopy(events);next(x for x in bad_count if x['event']=='progress')['score_count']=3
    bad_final=copy.deepcopy(events);bad_final.append({'event':'complete'})
    bad_payload=payload.replace('WL64','WL65',1)
    neg=[reject('nonhundred_progress_step',lambda:audit(a.root,bad_step,payload)),reject('wrong_parallel_score_count',lambda:audit(a.root,bad_count,payload)),reject('forged_completion',lambda:audit(a.root,bad_final,payload)),reject('payload_binding',lambda:audit(a.root,events,bad_payload))]
    verify=json.loads((a.root/'verify_prefix_gradient.json').read_text())
    assert verify['status']=='passed_float64_dft_correlation_gradient_check' and verify['max_finite_difference_error']<2e-7 and verify['fixed_logit_gradients_zero'] is True
    out={'status':'passed_interrupted_native_gradient_progress_audit','scope':'The rebooted production run supplies progress evidence through the listed step only. It has no terminal A/B/C rows, final event, exit receipt, exact factorization, or exclusion conclusion.','progress':main,'formula_regression':{'status':verify['status'],'max_finite_difference_error':verify['max_finite_difference_error'],'fixed_logit_gradients_zero':True},'negative_tests_rejected':neg,'source_sha256':{'production_wl':sha(a.root/'production.wl'),'payload':sha(a.root/'production-payload.py'),'stdout':sha(a.root/'production.stdout'),'auditor':sha(Path(__file__))},'interpreter':sys.version,'elapsed_seconds':time.monotonic()-beg}
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
