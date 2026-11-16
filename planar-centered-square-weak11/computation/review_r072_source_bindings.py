#!/usr/bin/env python3
"""R72 archive-only byte/provenance checks; not a mathematical rerun."""
import argparse
import datetime
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--packet',type=Path,action='append',required=True);ap.add_argument('--archive-root',type=Path,required=True);ap.add_argument('--report',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();r=a.archive_root;topic=r/'planar-centered-square-weak11';packets=[]
    def check(rel,pin,role):
        p=r/rel;assert p.is_file() and sha(p)==pin['sha256'],p
        if 'bytes' in pin:assert p.stat().st_size==pin['bytes'],p
        return dict(destination=rel,sha256=sha(p),bytes=p.stat().st_size,role=role)
    for p in a.packet:
        d=json.loads(p.read_text());rows=[]
        for pin in d['files']+d.get('existing_dependencies',[]):rows.append(check(pin.get('destination','planar-centered-square-weak11/'+pin.get('path','')),pin,'payload/dependency'))
        for pin in d.get('source_pins',[]):rows.append(check('planar-centered-square-weak11/'+pin['path'].split('/planar-centered-square-weak11/',1)[1],pin,'source binding'))
        for binding in d.get('source_bindings',[]):
            pin=binding['source'];rows.append(check('planar-centered-square-weak11/'+pin['path'],pin,'inherited/producer exact input'))
            assert sha(topic/binding['byte_identical_local_copy'])==pin['sha256']
        for key in ('method_source','borel_global_cap_proof'):
            if key in d:
                pin=d[key].get('source',d[key]);rows.append(check('planar-centered-square-weak11/'+pin['path'],pin,key))
        for pin in d.get('external_audit_results_and_receipts',[]):rows.append(check('planar-centered-square-weak11/'+pin['path'],pin,'independent actual run output/receipt'))
        packets.append(dict(manifest=p.name,manifest_sha256=sha(p),checked_rows=rows))
    proof1=topic/'proofs/r072-uniform-grid-borel-family-and-atomic-discontinuity.md';proof2=topic/'proofs/r072-maximizer-tie-breaking-reduction.md'
    assert sha(proof1)=='7b49ddee130e5deee6bbf17641554fbfe122dae9925f73c865574e7073eddf53'
    assert sha(proof2)=='673bf41ad02a5656e123f28db369e9a9a2fbcca63855cf0831678e16ab6fa981'
    own=[]
    for stem in ('r071-positive33','s25-only','aligned225','contracted-unique','uniform-grid-internal','whole-cell-certificate','borel-n5','borel-n17'):
        out=topic/'computation/results/round072-independent'/f'{stem}.json';receipt=out.with_name(f'{stem}-receipt.json');d=json.loads(out.read_text());rec=json.loads(receipt.read_text());assert rec['exit_status']==0 and rec['owned_pid_absent'] and not rec['stderr'];assert rec['output_sha256']==sha(out)
        cmd=rec['command'];rel=cmd[2].split('/planar-centered-square-weak11/',1)[1];assert sha(topic/rel)==d['source_sha256']==rec['source_sha256']
        def cmdpath(flag):return topic/cmd[cmd.index(flag)+1].split('/planar-centered-square-weak11/',1)[1]
        if 'input_sha256'in d and isinstance(d['input_sha256'],dict):
            for k,h in d['input_sha256'].items():assert sha(cmdpath('--'+k.replace('_','-')))==h
            assert d['independent_closure_and_area_helper_sha256']==sha(topic/'computation/review_r072_actual_halo_y_slabs.py')
        elif 'input_sha256'in d:assert sha(cmdpath('--input'))==d['input_sha256']
        if 'saved_result_sha256'in d:assert sha(cmdpath('--saved'))==d['saved_result_sha256'] and sha(cmdpath('--producer'))==d['producer_sha256']
        if 'summary_sha256'in d:assert sha(cmdpath('--summary'))==d['summary_sha256'] and sha(cmdpath('--certificate'))==d['certificate_sha256'] and sha(cmdpath('--producer'))==d['producer_sha256']
        own.append(dict(path=str(out.relative_to(r)),sha256=sha(out),receipt_sha256=sha(receipt),checker_sha256=d['source_sha256'],recorded_pid=rec['pid'],recorded_exit0_absent=True))
    unique=json.loads((topic/'computation/results/round072-independent/contracted-unique.json').read_text());old=json.loads((topic/'computation/results/round072-independent/r071-positive33.json').read_text())
    assert unique['maximum_count']==1 and unique['all_maximum_blocks'][0]['chosen_S']
    lo=F(unique['all_maximum_blocks'][0]['actual_Delta']['lower']);oldhi=F(next(j for j in old['all_maximum_blocks'] if j['chosen_S'])['actual_Delta']['upper'])
    assert lo>F(3,2) and lo>oldhi
    result=dict(status='PASS_R72_ARCHIVE_ONLY_SOURCE_AND_RECEIPT_BINDINGS',checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=sha(Path(__file__)),packets=packets,independent_actual_runs=own,author_proofs=[sha(proof1),sha(proof2)],report=dict(destination='planar-centered-square-weak11/rounds/072.md',sha256=sha(a.report),bytes=a.report.stat().st_size),exact_saved_area_comparisons=dict(unique_lower_above_3_over_2=True,unique_lower_above_R71_chosen_upper=True),scope='Byte/provenance and exact saved-bound comparisons only, all paths resolved inside explicit archive-root; no original workspace fallback or full mathematical re-execution. Proof reading, all finite geometry audits, n5 full leaf replay and n17 partial label replay retain separate scopes.')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status=result['status'],packets=len(packets),packet_rows=sum(len(p['checked_rows']) for p in packets),own_actual_run_chains=len(own))))


if __name__=='__main__':main()
