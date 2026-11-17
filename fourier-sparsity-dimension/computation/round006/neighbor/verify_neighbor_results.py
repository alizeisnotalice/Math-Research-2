#!/usr/bin/env python3
"""Compare complete quotient scans with different complements and bent criteria."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from verify_switch_closure import echelon, fwht


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--audit',type=Path,required=True);ap.add_argument('--results',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    folder=a.results;row=json.loads(a.audit.read_text())['rows'][1]
    assert Path(row['source_file']).name=='corrected-results.json' and row['row_index']==1
    C=[(1<<64)-1]+[sum(((x>>j)&1)<<x for x in range(64)) for j in range(6)]
    rel=[g for j,g in enumerate(row['relative_generators']) if j!=2]
    piv=echelon(C+rel);assert len(piv)==14
    def rem(x):
        for b in sorted(piv,reverse=True):
            if x>>b&1:x^=piv[b][0]
        return x
    raw=(folder/'remote-neighbor-execution.json').read_bytes();remote=json.loads(raw);r=remote['result'];assert remote['ssh_exit_code']==0
    for n,h in r['source_hashes'].items():assert hashlib.sha256((Path(__file__).parent/n).read_bytes()).hexdigest()==h
    intervals=sorted((x['summary']['start'],x['summary']['count']) for x in r['shards']);at=0
    for start,count in intervals:assert start==at;at+=count
    assert at==1<<28 and len(r['shards'])==14 and all(x['exit_code']==0 for x in r['shards'])
    rh={rem(int(x['relative_phase_hex'],16)) for x in r['found']};lh=set();summ=[]
    for i in range(2):
        s=json.loads((folder/f'local-derivative-{i}.summary.json').read_text());summ.append(s)
        for line in (folder/f'local-derivative-{i}.found.jsonl').read_text().splitlines():
            h=json.loads(line);value=rem(int(h['relative_phase_hex'],16));assert value==int(h['class_normal_form_hex'],16);lh.add(value)
    assert [(x['start'],x['count']) for x in summ]==[(0,1<<27),(1<<27,1<<27)]
    assert rh==lh=={0,rem(row['relative_generators'][2])}
    slices=[0]
    for g in rel:slices += [x^g for x in slices]
    q=[(((x&7)&row['permutation'][x>>3]).bit_count()&1)^((row['h_bits']>>(x>>3))&1) for x in range(64)]
    for h in r['found']:
        vv=int(h['relative_phase_hex'],16)
        for zz in slices:assert all(abs(x)==8 for x in fwht([1-2*(q[x]^(((zz^vv)>>x)&1)) for x in range(64)]))
    result={'status':'passed','remote_tested':at,'local_tested':sum(x['count'] for x in summ),'independent_criteria':['integer Walsh','balanced all nonzero derivatives'],'same_quotient_classes':[hex(x) for x in sorted(rh)],'only_extension_is_original_V':True,'returned_coset_slices_verified':256,'scope':'one specified W of dimension14 inside RM(3,6), not global A28','remote_receipt_sha256':hashlib.sha256(raw).hexdigest(),'local_summaries':summ,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'passed','classes':result['same_quotient_classes'],'all_candidates_each_machine':at}))

if __name__=='__main__':main()
