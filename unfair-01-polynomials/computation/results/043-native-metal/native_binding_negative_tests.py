#!/usr/bin/env python3
"""Five real tamper regressions for native coverage attestation binding."""
import copy,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; SRC=Path('/Volumes/data/CodexTemp/T2/unfair-01-polynomials/round043-native-task13-first')
WORK=HERE/'native-binding-negative-tests';WORK.mkdir(exist_ok=True)
base=json.loads((SRC/'coverage-ranges-attested.jsonl').read_text())
receipt=json.loads((SRC/'receipt.json').read_text());ledger=(SRC/'ledger.jsonl').read_text()
cases={
 'wrong_inputsha':lambda e,r,l:e.__setitem__('input_sha256','0'*64),
 'wrong_receiptsha':lambda e,r,l:e.__setitem__('receipt_sha256','0'*64),
 'wrong_original_ledger':lambda e,r,l:e.__setitem__('original_ledger_sha256','0'*64),
 'wrong_count':lambda e,r,l:e.__setitem__('count',e['count']+1),
 'false_coverage':lambda e,r,l:e.__setitem__('coverage_credit',False),
}
passed=[]
for name,edit in cases.items():
 d=WORK/name;d.mkdir(exist_ok=True);e=copy.deepcopy(base);r=copy.deepcopy(receipt);l=ledger;edit(e,r,l)
 (d/'coverage.jsonl').write_text(json.dumps(e,separators=(',',':'))+'\n');(d/'receipt.json').write_text(json.dumps(r));(d/'ledger.jsonl').write_text(l)
 cmd=[sys.executable,str(HERE/'build_handoff_manifest.py'),'--layer',f'bad={d}/coverage.jsonl','--out',str(d/'out.json')]
 got=subprocess.run(cmd,text=True,capture_output=True)
 if got.returncode==0:raise SystemExit(name+' accepted')
 passed.append(name)
out={'status':'passed','cases':passed};(WORK/'report.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
